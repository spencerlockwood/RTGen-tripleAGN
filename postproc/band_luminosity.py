"""
band_luminosity.py
==================

Luminosity in a chosen band for each member of a triple-AGN run: the host
galaxy's starlight (from the RADMC-3D images) and an estimate for the AGN
itself (from its bolometric luminosity, a bolometric correction and the dust
column along the camera's line of sight through the RT grid).

    IMPORTANT: the AGN are NOT light sources in the radiative transfer. The
    images contain stars (+ dust emission heated by stars) only. The AGN
    numbers below are an add-on estimate, not RT output. See
    docs/band_luminosities.md for the method and its caveats.

Inputs, all inside RUN_DIR:
    image_cont.out                                   dusty continuum image
    dust_obscuration_analysis/image_cont_nodust.out  dust-free image (optional)
    sph_interpolation/time_and_z.txt                 redshift
    sph_interpolation/bh_info.txt                    BH offsets (or, for older runs,
                                                     parsed from job_sph_interpolation.out)
    amr_grid.inp, dust_density.inp, dustopac.inp, dustkappa_*.inp, input/mass_gas.DAT
                                                     (only for the AGN line-of-sight column)

Band choices (exactly one):
    --band-rest LO HI     top-hat in REST-frame wavelength [um]
    --band-obs  LO HI     top-hat in OBSERVED-frame wavelength [um] (shifted to rest by 1/(1+z))
    --filter FILE         filter curve, 2 columns: wavelength [Angstrom], throughput
                          (format of the files in filters/); observed frame unless
                          --filter-frame rest
    --xray EMIN EMAX      X-ray band [keV], AGN only (the RT has no X-ray emission)

AGN bolometric luminosities (first one found is used):
    --lbol H N1 N2        erg/s
    bh_lbol_erg_s         in RUN_DIR/input/input_params.py (filled by scripts/make_run.sh)
    --catalog CSV         picks CSV from selection/find_dust_analogues.py; row matched by host BH id

Examples:
    python postproc/band_luminosity.py RUN --band-rest 0.5 0.6
    python postproc/band_luminosity.py RUN --filter filters/JWST/<file>.dat --catalog dust_analogues.csv
    python postproc/band_luminosity.py RUN --band-rest 0.48 0.54 --agn-kbol 9.26
    python postproc/band_luminosity.py RUN --xray 2 10 --lbol 8.3e43 8.9e40 1.7e40
    python postproc/band_luminosity.py RUN --band-rest 0.5 0.6 --json out.json
"""

import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_dust_obscuration import read_radmc3d_image, read_snapshot_redshift  # noqa: E402

C_UM = 2.99792458e14          # speed of light [um/s]
KPC_CM = 3.0857e21
MPC_CM = 3.0857e24
MSUN_G = 1.989e33
LSUN = 3.828e33
M_P = 1.6726e-24
X_H = 0.76                    # hydrogen mass fraction for N_H
# Romulus25 cosmology (Planck 2014), used only for D_L / apparent magnitudes
H0, OM = 67.77, 0.3086

# AGN bolometric corrections L_bol / (lambda L_lambda), Shen et al. (2008),
# from the Richards et al. (2006) mean quasar SED. Used automatically when the
# band's rest-frame pivot wavelength is within 15% of one of these.
OPT_BC = {0.1350: 3.81, 0.3000: 5.15, 0.5100: 9.26}
# Duras et al. (2020) hard X-ray (2-10 keV) correction, general AGN sample:
# K_X = a [1 + (log10(L_bol/Lsun) / b)^c]
DURAS = dict(a=10.96, b=11.93, c=17.79)
GAMMA_X = 1.8                 # photon index used to convert 2-10 keV to other X-ray bands

ROLES = ['Host', 'N1', 'N2']
trapz = getattr(np, 'trapezoid', None) or np.trapz


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def lum_distance_cm(z):
    zz = np.linspace(0, z, 2001)
    Ez = np.sqrt(OM * (1 + zz) ** 3 + 1 - OM)
    dc = (C_UM / 1e9) / H0 * trapz(1 / Ez, zz)          # Mpc (c in km/s)
    return (1 + z) * dc * MPC_CM


def read_bh_offsets(run):
    """BH offsets from the host [kpc, simulation axes] for Host, N1, N2."""
    f = os.path.join(run, 'sph_interpolation', 'bh_info.txt')
    if os.path.exists(f):
        rows = [l.split() for l in open(f) if l.strip() and not l.startswith('#')]
        return {r[0]: (int(r[1]), np.array(r[2:5], float)) for r in rows}
    log = os.path.join(run, 'sph_interpolation', 'job_sph_interpolation.out')
    pat = re.compile(r'(Host|N1|N2)\s+BH (\d+): offset from grid centre \(kpc\) = \[([^\]]+)\]')
    out = {}
    for line in open(log):
        m = pat.search(line)
        if m:
            out[m.group(1)] = (int(m.group(2)), np.array(m.group(3).split(), float))
    if len(out) != 3:
        raise RuntimeError(f'No bh_info.txt and could not parse 3 BH offsets from {log}')
    return out


def read_params_lbol(run):
    p = os.path.join(run, 'input', 'input_params.py')
    if not os.path.exists(p):
        return None
    m = re.search(r'^bh_lbol_erg_s\s*=\s*(.+)$', open(p).read(), re.M)
    if not m or m.group(1).strip().startswith(('None', '@')):
        return None
    return [float(v) for v in re.findall(r'[-+0-9.eE]+', m.group(1))][:3]


def catalog_lbol(csv, host_iord, iords):
    d = pd.read_csv(csv)
    row = d[d['BH_id'].astype(np.int64) == host_iord]
    if row.empty:
        raise RuntimeError(f'host BH {host_iord} not found in {csv}')
    row = row.iloc[0]
    ids = [int(row[c]) for c in ('BH_id', 'Neighbour1BH_id', 'Neighbour2BH_id')]
    if ids != list(iords):
        print(f'  WARNING: catalogue member ids {ids} differ from run ids {list(iords)}')
    return [float(row[c]) for c in ('BH_Lbol', 'Neighbour1BH_Lbol', 'Neighbour2BH_Lbol')]


def band_weights(args, z):
    """Return (lam_rest_fine [um], T_fine, label). T is the throughput on a
    fine rest-frame wavelength grid."""
    if args.band_rest:
        lo, hi = args.band_rest; frame = 'rest'
        lab = f'rest-frame {lo:g}-{hi:g} um (top-hat)'
        lam = np.geomspace(lo, hi, 2000); T = np.ones_like(lam)
    elif args.band_obs:
        lo, hi = args.band_obs; frame = 'obs'
        lab = f'observed-frame {lo:g}-{hi:g} um (top-hat) = rest {lo/(1+z):.4g}-{hi/(1+z):.4g} um'
        lam = np.geomspace(lo, hi, 2000) / (1 + z); T = np.ones_like(lam)
    else:
        f = np.genfromtxt(args.filter)
        lam_f, T_f = f[:, 0] / 1e4, f[:, 1]
        frame = args.filter_frame
        good = T_f > 0
        lam_f, T_f = lam_f[good], T_f[good]
        if frame == 'obs':
            lam_f = lam_f / (1 + z)
        lam = np.geomspace(lam_f.min(), lam_f.max(), 3000)
        T = np.interp(lam, lam_f, T_f)
        lab = f'{os.path.basename(args.filter)} ({frame} frame; rest {lam.min():.4g}-{lam.max():.4g} um)'
    return lam, T, lab, frame


def integrate_band(lam_cam, Lnu_cam, lam, T):
    """Band quantities of a spectrum L_nu(lam_cam) [erg/s/Hz] sampled on the
    camera grid, for throughput T(lam) on a fine rest-frame grid.

    Returns
      L_band   = int L_nu T dnu                 [erg/s]  energy in the band
      Lnu_mean = int L_nu T dnu/nu / int T dnu/nu   [erg/s/Hz] (AB / photon-counting mean)
    """
    good = Lnu_cam > 0
    logL = np.interp(np.log(lam), np.log(lam_cam[good]), np.log(Lnu_cam[good]))
    Lnu = np.exp(logL)
    nu = C_UM / lam
    # lam ascending -> nu descending; integrate in ln(lam): dnu = nu dln(lam) (abs)
    lnl = np.log(lam)
    L_band = trapz(Lnu * T * nu, lnl)
    Lnu_mean = trapz(Lnu * T, lnl) / trapz(T, lnl)
    return L_band, Lnu_mean


def pivot_um(lam, T):
    return np.sqrt(trapz(T * lam, lam) / trapz(T / lam, lam))


def agn_kbol_optical(lam_piv, override):
    if override is not None:
        return override, 'user --agn-kbol'
    for l0, k in OPT_BC.items():
        if abs(lam_piv / l0 - 1) < 0.15:
            return k, f'Shen+08 BC at {l0*1e4:.0f} A'
    return None, None


def kx_duras(lbol):
    l = np.log10(lbol / LSUN)
    return DURAS['a'] * (1 + (l / DURAS['b']) ** DURAS['c'])


def xray_band_fraction(e1, e2, gamma=GAMMA_X):
    """L(e1-e2) / L(2-10 keV) for a power law with photon index gamma."""
    s = 2 - gamma
    f = (lambda a, b: np.log(b / a)) if abs(s) < 1e-9 else (lambda a, b: (b ** s - a ** s) / s)
    return f(e1, e2) / f(2.0, 10.0)


# --------------------------------------------------------------------------
# line-of-sight dust / gas column from each BH to the observer
# --------------------------------------------------------------------------
def los_columns(run, offsets_kpc, lam_eval):
    """Dust optical depth tau(lam_eval) and N_H from each BH to the observer.

    Geometry (verified against the images): sph_interpolation.py writes the
    grid with simulation (x, z, y) -> RADMC (x, y, z), and the image command
    uses incl 90, so the camera looks along RADMC +y = simulation +z from the
    -y side. The column is the sum over cells with the BH's (ix, iz) and
    RADMC y below the BH (half of the BH's own cell included).
    """
    g = open(os.path.join(run, 'amr_grid.inp')).read().split()
    nx, ny, nz = map(int, g[7:10])
    walls = np.array(g[10:10 + nx + ny + nz + 3], float)
    xw, yw, zw = walls[:nx + 1], walls[nx + 1:nx + ny + 2], walls[nx + ny + 2:]
    dyc = np.diff(yw)

    def cell(w, v):
        return int(np.clip(np.searchsorted(w, v) - 1, 0, len(w) - 2))

    # species and opacities from dustopac.inp
    op = [l.split()[0] for l in open(os.path.join(run, 'dustopac.inp'))]
    nspec = int(op[1])
    names = [op[i] for i in range(len(op)) if op[i - 1] == '0' and op[i - 2] == '1' and i > 3][:nspec]
    kext = []
    for nm in names:
        L = [l.split() for l in open(os.path.join(run, f'dustkappa_{nm}.inp'))
             if l.strip() and not l.lstrip().startswith('#')]
        a = np.array([r[:3] for r in L[2:]], float)
        k = a[:, 1] + a[:, 2]                    # absorption + scattering = extinction [cm^2/g dust]
        kext.append(np.exp(np.interp(np.log(lam_eval), np.log(a[:, 0]), np.log(k))))

    rho = pd.read_csv(os.path.join(run, 'dust_density.inp'), skiprows=3, header=None,
                      sep=r'\s+', engine='c').values.ravel().reshape(nspec, nz, ny, nx)
    gas_path = os.path.join(run, 'input', 'mass_gas.DAT')
    mgas = np.fromfile(gas_path, dtype=np.float64).reshape(nz, ny, nx) if os.path.exists(gas_path) else None
    vol = np.diff(xw)[0] * dyc[0] * np.diff(zw)[0]

    out = {}
    for role, off in offsets_kpc.items():
        X, Y, Z = off[0] * KPC_CM, off[2] * KPC_CM, off[1] * KPC_CM     # sim (x,y,z) -> RADMC (x, z, y)
        ix, iy, iz = cell(xw, X), cell(yw, Y), cell(zw, Z)
        frac = np.ones(iy + 1); frac[-1] = (Y - yw[iy]) / dyc[iy]
        sig = [np.sum(rho[s, iz, :iy + 1, ix] * dyc[:iy + 1] * frac) for s in range(nspec)]   # g/cm^2
        tau = sum(s * k for s, k in zip(sig, kext))
        NH = (np.sum(mgas[iz, :iy + 1, ix] * frac) * MSUN_G / vol * dyc[0] * X_H / M_P
              if mgas is not None else np.nan)
        out[role] = dict(tau=tau, Sigma_dust=float(sum(sig)), N_H=float(NH))
    return out


# --------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run')
    b = ap.add_mutually_exclusive_group(required=True)
    b.add_argument('--band-rest', type=float, nargs=2, metavar=('LO_UM', 'HI_UM'))
    b.add_argument('--band-obs', type=float, nargs=2, metavar=('LO_UM', 'HI_UM'))
    b.add_argument('--filter')
    b.add_argument('--xray', type=float, nargs=2, metavar=('EMIN_KEV', 'EMAX_KEV'))
    ap.add_argument('--filter-frame', choices=['obs', 'rest'], default='obs')
    ap.add_argument('--aperture-kpc', type=float, default=1.0,
                    help='radius of the circle around each BH (pixels go to the nearest BH). Default 1 kpc')
    ap.add_argument('--lbol', type=float, nargs=3, metavar=('HOST', 'N1', 'N2'))
    ap.add_argument('--catalog')
    ap.add_argument('--agn-kbol', type=float, default=None,
                    help='L_bol / (nu L_nu) of the AGN in this band; overrides the built-in corrections')
    ap.add_argument('--no-agn-dust', action='store_true', help='skip the line-of-sight dust column (faster)')
    ap.add_argument('--json', help='also write the results to this JSON file')
    a = ap.parse_args(argv)

    run = os.path.abspath(a.run)
    z = read_snapshot_redshift(os.path.join(run, 'sph_interpolation', 'time_and_z.txt'))
    DL = lum_distance_cm(z)
    bh = read_bh_offsets(run)
    iords = [bh[r][0] for r in ROLES]
    offs = {r: bh[r][1] for r in ROLES}

    lbol, src = a.lbol, '--lbol'
    if lbol is None:
        lbol, src = read_params_lbol(run), 'input_params.py'
    if lbol is None and a.catalog:
        lbol, src = catalog_lbol(a.catalog, iords[0], iords), a.catalog
    if lbol is None:
        src = None

    print(f'Run: {run}')
    print(f'z = {z:.4f}   D_L = {DL / MPC_CM:.1f} Mpc')
    for r in ROLES:
        o = offs[r]
        print(f'  {r:4s} BH {bh[r][0]}: offset (x, y, z) = ({o[0]:+.2f}, {o[1]:+.2f}, {o[2]:+.2f}) kpc; '
              f'projected (image x, y) = ({o[0]:+.2f}, {o[1]:+.2f})')
    res = dict(run=run, z=z, D_L_Mpc=DL / MPC_CM, members={r: {} for r in ROLES})

    # ---------------- X-ray: AGN only ----------------
    if a.xray:
        e1, e2 = a.xray
        print(f'\nBand: {e1:g}-{e2:g} keV (AGN only; the RT contains no X-ray emission)')
        if lbol is None:
            sys.exit('Need AGN L_bol: pass --lbol, --catalog, or set bh_lbol_erg_s in input_params.py')
        cols = None if a.no_agn_dust else los_columns(run, offs, np.array([0.55]))
        print(f'AGN L_bol from {src}. K_X: Duras+20 (2-10 keV), Gamma = {GAMMA_X} for other bands.')
        print(f'{"":5s} {"L_bol":>10s} {"K_X(2-10)":>9s} {"L_X,int":>10s} {"N_H(los)":>10s}')
        for r, L in zip(ROLES, lbol):
            kx = kx_duras(L)
            LX = L / kx * xray_band_fraction(e1, e2)
            NH = cols[r]['N_H'] if cols else np.nan
            print(f'{r:5s} {L:10.3e} {kx:9.2f} {LX:10.3e} {NH:10.3e}')
            res['members'][r] = dict(iord=bh[r][0], Lbol=L, K_X_2_10=kx, L_X_intrinsic=LX, N_H_los=NH)
        print('L_X,int is intrinsic (unabsorbed). N_H is the resolved gas column from the BH to the\n'
              'observer through the RT grid (no torus); apply your own absorption model for observed L_X.')
        if a.json:
            json.dump(res, open(a.json, 'w'), indent=1, default=float)
        return res

    # ---------------- images ----------------
    lam, T, lab, frame = band_weights(a, z)
    lam_p = pivot_um(lam, T)
    print(f'\nBand: {lab}; rest-frame pivot {lam_p:.4g} um')

    dusty = read_radmc3d_image(os.path.join(run, 'image_cont.out'))
    fclear = os.path.join(run, 'dust_obscuration_analysis', 'image_cont_nodust.out')
    clear = read_radmc3d_image(fclear) if os.path.exists(fclear) else None
    lam_cam = dusty['lam']
    if lam.min() < lam_cam.min() or lam.max() > lam_cam.max():
        sys.exit(f'Band (rest {lam.min():.3g}-{lam.max():.3g} um) is outside the camera grid')
    dlog = np.median(np.diff(np.log10(lam_cam)))
    if np.log10(lam.max() / lam.min()) < 2 * dlog:
        print(f'  WARNING: band narrower than 2 camera samples ({dlog:.3f} dex each): '
              'result is an interpolation; re-render the image with finer wavelengths for narrow bands')
    if lam.max() > 160.:
        print('  WARNING: band reaches rest > 160 um, where the Starburst99 input is extrapolated '
              '(the dust-free numbers there are not physical)')

    nx, ny = dusty['nx'], dusty['ny']
    px, py = dusty['sizepix_x'] / KPC_CM, dusty['sizepix_y'] / KPC_CM
    xc = (np.arange(nx) + 0.5 - nx / 2) * px          # image x = simulation x  (host at 0)
    yc = (np.arange(ny) + 0.5 - ny / 2) * py          # image y = simulation y
    XX, YY = np.meshgrid(xc, yc)
    d2 = np.stack([(XX - offs[r][0]) ** 2 + (YY - offs[r][1]) ** 2 for r in ROLES])
    nearest = np.argmin(d2, axis=0)
    masks = {r: (nearest == k) & (d2[k] <= a.aperture_kpc ** 2) for k, r in enumerate(ROLES)}
    masks['Total'] = np.ones_like(nearest, bool)
    for i in range(3):
        for j in range(i + 1, 3):
            s = np.hypot(*(offs[ROLES[i]][:2] - offs[ROLES[j]][:2]))
            if s < 2 * a.aperture_kpc:
                print(f'  note: {ROLES[i]}-{ROLES[j]} projected separation {s:.2f} kpc < 2 x aperture; '
                      'shared pixels go to the nearer BH')

    pixA = dusty['sizepix_x'] * dusty['sizepix_y']
    def spectra(im):
        # L_nu [erg/s/Hz] = 4 pi * sum(I_nu) * pixel area  (isotropic-equivalent, this sightline)
        return {r: 4 * np.pi * pixA * im['image'][:, m].sum(axis=1) for r, m in masks.items()}
    Sd = spectra(dusty)
    Sc = spectra(clear) if clear is not None else None

    # ---------------- AGN ----------------
    kbol, ksrc = agn_kbol_optical(lam_p, a.agn_kbol)
    cols = None
    if lbol is not None and kbol is not None and not a.no_agn_dust:
        print('  reading dust_density.inp for the AGN line-of-sight column ...')
        cols = los_columns(run, offs, lam)
    lnl = np.log(lam)
    wT = trapz(T, lnl)   # int T dnu/nu

    print(f'\nStars (+dust emission) from the RT images, aperture r = {a.aperture_kpc:g} kpc around each BH')
    hdr = f'{"":6s} {"L_band dusty":>13s} {"L_band nodust":>13s} {"A_band":>7s} {"nuLnu dusty":>12s} {"M_AB dusty":>10s}'
    if frame == 'obs':
        hdr += f' {"m_AB dusty":>10s}'
    print(hdr + '\n' + ' ' * 7 + '[erg/s]       [erg/s]       [mag]   [erg/s]      [rest-band]' + ('  [observed]' if frame == 'obs' else ''))
    for r in ROLES + ['Total']:
        Lb, Lnu = integrate_band(lam_cam, Sd[r], lam, T)
        rec = dict(L_band_stars_dusty=Lb, Lnu_mean_stars_dusty=Lnu, nuLnu_pivot_stars_dusty=Lnu * C_UM / lam_p)
        Lb0 = Lnu0 = np.nan
        if Sc is not None:
            Lb0, Lnu0 = integrate_band(lam_cam, Sc[r], lam, T)
            rec.update(L_band_stars_nodust=Lb0, A_band_stars=-2.5 * np.log10(Lnu / Lnu0))
        M_AB = -2.5 * np.log10(Lnu / (4 * np.pi * (10 * 3.0857e18) ** 2)) - 48.60
        rec['M_AB_stars_dusty'] = M_AB
        line = (f'{r:6s} {Lb:13.4e} {Lb0:13.4e} {rec.get("A_band_stars", np.nan):7.2f} '
                f'{Lnu * C_UM / lam_p:12.4e} {M_AB:10.2f}')
        if frame == 'obs':
            Fnu_obs = (1 + z) * Lnu / (4 * np.pi * DL ** 2)
            rec['m_AB_stars_dusty'] = -2.5 * np.log10(Fnu_obs) - 48.60
            rec['F_band_stars_dusty_cgs'] = Lb / (4 * np.pi * DL ** 2)
            line += f' {rec["m_AB_stars_dusty"]:10.2f}'
        print(line)
        res['members'].setdefault(r, {}).update(rec)

    if lbol is None:
        print('\nAGN: no L_bol given (use --lbol, --catalog, or bh_lbol_erg_s in input_params.py)')
    elif kbol is None:
        print(f'\nAGN: no built-in bolometric correction for a rest-frame pivot of {lam_p:.3g} um. '
              'Pass --agn-kbol K (= L_bol / nu L_nu in this band).')
    else:
        print(f'\nAGN estimate (NOT from the RT): L_bol from {src}; nuL_nu = L_bol / K, K = {kbol:g} ({ksrc});')
        print('flat nuL_nu across the band; dust extinction (abs+scat) along the camera sightline from the BH')
        print(f'{"":6s} {"L_bol":>10s} {"nuLnu,int":>10s} {"L_band,int":>10s} {"tau_band":>8s} {"L_band,obs":>10s} '
              f'{"Sigma_d":>9s} {"N_H(los)":>10s} {"AGN frac":>8s}')
        for r, L in zip(ROLES, lbol):
            nuLnu = L / kbol
            Lband_int = nuLnu * wT          # = int (nuL_nu) T dln(nu) for flat nuL_nu
            if cols:
                trans = trapz(T * np.exp(-cols[r]['tau']), lnl) / wT
                tau_b = -np.log(max(trans, 1e-300))
            else:
                trans, tau_b = 1.0, 0.0
            Lband_obs = Lband_int * trans
            star = res['members'][r]['L_band_stars_dusty']
            frac = Lband_obs / (Lband_obs + star)
            c = cols[r] if cols else dict(Sigma_dust=np.nan, N_H=np.nan)
            print(f'{r:6s} {L:10.3e} {nuLnu:10.3e} {Lband_int:10.3e} {tau_b:8.2f} {Lband_obs:10.3e} '
                  f'{c["Sigma_dust"]:9.2e} {c["N_H"]:10.3e} {frac:8.3f}')
            res['members'][r].update(Lbol=L, K_band=kbol, nuLnu_agn_intrinsic=nuLnu,
                                     L_band_agn_intrinsic=Lband_int, tau_band_agn=tau_b,
                                     L_band_agn_attenuated=Lband_obs, Sigma_dust_los=c['Sigma_dust'],
                                     N_H_los=c['N_H'], agn_fraction=frac,
                                     L_band_total=Lband_obs + star)
        print('Member total = stars (dusty, in aperture) + AGN (attenuated). AGN frac = AGN / total.')

    if a.json:
        json.dump(res, open(a.json, 'w'), indent=1, default=float)
        print(f'Wrote {a.json}')
    return res


if __name__ == '__main__':
    main()
