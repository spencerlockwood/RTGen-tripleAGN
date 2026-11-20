"""
plot_dust_obscuration.py
========================

Dusty vs. dust-free continuum comparison for the RTGen triple-AGN runs.
Same analysis as cell 3 of RTGen-public_bh1/analysis.ipynb:

  * integrated flux (sum over pixels) of the dusty and dust-free images,
    vs. observed-frame wavelength  lam_obs = lam_rest * (1+z)
  * dust attenuation  A = -2.5 log10(F_dusty / F_dustfree)  [mag]
  * spatial attenuation map at the wavelength closest to observed V (0.55 um)

For each run directory it needs
  image_cont.out                                    (submit_master_pipeline.sh)
  dust_obscuration_analysis/image_cont_nodust.out   (submit_dust_obscuration.sh)
  sph_interpolation/time_and_z.txt                  (rewrite_sim.py)

and writes rtgen_dust_obscuration.png / rtgen_dust_obscuration_map.png into
that run directory. With more than one run it also writes a comparison of the
attenuation curves, dust_attenuation_all_systems.png, in the current directory.

Usage:
    python postproc/plot_dust_obscuration.py RUN_DIR [RUN_DIR ...]
    cd RUNS_DIR && python /path/to/repo/postproc/plot_dust_obscuration.py   # every sys*_snap* run here
"""

import os
import sys
import glob

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.getcwd()   # default place to look for runs and to write the comparison plot
KPC_CM = 3.0857e21
SB99_LAM_MAX = 160.0  # um, rest frame: longest wavelength in the Starburst99 spectra


def read_radmc3d_image(fname):
    """Parse a RADMC-3D ASCII image.out file (pandas for the body: much
    faster than np.loadtxt on these ~200 MB files)."""
    with open(fname, 'r') as f:
        iformat = int(f.readline().split()[0])
        nx, ny = (int(v) for v in f.readline().split())
        nf = int(f.readline().split()[0])
        sizepix_x, sizepix_y = (float(v) for v in f.readline().split())
        lam = np.array([float(f.readline().split()[0]) for _ in range(nf)])
    nhead = 4 + nf
    body = pd.read_csv(fname, skiprows=nhead, header=None, sep=r'\s+',
                       engine='c', skip_blank_lines=True).values.ravel()
    stokes = iformat in (3, 4)
    image = body.reshape(nf, ny, nx, 4)[..., 0] if stokes else body.reshape(nf, ny, nx)
    return {'nx': nx, 'ny': ny, 'nf': nf, 'sizepix_x': sizepix_x,
            'sizepix_y': sizepix_y, 'lam': lam, 'image': image}


def read_snapshot_redshift(fname):
    """Redshift of the snapshot, as written by rewrite_sim.py."""
    with open(fname) as f:
        for line in f:
            if line.strip() and not line.strip().startswith('#'):
                return float(line.split()[2])
    raise ValueError(f'Could not find a data line in {fname}')


def analyse(run):
    name = os.path.basename(os.path.normpath(run))
    im_dusty = read_radmc3d_image(os.path.join(run, 'image_cont.out'))
    im_clear = read_radmc3d_image(os.path.join(run, 'dust_obscuration_analysis', 'image_cont_nodust.out'))
    assert np.allclose(im_dusty['lam'], im_clear['lam']), \
        f'{name}: wavelength grids differ between the dusty and dust-free images'

    z = read_snapshot_redshift(os.path.join(run, 'sph_interpolation', 'time_and_z.txt'))
    lam_rest = im_dusty['lam']
    lam = lam_rest * (1.0 + z)   # observed frame; the flux ratio needs no distance correction

    flux_dusty = im_dusty['image'].sum(axis=(1, 2))
    flux_clear = im_clear['image'].sum(axis=(1, 2))

    # Starburst99 spectra end at 160 um (rest). problem_setup_dust_continuum.py
    # extrapolates them out to 1e4 um along the slope of the last tabulated
    # segment, which rises, so beyond SB99_LAM_MAX the dust-free curve (and A)
    # is extrapolation, not physics. The final camera point (1e4 um) has
    # essentially zero dust-free flux (~1e-35 of the peak) and is dropped.
    flux_dusty = flux_dusty.astype(float); flux_clear = flux_clear.astype(float)
    dead = flux_clear < 1e-20 * flux_clear.max()
    flux_dusty[dead] = np.nan; flux_clear[dead] = np.nan
    extrap = lam_rest > SB99_LAM_MAX

    A_mag = np.full_like(flux_dusty, np.nan)
    ok = (flux_clear > 0) & (flux_dusty > 0)
    A_mag[ok] = -2.5 * np.log10(flux_dusty[ok] / flux_clear[ok])

    # Bolometric fluxes. The images are I_nu, so integrate over FREQUENCY:
    # integrating I_nu d(lambda) (as the bh1 notebook did) weights each point by
    # lambda^2/c, letting the far-IR dust emission dominate and making the dusty
    # run look many times brighter than the dust-free one.
    trapz = getattr(np, 'trapezoid', None) or np.trapz
    nu = 2.998e14 / lam_rest                       # Hz (rest frame; a constant (1+z) cancels in ratios)
    def bol(f, sel=slice(None)):
        f, x = f[sel], nu[sel]
        g = np.isfinite(f)
        return trapz(f[g][::-1], x[g][::-1])  # nu decreases with lam, so reverse
    stellar = lam_rest < 8.0                       # rest-frame UV/optical/NIR: direct starlight
    L_dusty_bol, L_clear_bol = bol(flux_dusty), bol(flux_clear)
    absorbed = 1 - bol(flux_dusty, stellar) / bol(flux_clear, stellar)
    i_v = np.argmin(np.abs(lam - 0.55))
    print(f'{name}: z = {z:.3f}')
    print(f'   Bolometric flux, int I_nu dnu (dusty):     {L_dusty_bol:.6e}')
    print(f'   Bolometric flux, int I_nu dnu (dust-free): {L_clear_bol:.6e}')
    print(f'   Energy check, dusty/dust-free bolometric (~1 expected): {L_dusty_bol / L_clear_bol:.3f}')
    print(f'   Fraction of starlight (rest < 8 um) absorbed by dust: {absorbed * 100:.1f}%')
    print(f'   A at observed {lam[i_v]:.3f} um (rest {lam_rest[i_v]:.3f} um): {A_mag[i_v]:.2f} mag')

    # --- flux + attenuation vs observed wavelength ---
    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    axes[0].plot(lam, flux_dusty, label='With dust', color='firebrick')
    axes[0].plot(lam, flux_clear, label='Dust-free (intrinsic)', color='steelblue', linestyle='--')
    axes[0].set_xscale('log'); axes[0].set_yscale('log')
    axes[0].set_ylabel(r'Integrated flux [arb. units, $\sum_{pix}$ I$_\nu$]')
    axes[0].legend()
    axes[0].set_title(f'Dusty vs. Intrinsic (Dust-Free) Continuum Flux (Observed Frame, z={z:.3f})')
    axes[1].plot(lam, A_mag, color='darkorange')
    axes[1].axhline(0, color='gray', linewidth=0.8)
    axes[1].set_xscale('log')
    axes[1].set_xlabel(r'Observed-frame wavelength [$\mu$m]')
    axes[1].set_ylabel('Dust attenuation A [mag]')
    axes[1].set_title(f'Dust Attenuation vs. Observed-Frame Wavelength (z={z:.3f})')
    x0 = SB99_LAM_MAX * (1 + z)
    for ax in axes:
        ax.axvspan(x0, np.nanmax(lam[~dead]) * 1.05, color='0.85', zorder=0)
    axes[0].text(x0 * 1.1, 0.04, 'dust-free = extrapolated\nstellar template', transform=axes[0].get_xaxis_transform(),
                 fontsize=8, color='0.35', va='bottom')
    fig.tight_layout()
    fig.savefig(os.path.join(run, 'rtgen_dust_obscuration.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)

    # --- spatial attenuation map at observed V band ---
    img_d = im_dusty['image'][i_v]
    img_c = im_clear['image'][i_v]
    A_map = np.full_like(img_d, np.nan)
    m = (img_c > 0) & (img_d > 0)
    A_map[m] = -2.5 * np.log10(img_d[m] / img_c[m])
    half_x = im_dusty['nx'] * im_dusty['sizepix_x'] / 2.0 / KPC_CM
    half_y = im_dusty['ny'] * im_dusty['sizepix_y'] / 2.0 / KPC_CM
    fig2, ax2 = plt.subplots(figsize=(6, 5))
    im_show = ax2.imshow(A_map, origin='lower', cmap='inferno', vmin=0,
                         extent=[-half_x, half_x, -half_y, half_y])
    ax2.set_title(f'Spatial Dust Attenuation Map @ {lam[i_v]:.2f} $\\mu$m (observed frame, z={z:.3f})')
    ax2.set_xlabel('kpc'); ax2.set_ylabel('kpc')
    fig2.colorbar(im_show, ax=ax2, fraction=0.046, pad=0.04, label='A [mag]')
    fig2.tight_layout()
    fig2.savefig(os.path.join(run, 'rtgen_dust_obscuration_map.png'), dpi=150, bbox_inches='tight')
    plt.close(fig2)

    return {'name': name, 'z': z, 'lam': lam[~extrap], 'A_mag': A_mag[~extrap]}


def main(args):
    if args:
        runs = [a if os.path.isdir(a) else os.path.join(HERE, a) for a in args]
    else:
        runs = sorted(glob.glob(os.path.join(HERE, 'sys*_snap*')))
    results = []
    for run in runs:
        if not os.path.exists(os.path.join(run, 'dust_obscuration_analysis', 'image_cont_nodust.out')):
            print(f'{os.path.basename(run)}: no dust-free image yet, skipping')
            continue
        results.append(analyse(run))

    if len(results) > 1:
        fig, ax = plt.subplots(figsize=(8, 5))
        for r in results:
            ax.plot(r['lam'], r['A_mag'], label=f"{r['name']}  (z={r['z']:.2f})")
        ax.axhline(0, color='gray', linewidth=0.8)
        ax.set_xscale('log')
        ax.set_xlabel(r'Observed-frame wavelength [$\mu$m]')
        ax.set_ylabel('Dust attenuation A [mag]')
        ax.set_title(f'Dust attenuation, triple-AGN systems (rest-frame $\\lambda$ < {SB99_LAM_MAX:.0f} $\\mu$m)')
        ax.legend(fontsize=8)
        fig.tight_layout()
        out = os.path.join(HERE, 'dust_attenuation_all_systems.png')
        fig.savefig(out, dpi=150, bbox_inches='tight')
        print(f'Wrote {out}')


if __name__ == '__main__':
    main(sys.argv[1:])
