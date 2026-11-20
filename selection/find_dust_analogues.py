
"""
find_dust_analogues.py
======================

Picks ~5 triple AGN from the Romulus25 physical catalogues that most closely
resemble UGC 2369S (Ding et al. 2025), spread across redshift, as targets for
the dust / inclination study.

THE TEMPLATE: UGC 2369S
-----------------------
Southern component of the UGC 2369 interacting pair, a LIRG at z = 0.0318
(log L_IR/L_sun = 11.57) with three cores inside ~2 kpc. Properties from
Ding et al. Table 2, with members ordered by black hole mass so they line up
with the catalogue's Host / Neighbour1 / Neighbour2 ranking:

  core   M_BH (Msun)   M* (Msun)   SFR (Msun/yr)   E(B-V)   sigma* (km/s)
  N      4.73e8        2.32e10     14.70           0.78     220
  SE     1.52e8        6.5e9        3.55           0.69     170
  SW     9.08e6        4.8e9        2.63           0.46      89

  L(2-8 keV)  : N 1.33e40 (detected), SE < 3.97e39, SW < 3.05e39 erg/s
                (Table 2; SE and SW are star-forming in every diagnostic, and
                the paper's conclusion rules out AGN in them)
  separations : projected. Inner binary N-SW a_in = 1.18 kpc, outer intruder
                SE at R_out = 1.78 kpc (Sec. 4.2). The inner binary is
                dominated by N, so R_out ~ the N-SE distance.
  obscuration : N_H ~ 1e25.5 cm^-2 towards the N core, A_V approaching 10 mag
  cold H2     : (0.4 - 2.1)e9 Msun for the system

Caveats on the template, which is why the defaults are what they are:
* M* is measured in a 1.5" aperture, i.e. ~0.5 kpc RADIUS at this redshift.
  The catalogue's BH_Mstar_1 is within 1 kpc radius, so simulated values are
  expected to run higher. Treated as an order-of-magnitude match only.
* The luminosity is deeply uncertain. Observed L(2-8 keV) = 1.33e40 erg/s for
  the N core, but the source is Compton-thick and the radio-predicted
  L(2-10 keV) is ~1.3e43 erg/s (Ding et al. Table 3), a factor of ~1000
  higher. Luminosity is therefore scored against the RANGE
  TEMPLATE['Lbol_range'] rather than a single value: the Host's
  L_bol costs nothing anywhere inside the bracket, and is penalised by its
  squared log10 distance to the nearer edge outside it. Use --no-lbol to drop
  the term.
* Object 2 (SW) may be a young in-situ star-forming clump rather than a core
  with its own black hole, so its 9e6 Msun entry is the least secure.

WHAT IS MATCHED
---------------
Three terms, member by member (Host = N, Neighbour1 = SE, Neighbour2 = SW):

  score = T_mass + T_lum + T_sep,    each T = mean over members of (r / sigma)^2

  T_mass : r = log10(M_BH,sim / M_BH,obs) for each of the three BHs.
           sigma = 0.3 dex, the scatter of the M_BH-sigma relation the
           observed masses come from.
  T_lum  : Host L_bol against the N-core range (zero inside it). Neighbour1
           and Neighbour2 against the SE / SW upper limits: zero below, r =
           log10(L / limit) above. sigma = 0.3 dex.
  T_sep  : projected (2D) separations, Host-N1 vs 1.78 kpc and Host-N2 vs
           1.18 kpc. sigma = 0.15 dex for the unknown viewing angle.

A score of 1 means every quantity is off by about one sigma on average; lower
is better. Members with a missing value are left out of that term's mean.
Stellar mass, SFR and gas are NOT scored; they are printed next to each pick
as a check that it is gas-rich enough for the dust study.

Companion L_bol limits (--companion-limits):
  xray  (default): the observed X-ray limits x 20, i.e. SE < 7.9e40 and
                   SW < 6.1e40 erg/s. Assumes the companions are not hidden
                   behind Compton-thick gas.
  radio          : both companions are >= 19x fainter than N at 33 GHz, which
                   dust does not affect; scaling the N core's radio-predicted
                   luminosity gives L_bol < 1.3e43 erg/s for each.

The previous five-term weighted score (BH mass, M*, 3D separation, SFR,
brightest-member L_bol) is still available with --legacy.

Because the systems are wanted "at different redshifts", the selection takes
the best-scoring system in each of N redshift bins rather than the N best
overall, which would otherwise cluster at one epoch. Pass --z-targets to snap
to specific redshifts instead.

Romulus25's lowest snapshot is z ~ 0.05, so the z = 0.0318 of UGC 2369S itself
is not reachable; the lowest bin is the closest available analogue.

Usage
-----
Command line:
    python find_dust_analogues.py
    python find_dust_analogues.py --n-systems 5 --zmax 1.5 --max-sep-kpc 5
    python find_dust_analogues.py --z-targets 0.06 0.25 0.5 1.0 2.0
    python find_dust_analogues.py --companion-limits radio    # looser companion L_bol limits
    python find_dust_analogues.py --legacy                    # previous five-term score
Jupyter: paste the file into a cell (edit NOTEBOOK_SETTINGS at the bottom), or
    from find_dust_analogues import select;  picks, cat = select(n_systems=5)
"""

import os
import sys
import glob
import argparse
from types import SimpleNamespace

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# UGC 2369S, ordered by BH mass so index 0 = Host, 1 = Neighbour1, 2 = Neighbour2
TEMPLATE = dict(
    name='UGC 2369S (Ding et al. 2025)',
    z=0.0318,
    core_labels=['N', 'SE', 'SW'],
    BH_mass=[4.73e8, 1.52e8, 9.08e6],        # Msun, M_BH-sigma relation
    Mstar=[2.32e10, 6.5e9, 4.8e9],           # Msun, 1.5" aperture (~0.5 kpc radius)
    SFR=[14.70, 3.55, 2.63],                 # Msun/yr per core
    EBV=[0.78, 0.69, 0.46],
    sigma_star=[220.3, 170.0, 89.4],         # km/s
    sep_inner_kpc=1.18,                      # projected; inner binary N-SW
    sep_outer_kpc=1.78,                      # projected; SE intruder to the inner binary
    # luminosity bracket: observed X-ray with k_bol ~ 20, to radio-predicted
    Lbol_range=(2.7e41, 2.6e44),             # erg/s; scored as a range, no penalty inside
    LX_28keV=[1.33e40, 3.97e39, 3.05e39],    # erg/s; SE and SW are upper limits
    k_bol=20,                                # bolometric correction used throughout
    # companion (SE, SW) L_bol upper limits, see --companion-limits
    Lbol_companion_max=dict(xray=(3.97e39 * 20, 3.05e39 * 20),
                            radio=(1.3e43, 1.3e43)),
    M_H2_cold=(0.39e9, 2.07e9),              # Msun, system total
    NH_primary=10 ** 25.5,                   # cm^-2 towards the N core
    L_IR_Lsun=10 ** 11.57,
)

# Uncertainty (dex) each residual is divided by in the default score.
SIGMA = dict(mass=0.3, lum=0.3, sep=0.15)

# Weights for the old five-term score (--legacy only).
WEIGHTS = dict(bh_mass=1.0, mstar=0.5, separation=1.0, sfr=0.4, lbol=1.0)

BH_COLS = ['BH_mass', 'Neighbour1BH_mass', 'Neighbour2BH_mass']
LBOL_COLS = ['BH_Lbol', 'Neighbour1BH_Lbol', 'Neighbour2BH_Lbol']
SEP_COLS = ['Separation_Host_N1_3D_kpc', 'Separation_Host_N2_3D_kpc',
            'Separation_N1_N2_3D_kpc']
SEP2D_COLS = ['Separation_Host_N1_2D_kpc', 'Separation_Host_N2_2D_kpc']  # vs outer, inner
TERM_COLS = ['T_mass', 'T_lum', 'T_sep']
ID_COLS = ['BH_id', 'Neighbour1BH_id', 'Neighbour2BH_id']

DEFAULTS = dict(
    # Directory holding TripleAGN-Catalog-*.pkl. Set TRIPLEAGN_CATALOG_DIR or pass --phys-dir.
    phys_dir=os.environ.get('TRIPLEAGN_CATALOG_DIR', '.'),
    threshold=1e40,
    n_systems=5,
    z_targets=None,
    z_tol=0.15,
    zmin=None,
    zmax=None,
    aperture_kpc=1,
    max_sep_kpc=None,
    score_lbol=True,
    companion_limits='xray',
    legacy=False,
    out=None,
    quiet=False,
)


def _member_cols(base, r):
    """e.g. ('Mstar', 1) -> the three columns holding BH_Mstar_1 per member."""
    return [f'BH_{base}_{r}', f'Neighbour1BH_{base}_{r}', f'Neighbour2BH_{base}_{r}']


# Rough spherical column density from the mean gas density in a sphere of
# radius r: N_H ~ rho * r / (mu m_H). Strictly a kpc-scale average, so it will
# fall far below the ~1e25 cm^-2 measured on the ~50 pc nuclear obscurer in
# UGC 2369S. Useful for ranking, not for predicting an observed N_H.
_MSUN_PER_KPC3_TO_G_PER_CM3 = 1.989e33 / (3.086e21) ** 3
_KPC_CM = 3.086e21
_MU_MH = 1.4 * 1.6726e-24


def estimate_NH(gas_density_msun_kpc3, radius_kpc):
    rho = np.asarray(gas_density_msun_kpc3, float) * _MSUN_PER_KPC3_TO_G_PER_CM3
    return rho * radius_kpc * _KPC_CM / _MU_MH


def _in_notebook():
    try:
        from IPython import get_ipython
        return get_ipython() is not None
    except ImportError:
        return False


# ----------------------------------------------------------------------
def load_catalogues(phys_dir, threshold):
    """Concatenate every per-snapshot triple AGN catalogue at this threshold."""
    tag = '0' if threshold == 0 else f'{threshold:.0e}'.replace('+', '')
    d = os.path.join(phys_dir, f'{tag}lum')
    paths = sorted(glob.glob(os.path.join(d, 'TripleAGN-Catalog-*.pkl')))
    if not paths:
        raise FileNotFoundError(
            f'No catalogues in {d}. Check --phys-dir and --threshold '
            '(the directory is named <threshold>lum, e.g. 1e40lum).')
    frames, faulty = [], []
    for p in paths:
        df = pd.read_pickle(p)
        if len(df) == 0 or not all(c in df.columns for c in BH_COLS):
            continue
        # faulty snapshot: tangos had no BH mass/mdot, so Halomap wrote zeros
        if float(np.mean(~(df[BH_COLS].to_numpy(float) > 0))) > 0.5:
            faulty.append(os.path.basename(p))
            continue
        frames.append(df)
    if faulty:
        print(f'Excluded {len(faulty)} faulty snapshot(s): {", ".join(faulty)}')
    if not frames:
        raise RuntimeError(f'Every catalogue in {d} was empty or faulty.')
    cat = pd.concat(frames, ignore_index=True)
    if 'snapshot' in cat.columns:
        cat['snapshot'] = cat['snapshot'].astype(str).str.zfill(6)
    print(f'Loaded {len(cat)} triple AGN from {len(frames)} snapshots '
          f'(threshold {tag}), z = {cat.z.min():.3f} - {cat.z.max():.3f}')
    return cat


def _log_resid(sim, lo, hi=None):
    """log10 residual of sim against a value (hi None) or a range [lo, hi],
    zero inside the range. lo or hi may be None for a one-sided limit."""
    with np.errstate(divide='ignore', invalid='ignore'):
        logsim = np.log10(np.where(sim > 0, sim, np.nan))
        if hi is None and lo is not None:
            return logsim - np.log10(lo)
        r = np.zeros_like(logsim)
        if lo is not None:
            r = r + np.minimum(logsim - np.log10(lo), 0.0)
        if hi is not None:
            r = r + np.maximum(logsim - np.log10(hi), 0.0)
        return np.where(np.isfinite(logsim), r, np.nan)


def _mean_sq(resid, sigma):
    """Mean of (r/sigma)^2 over the members that have a value; NaN if none do."""
    ok = np.isfinite(resid)
    n = ok.sum(axis=1)
    sq = np.where(ok, (resid / sigma) ** 2, 0.0).sum(axis=1)
    return np.where(n > 0, sq / np.maximum(n, 1), np.nan)


def score_systems(cat, args):
    """Member-by-member comparison of BH mass, L_bol and projected separation.
    Returns the total score and writes the three terms into cat[TERM_COLS]."""
    t = TEMPLATE
    m = cat[BH_COLS].to_numpy(float)
    r_mass = np.column_stack([_log_resid(m[:, i], t['BH_mass'][i]) for i in range(3)])
    cat['T_mass'] = _mean_sq(r_mass, SIGMA['mass'])

    terms = ['T_mass']
    if args.score_lbol:
        L = cat[LBOL_COLS].to_numpy(float)
        lo, hi = t['Lbol_range']
        lim1, lim2 = t['Lbol_companion_max'][args.companion_limits]
        r_lum = np.column_stack([_log_resid(L[:, 0], lo, hi),
                                 _log_resid(L[:, 1], None, lim1),
                                 _log_resid(L[:, 2], None, lim2)])
        cat['T_lum'] = _mean_sq(r_lum, SIGMA['lum'])
        terms.append('T_lum')
    else:
        cat['T_lum'] = 0.0

    s = cat[SEP2D_COLS].to_numpy(float)
    r_sep = np.column_stack([_log_resid(s[:, 0], t['sep_outer_kpc']),
                             _log_resid(s[:, 1], t['sep_inner_kpc'])])
    cat['T_sep'] = _mean_sq(r_sep, SIGMA['sep'])
    terms.append('T_sep')

    lim1, lim2 = t['Lbol_companion_max'][args.companion_limits]
    print(f'Scored on (each member, residual/sigma squared): BH mass '
          f'(sigma {SIGMA["mass"]} dex), '
          + (f'L_bol (Host in [{lo:.1e}, {hi:.1e}], N1 < {lim1:.1e}, N2 < {lim2:.1e}; '
             f'{args.companion_limits} limits; sigma {SIGMA["lum"]} dex), '
             if args.score_lbol else 'L_bol NOT scored, ')
          + f'projected separation (Host-N1 vs {t["sep_outer_kpc"]}, Host-N2 vs '
          f'{t["sep_inner_kpc"]} kpc; sigma {SIGMA["sep"]} dex)')
    return cat[terms].sum(axis=1, skipna=False).to_numpy()


def score_systems_legacy(cat, args):
    """Weighted sum of mean squared log10 residuals against the template."""
    score = np.zeros(len(cat))
    used, skipped = [], []

    def add(term, sim, target, weight, target_hi=None):
        """target_hi given: [target, target_hi] is a range, and the residual is
        the log10 distance to the nearer edge (zero anywhere inside it)."""
        nonlocal score
        if weight <= 0:
            return
        sim = np.atleast_2d(np.asarray(sim, float).T).T
        if sim.ndim == 1:
            sim = sim[:, None]
        target = np.atleast_1d(np.asarray(target, float))[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            logsim = np.log10(np.where(sim > 0, sim, np.nan))
            if target_hi is None:
                resid = logsim - np.log10(target)
            else:
                hi = np.atleast_1d(np.asarray(target_hi, float))[None, :]
                resid = (np.minimum(logsim - np.log10(target), 0.0)
                         + np.maximum(logsim - np.log10(hi), 0.0))
        ok = np.isfinite(resid)
        n_ok = ok.sum(axis=1)
        sq = np.where(ok, resid ** 2, 0.0).sum(axis=1)
        # mean over the quantities that landed, so a missing member value does
        # not NaN out a system or make its score artificially small
        contrib = np.where(n_ok > 0, weight * sq / np.maximum(n_ok, 1), np.nan)
        score = score + contrib
        used.append(f'{term} (w={weight})')

    # --- the three BH masses ---
    add('BH masses', cat[BH_COLS].to_numpy(float), TEMPLATE['BH_mass'],
        WEIGHTS['bh_mass'])

    # --- local stellar mass around each member ---
    mcols = _member_cols('Mstar', args.aperture_kpc)
    if all(c in cat.columns for c in mcols):
        add('stellar masses', cat[mcols].to_numpy(float), TEMPLATE['Mstar'],
            WEIGHTS['mstar'])
    else:
        skipped.append(f'stellar masses (no {mcols[0]} column; '
                       'was the builder run with --skip-env?)')

    # --- largest pairwise separation ---
    if all(c in cat.columns for c in SEP_COLS):
        add('max separation', cat[SEP_COLS].to_numpy(float).max(axis=1),
            TEMPLATE['sep_outer_kpc'], WEIGHTS['separation'])
    else:
        skipped.append('separations')

    # --- total SFR over the three members ---
    scols = _member_cols('SFR', args.aperture_kpc)
    if all(c in cat.columns for c in scols):
        tot = cat[scols].to_numpy(float).sum(axis=1)
        add('total SFR', tot, sum(TEMPLATE['SFR']), WEIGHTS['sfr'])
    else:
        skipped.append('SFR')

    # --- luminosity of the brightest member, against the UGC 2369S bracket ---
    # The bracket is for the N core (the only one with a luminosity estimate),
    # so it is compared with the brightest simulated member. Range, not a point:
    # see the module docstring for why the true value is uncertain by ~1000x.
    if args.score_lbol:
        if all(c in cat.columns for c in LBOL_COLS):
            lo, hi = TEMPLATE['Lbol_range']
            add(f'brightest L_bol in [{lo:.1e}, {hi:.1e}]',
                cat[LBOL_COLS].to_numpy(float).max(axis=1), lo, WEIGHTS['lbol'],
                target_hi=hi)
        else:
            skipped.append('L_bol (no luminosity columns)')

    print('Scored on: ' + ', '.join(used))
    for s in skipped:
        print(f'  NOT scored: {s}')
    return score


def pick_spread_in_z(cat, args):
    """Best-scoring system in each of n redshift bins (or nearest each target)."""
    if args.z_targets:
        picks = []
        for zt in args.z_targets:
            sub = cat[np.abs(cat.z - zt) <= args.z_tol]
            if sub.empty:
                print(f'  no system within dz={args.z_tol} of z={zt}; skipped')
                continue
            row = sub.loc[sub.score.idxmin()].copy()
            row['z_target'] = zt
            row['n_candidates_in_bin'] = len(sub)
            picks.append(row)
        return pd.DataFrame(picks).reset_index(drop=True)

    n = args.n_systems
    edges = np.linspace(cat.z.min(), cat.z.max() + 1e-9, n + 1)
    b = np.clip(np.digitize(cat.z, edges) - 1, 0, n - 1)
    picks = []
    for k in range(n):
        sub = cat[b == k]
        if sub.empty:
            print(f'  z bin {edges[k]:.2f}-{edges[k+1]:.2f} is empty; skipped')
            continue
        row = sub.loc[sub.score.idxmin()].copy()
        row['z_target'] = 0.5 * (edges[k] + edges[k + 1])
        row['n_candidates_in_bin'] = len(sub)
        picks.append(row)
    return pd.DataFrame(picks).reset_index(drop=True)


def report(picks, args):
    t = TEMPLATE
    print('\n' + '=' * 78)
    print(f'SELECTED {len(picks)} SYSTEMS  (best match to {t["name"]} per redshift bin)')
    print('=' * 78)

    mcols = _member_cols('Mstar', args.aperture_kpc)
    scols = _member_cols('SFR', args.aperture_kpc)
    gcols = _member_cols('GasDensity', args.aperture_kpc)
    zcols = _member_cols('GasMetallicity', args.aperture_kpc)

    for i, (_, r) in enumerate(picks.iterrows(), 1):
        snap = r.get('snapshot', '?')
        print(f'\n--- System {i}:  snapshot {snap},  z = {r.z:.3f}   '
              f'(score {r.score:.3f}, best of {int(r.get("n_candidates_in_bin", 1))} '
              'in this z bin) ---')
        if not args.legacy and all(c in r.index for c in TERM_COLS):
            print(f'  score terms : mass {r.T_mass:.2f} + L_bol {r.T_lum:.2f} + '
                  f'separation {r.T_sep:.2f}')
        print('  BH ids      : ' + ', '.join(
            str(int(r[c])) for c in ID_COLS if c in r.index and pd.notna(r[c])))
        if 'pynbody_haloid' in r.index and pd.notna(r.pynbody_haloid):
            print(f'  host halo   : {int(r.pynbody_haloid)}')
        print(f'  position    : ({r.BH_x:.1f}, {r.BH_y:.1f}, {r.BH_z:.1f}) proper kpc '
              '(host BH, box frame)')

        print(f'  {"":14s}{"Host":>12s}{"N1":>12s}{"N2":>12s}   |  UGC 2369S (N / SE / SW)')
        print(f'  {"M_BH":14s}' + ''.join(f'{r[c]:12.2e}' for c in BH_COLS) +
              '   |  ' + ' / '.join(f'{v:.1e}' for v in t['BH_mass']))
        lim1, lim2 = t['Lbol_companion_max'][args.companion_limits]
        print(f'  {"L_bol":14s}' + ''.join(f'{r[c]:12.2e}' for c in LBOL_COLS) +
              f'   |  {t["Lbol_range"][0]:.0e}-{t["Lbol_range"][1]:.0e} / <{lim1:.0e} / <{lim2:.0e}')
        if all(c in r.index for c in mcols):
            print(f'  {f"M* (<{args.aperture_kpc}kpc)":14s}' +
                  ''.join(f'{r[c]:12.2e}' for c in mcols) +
                  '   |  ' + ' / '.join(f'{v:.1e}' for v in t['Mstar']))
        if all(c in r.index for c in scols):
            print(f'  {f"SFR (<{args.aperture_kpc}kpc)":14s}' +
                  ''.join(f'{r[c]:12.2f}' for c in scols) +
                  '   |  ' + ' / '.join(f'{v:.1f}' for v in t['SFR']))

        if all(c in r.index for c in SEP_COLS):
            s = [r[c] for c in SEP_COLS]
            print(f'  separations : Host-N1 {s[0]:.2f}, Host-N2 {s[1]:.2f}, '
                  f'N1-N2 {s[2]:.2f} kpc   |  {t["sep_inner_kpc"]:.2f} / '
                  f'{t["sep_outer_kpc"]:.2f} kpc (inner / outer)  [3D]')
        if all(c in r.index for c in SEP2D_COLS):
            print(f'  projected   : Host-N1 {r[SEP2D_COLS[0]]:.2f} (vs {t["sep_outer_kpc"]:.2f}), '
                  f'Host-N2 {r[SEP2D_COLS[1]]:.2f} (vs {t["sep_inner_kpc"]:.2f}) kpc'
                  + ('  [scored]' if not args.legacy else ''))

        # dust-relevant diagnostics
        if all(c in r.index for c in gcols):
            nh = estimate_NH([r[c] for c in gcols], args.aperture_kpc)
            print('  gas density : ' + ', '.join(f'{r[c]:.2e}' for c in gcols) +
                  ' Msun/kpc^3')
            print('  est. N_H    : ' + ', '.join(f'{v:.1e}' for v in nh) +
                  f' cm^-2  ({args.aperture_kpc} kpc average; '
                  f'UGC 2369S nuclear ~{t["NH_primary"]:.0e})')
        if all(c in r.index for c in zcols):
            print('  gas Z       : ' + ', '.join(f'{r[c]:.4f}' for c in zcols) +
                  '  (absolute metal mass fraction; dust traces this)')
        if 'Halo_M200' in r.index and pd.notna(r.get('Halo_M200')):
            print(f'  host M200   : {r.Halo_M200:.2e} Msun')
        if 'n_BH_within_R_of_centre' in r.index and pd.notna(r.get('n_BH_within_R_of_centre')):
            n_bh = int(r.n_BH_within_R_of_centre)
            print(f'  BHs within 30 kpc of the centre: {n_bh}'
                  + ('  (isolated triple)' if n_bh == 3 else
                     f'  ({n_bh - 3} extra faint BH(s) nearby)'))

    print('\n' + '=' * 78)
    print('NOTES')
    print('=' * 78)
    print(f'* UGC 2369S is at z = {t["z"]}, below Romulus25\'s lowest snapshot '
          '(z ~ 0.05),')
    print('  so the lowest-z pick is the closest reachable analogue.')
    print('* Template M* is a 1.5" aperture (~0.5 kpc radius); the simulated values')
    print(f'  are within {args.aperture_kpc} kpc radius, so expect them to run higher.')
    if not args.legacy:
        print('* M*, SFR and gas are shown for checking, not scored. UGC 2369S is a LIRG')
        print('  (total SFR ~21 Msun/yr), so a gas-poor pick is a poor target for the dust study.')
        print('* Companion L_bol limits are the ' + args.companion_limits + ' ones '
              '(see --companion-limits). The paper rules out AGN in SE and SW.')
        print("* Projected separations use the catalogue's fixed 2D projection; a real")
        print('  viewing angle can shorten them further, hence the 0.15 dex allowance.')
    elif args.score_lbol:
        lo, hi = t['Lbol_range']
        print(f'* L_bol is scored as a range ({lo:.1e} - {hi:.1e} erg/s, brightest member):')
        print('  UGC 2369S is Compton-thick and its observed vs radio-predicted X-ray')
        print('  luminosities differ by ~1000x, so only values outside that bracket cost.')
    else:
        print('* L_bol was NOT scored (--no-lbol).')
    print('* The estimated N_H is a kpc-scale spherical average and will sit far below')
    print('  the ~1e25 cm^-2 measured on the ~50 pc nuclear obscurer. Use it to rank')
    print('  systems by how gas-rich they are, not as a predicted column.')


def select(**kwargs):
    """Find the analogues. Returns (picks, scored catalogue)."""
    unknown = set(kwargs) - set(DEFAULTS)
    if unknown:
        raise TypeError(f'Unknown option(s): {sorted(unknown)}. Valid: {sorted(DEFAULTS)}')
    args = SimpleNamespace(**{**DEFAULTS, **kwargs})

    cat = load_catalogues(args.phys_dir, args.threshold)

    if args.zmin is not None:
        cat = cat[cat.z >= args.zmin]
    if args.zmax is not None:
        cat = cat[cat.z <= args.zmax]
    if args.max_sep_kpc is not None and all(c in cat.columns for c in SEP_COLS):
        keep = cat[SEP_COLS].to_numpy(float).max(axis=1) <= args.max_sep_kpc
        print(f'Compactness cut (max separation <= {args.max_sep_kpc} kpc): '
              f'{int(keep.sum())}/{len(cat)} kept')
        cat = cat[keep]
    if cat.empty:
        raise SystemExit('No systems left after the redshift / compactness cuts.')

    cat = cat.copy()
    cat['score'] = (score_systems_legacy if args.legacy else score_systems)(cat, args)
    n_before = len(cat)
    cat = cat[np.isfinite(cat.score)]
    if len(cat) < n_before:
        print(f'  dropped {n_before - len(cat)} system(s) with no scorable values')
    if cat.empty:
        raise SystemExit('No systems could be scored.')

    picks = pick_spread_in_z(cat, args)
    if picks.empty:
        raise SystemExit('No systems selected.')
    if not args.quiet:
        report(picks, args)

    out = args.out or os.path.join(args.phys_dir, 'dust_analogues.csv')
    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    picks.to_csv(out, index=False)
    print(f'\nSaved the full rows (every catalogue column) to {out}')
    print('Columns you need to re-load each system: snapshot, BH_id, '
          'Neighbour1BH_id, Neighbour2BH_id, BH_x/y/z, pynbody_haloid.')
    return picks, cat


def main():
    p = argparse.ArgumentParser(
        description='Find Romulus25 triple AGN resembling UGC 2369S, for the dust study')
    p.add_argument('--phys-dir', default=DEFAULTS['phys_dir'])
    p.add_argument('--threshold', type=float, default=DEFAULTS['threshold'],
                   help='Which catalogue to search (default 1e40)')
    p.add_argument('--n-systems', type=int, default=DEFAULTS['n_systems'],
                   help='How many systems to pick, one per redshift bin (default 5)')
    p.add_argument('--z-targets', type=float, nargs='+', default=None,
                   help="Snap to these redshifts instead of binning")
    p.add_argument('--z-tol', type=float, default=DEFAULTS['z_tol'],
                   help='Redshift tolerance when using --z-targets (default 0.15)')
    p.add_argument('--zmin', type=float, default=None)
    p.add_argument('--zmax', type=float, default=None)
    p.add_argument('--aperture-kpc', type=int, choices=[1, 5], default=1,
                   help='Which local-environment aperture to use (default 1 kpc)')
    p.add_argument('--max-sep-kpc', type=float, default=None,
                   help='Optional hard cut on the largest pairwise separation. '
                        'UGC 2369S spans ~1.8 kpc, so e.g. 5 keeps only compact systems.')
    p.add_argument('--no-lbol', dest='score_lbol', action='store_false',
                   help='Do not score on luminosity')
    p.add_argument('--companion-limits', choices=['xray', 'radio'],
                   default=DEFAULTS['companion_limits'],
                   help='L_bol upper limits for the two companions: observed X-ray x 20 '
                        '(default) or the looser radio-scaled limit')
    p.add_argument('--legacy', action='store_true',
                   help='Use the previous five-term weighted score')
    # Old flag, kept so earlier commands still run: luminosity is now on by default.
    p.add_argument('--match-lbol', dest='match_lbol_legacy', action='store_true',
                   help=argparse.SUPPRESS)
    p.add_argument('--out', default=None, help='Write the selected systems to this CSV')
    opts = vars(p.parse_args())
    opts.pop('match_lbol_legacy')
    select(**opts)


if __name__ == '__main__':
    # Command-line options are honoured even inside Jupyter (%run script.py --zmax 1).
    # NOTEBOOK_SETTINGS is used only when the file is pasted into a cell / imported,
    # i.e. when it was not launched as a script. (The kernel's own sys.argv contains
    # unrelated flags, so we compare the launched file name instead.)
    _script = os.path.basename(globals().get('__file__', ''))
    _ran_as_script = bool(_script) and os.path.basename(sys.argv[0]) == _script
    if _in_notebook() and not _ran_as_script:
        NOTEBOOK_SETTINGS = dict(
            phys_dir=DEFAULTS['phys_dir'],
            threshold=1e40,
            n_systems=5,
            # z_targets=[0.06, 0.25, 0.5, 1.0, 2.0],   # specific target redshifts
            max_sep_kpc=None,   # e.g. 5 to keep only compact (UGC 2369S-like) systems
            score_lbol=True,   # False to ignore luminosity
            companion_limits='xray',   # or 'radio'
            legacy=False,   # True for the previous five-term score
        )
        picks, catalogue = select(**NOTEBOOK_SETTINGS)
    else:
        main()