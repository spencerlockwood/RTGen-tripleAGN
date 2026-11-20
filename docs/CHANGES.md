# Changes relative to upstream RTGen

Upstream: <https://github.com/francescosinigaglia/RTGen-public> (Sinigaglia et al. 2025).
This fork adapts it to **cosmological** Romulus25 snapshots and to **triple-AGN** systems.
It also adds a dust-free comparison run and band-luminosity tools. Everything
in the RT itself (opacities, Starburst99 sources, RADMC-3D calls, H2/HI split,
CO lines) is unchanged unless listed here.

## Pipeline

| File | Change | Why |
|---|---|---|
| `sph_interpolation/rewrite_sim.py` | Loads the snapshot lazily and extracts one AHF halo (`priority=['AHFCatalogue']`, `halo_numbers='v1'`) | The full 25 Mpc box does not fit in memory. `'v1'` matches tangos `halo_number`; 0-based numbering is off by one |
| | Age of the universe and z read from the snapshot | The upstream code assumed an isolated run starting at z = 4 and a hard-coded `tsnap` |
| | Centres on the **host BH** (by iord), with periodic minimum-image unwrapping | BHs need not sit at the halo centre; haloes can straddle the box edge |
| | Checks that all three `bh_iords` are in the halo; prints each BH's offset; writes `bh_info.txt` | Catch a wrong halo early; positions for apertures |
| | Excludes BH particles (`tform < 0`) from the stars | Otherwise each BH (up to ~5e8 M☉) is treated as an old stellar population at the nucleus |
| `sph_interpolation/sph_interpolation.py` | Smoothing lengths from one batched, multi-threaded `KDTree.query` | The per-particle Python loop was the bottleneck |
| | Stars use the fixed softening length (`softlength`) instead of gas smoothing lengths | Stars are collisionless; the upstream line passed the gas `hharr` to the stellar deposit |
| | Timing and diagnostic prints | |
| `input/input_params.py` | Placeholders for snapshot, halo, BH iords, host position, L_bol, redshift; `tgas_eq_tdust = 0`; `nthreads = 64` | Filled in per system by `scripts/make_run.sh` |
| `*/submit_*.sh`, `submit_master_pipeline.sh` | Source `rtgen_env.sh` (modules + venv) instead of `conda run -n numba-env`; resources for a 64-core node (stage 1: 192 GB, 3 h; stage 3: 8 h) | Portable across clusters |

## Added

| File | Purpose |
|---|---|
| `make_zero_dust_density.py` | Dust-free density file (1e-60 per species) |
| `submit_dust_obscuration.sh` | Stage 4: dust-free image with the same camera and dust temperature, then the plots and default band luminosities |
| `postproc/plot_dust_obscuration.py` | Dusty vs. dust-free flux, attenuation curves and maps; bolometric check integrated over ν |
| `postproc/band_luminosity.py` | Host-galaxy and AGN band luminosities per member (see `docs/band_luminosities.md`) |
| `selection/find_dust_analogues.py` | Picks UGC 2369S-like triple AGN from the Romulus25 catalogue |
| `scripts/` | `config.example.sh`, `make_run.sh`, `runs_from_picks.py`, `submit_all.sh`, `check_runs.sh`, `check_halo.py` |
| `docs/` | Step-by-step guide, band-luminosity guide, this file, the original README |

## Removed

* `.DS_Store` files, editor backups (`*~`) and a committed `__pycache__/*.pyc`.
