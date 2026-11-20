# RTGen-tripleAGN

Dust radiative transfer of **triple-AGN systems in the Romulus25 cosmological
simulation**. It is built on the [RTGen](https://github.com/francescosinigaglia/RTGen-public)
pipeline (Sinigaglia et al. 2025) and the [RADMC-3D](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/)
code.

For each system, it produces:

* a **dusty** continuum image cube, 0.01–10⁴ µm (stars + dust emission), plus the H2/HI split and a CO line cube;
* a **dust-free** image cube from the same camera, which gives the intrinsic starlight;
* **dust attenuation** curves and maps (dusty vs. dust-free), against observed-frame wavelength;
* **band luminosities** for each member's host galaxy, plus an AGN estimate
  along the same line of sight, in any top-hat band, filter curve or X-ray band.

> The AGN are **not** light sources in the radiative transfer; only stars
> emit. See [docs/band_luminosities.md](docs/band_luminosities.md) for how AGN band
> luminosities are estimated and what that means.

## Documentation

| Document | Contents |
|---|---|
| [docs/pipeline_guide.md](docs/pipeline_guide.md) | **Step-by-step guide**: setup, choosing systems, creating and submitting runs, monitoring, plots, gotchas |
| [docs/band_luminosities.md](docs/band_luminosities.md) | Host-galaxy and AGN luminosity in a chosen band, from the output |
| [docs/CHANGES.md](docs/CHANGES.md) | What differs from upstream RTGen, and why |
| [docs/RTGen_original_README.md](docs/RTGen_original_README.md) | The upstream README |

## Quick start

```bash
git clone https://github.com/spencerlockwood/RTGen-tripleAGN.git && cd RTGen-tripleAGN
cp scripts/config.example.sh scripts/config.sh && $EDITOR scripts/config.sh   # paths, modules, venv
source scripts/config.sh

# 1. pick UGC 2369S-like triples (one per redshift bin)
python selection/find_dust_analogues.py --out $RUNS_DIR/picks.csv

# 2. one run directory per system
python scripts/runs_from_picks.py $RUNS_DIR/picks.csv | bash

# 3. submit the 4-stage SLURM chain per system (about 3–4 h each)
scripts/submit_all.sh
scripts/check_runs.sh

# 4. results
python postproc/plot_dust_obscuration.py $RUNS_DIR/sys*_snap*
python postproc/band_luminosity.py $RUNS_DIR/sys1_snap005529 --band-rest 0.5 0.6
python postproc/band_luminosity.py $RUNS_DIR/sys1_snap005529 --filter filters/HST/HST_WFC3_IR.F160W.dat
```

## Requirements

* [RADMC-3D](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/installation.html),
  with `radmc3d` on your `PATH`.
* Python 3.10+:
  * numpy < 2, scipy, numba, pandas, matplotlib, astropy;
  * pynbody ≥ 2.
  * Tested with numpy 1.26, scipy 1.12, numba 0.65, pynbody 2.3.1.
* A SLURM cluster:
  * stage 1 asks for 64 cores and 192 GB (it uses about 30 GB);
  * the RT stages ask for 32 cores and 20 GB;
  * each run takes about 9 GB of disk.
* Romulus25 snapshots with AHF halo catalogues, and the triple-AGN catalogue
  (`TripleAGN-Catalog-*.pkl`) for the selection step.

## Repository layout

```
input/              input_params.py (per-system placeholders), constants, opacities, molecular data
sph_interpolation/  stage 1: halo extraction, host-BH centring, SPH deposit onto the 256³ grid
stars_model/        stage 2: Starburst99 / Padova stellar spectra
dust_transfer/      stage 3: dust RT set-up (called by submit_master_pipeline.sh)
H2_HI_splitting/    stage 3: H2/HI split from the 1000 Å radiation field
line_transfer/      stage 3: CO line transfer set-up
submit_master_pipeline.sh, submit_dust_obscuration.sh, make_zero_dust_density.py
postproc/           plot_dust_obscuration.py, band_luminosity.py (+ original SED helpers)
filters/            HST / JWST / Euclid filter curves
compute_opacity/    Mie opacity code
selection/          find_dust_analogues.py
scripts/            config + helpers: make_run.sh, runs_from_picks.py, submit_all.sh, check_runs.sh, check_halo.py
docs/               documentation
```

## Citing

If you use this pipeline, please cite:

* RADMC-3D: [Dullemond et al. 2012, ASCL 1202.015](https://ui.adsabs.harvard.edu/abs/2012ascl.soft02015D/abstract)
* RTGen: [Sinigaglia et al. 2025, A&A](https://ui.adsabs.harvard.edu/abs/2024arXiv241208609S/abstract)
* Romulus25: Tremmel et al. 2017, MNRAS 470, 1121
* Starburst99: Leitherer et al. 1999, ApJS 123, 3

The original RTGen code is by F. Sinigaglia. The modifications are listed in
[docs/CHANGES.md](docs/CHANGES.md). The code and parameters come with no
warranty: always check the physical meaning of your results.
