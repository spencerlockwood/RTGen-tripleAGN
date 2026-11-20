# Step-by-step guide: running RTGen on Romulus25 triple AGN

This guide takes you from a triple AGN in the Romulus25 catalogue to:

* a dusty radiative-transfer continuum image (RADMC-3D), plus the H2/HI split and the CO line image;
* a **dust-free** image of the same galaxy, seen from the same camera;
* plots of dusty vs. dust-free flux and dust attenuation against observed-frame wavelength;
* band luminosities of each member's host galaxy and AGN ([band_luminosities.md](band_luminosities.md)).

All paths below are relative to the repository root unless they start with `$RUNS_DIR`.

---

## 0. How the pieces fit together

```
RTGen-tripleAGN/                  (this repository = the pipeline template)
├── input/input_params.py         parameters; @PLACEHOLDERS@ are filled per system by make_run.sh
├── sph_interpolation/            stage 1: cut the halo out of the snapshot, put gas/stars on a 256³ grid
├── stars_model/                  stage 2: Starburst99 stellar source spectra
├── dust_transfer/  H2_HI_splitting/  line_transfer/  submit_master_pipeline.sh
│                                 stage 3: RADMC-3D dust RT + image, H2/HI split, CO line image
├── make_zero_dust_density.py, submit_dust_obscuration.sh
│                                 stage 4: dust-free image + plots + default band luminosities
├── postproc/                     plot_dust_obscuration.py, band_luminosity.py (+ original helpers)
├── filters/                      HST / JWST / Euclid filter curves
├── compute_opacity/              dust opacity (Mie) code
├── selection/find_dust_analogues.py   picks the systems to run
├── scripts/                      config + helpers (make_run.sh, submit_all.sh, ...)
└── docs/                         this guide and the band-luminosity guide

$RUNS_DIR/                        (on scratch: about 9 GB per run)
├── systems.txt                   one line per run created by make_run.sh
└── sys<N>_snap<SNAP>/            one full copy of the pipeline per system
```

Each run directory is a full copy of the pipeline with the system's values
filled in. Runs are independent: you can delete, archive or re-run one without
touching the others.

---

## 1. One-time setup

### 1a. Clone

```bash
git clone https://github.com/<your-account>/RTGen-tripleAGN.git
cd RTGen-tripleAGN
```

The stellar-track tables in `stars_model/` make the clone about 3 GB.

### 1b. RADMC-3D

Install RADMC-3D and put `radmc3d` on your `PATH`. Follow the
[installation guide](https://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/manual_radmc3d/installation.html).

### 1c. Python environment

You need Python 3.10+ with:

* numpy **< 2**;
* scipy, numba, pandas, matplotlib, astropy;
* **pynbody ≥ 2**, for reading the snapshots and AHF halos.

For example:

```bash
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5      # or your cluster's equivalent
python -m venv ~/.venvs/rtgen
source ~/.venvs/rtgen/bin/activate
pip install "numpy<2" scipy numba pandas matplotlib astropy "pynbody>=2"
```

Tested with numpy 1.26, scipy 1.12, numba 0.65, pynbody 2.3.1.

> **Don't mix numpy 2 with pynbody 1.x / numba 0.58.** That environment
> crashed in stage 1. `rewrite_sim.py` needs the pynbody 2 halo API
> (`halo_numbers='v1'`).

### 1d. Configure

```bash
cp scripts/config.example.sh scripts/config.sh
$EDITOR scripts/config.sh
```

| Variable | What to set |
|---|---|
| `RUNS_DIR` | Where run directories are created. Use scratch: about 9 GB each |
| `SIM_DIR`, `SIM_PREFIX` | Snapshot directory, and the file name without the snapshot number, e.g. `cosmo25p.768sg1bwK1BHe75` |
| `TRIPLEAGN_CATALOG_DIR` | Directory with `TripleAGN-Catalog-*.pkl` |
| `SIM_H`, `SIM_BOX_KPC` | h and the comoving box size, used by `check_halo.py` |
| `RTGEN_ENV_SETUP` | The lines every job runs first: `module load …` and `source …/activate`. Copied into each run as `rtgen_env.sh` |
| `SBATCH_ACCOUNT` | Optional, if your cluster needs `--account` |

`scripts/config.sh` is git-ignored, so your paths stay out of the repository.
Then:

```bash
source scripts/config.sh     # do this in each new shell before using the scripts below
```

The scripts in `scripts/` source the config themselves. The `source` above is
for the Python tools that read `$TRIPLEAGN_CATALOG_DIR` and `$SIM_*`.

---

## 2. Choose the systems

```bash
python selection/find_dust_analogues.py --out $RUNS_DIR/picks.csv
```

This scores every triple against UGC 2369S (Ding et al. 2025) and picks the
best match in each of 5 redshift bins. It prints a report and writes the picks
to the CSV. Without `--out` it writes `$TRIPLEAGN_CATALOG_DIR/dust_analogues.csv`
and **overwrites** it each time.

How the score works:

* **Formula:** `score = T_mass + T_lum + T_sep`. Each term is the mean over
  members of `(log10 residual / σ)²`. Lower is better. σ sets how far off a value
  can be before it costs as much as one unit of score.
* **Members:** Host = N core, Neighbour 1 = SE, Neighbour 2 = SW (ordered by BH mass).
* **Mass:** each BH mass against Table 2 of the paper. σ = 0.3 dex.
* **Luminosity:** host L_bol must lie in 2.7e41–2.6e44 erg/s (no penalty inside).
  The companions are penalised only above 7.9e40 (SE) and 6.1e40 (SW). σ = 0.3 dex.
* **Separation:** projected Host–N1 against 1.78 kpc and Host–N2 against 1.18 kpc. σ = 0.15 dex.
* **Not scored:** stellar mass, SFR and gas are printed next to each pick, as a check.

Options:

| Option | Effect |
|---|---|
| `--companion-limits radio` | Looser companion L_bol limit (1.3e43), not assuming the companions are unobscured |
| `--z-targets 0.1 0.5 1 2` | Pick the best system near these redshifts instead of using 5 bins |
| `--n-systems N`, `--zmin`, `--zmax` | Number of bins and redshift range |
| `--max-sep-kpc 5` | Keep only compact systems |
| `--legacy` | Old five-term weighted score |
| `--phys-dir DIR` | Catalogue directory; overrides `$TRIPLEAGN_CATALOG_DIR` |

You can use any other triple instead. Each run needs:

* the snapshot and its redshift;
* the halo number;
* the three BH iords;
* the host position.

---

## 3. Create the run directories

### 3a. Generate the commands from the picks

```bash
python scripts/runs_from_picks.py $RUNS_DIR/picks.csv           # prints one make_run.sh line per pick
python scripts/runs_from_picks.py $RUNS_DIR/picks.csv | bash    # ... and runs them
```

Each line looks like this:

```bash
scripts/make_run.sh 2 003072 1.185 118727 1981808889 1981808959 1981838042 "-3864.7, -1807.4, 4519.8" 8.2678e+43 8.8571e+40 1.6509e+40
#                   |  |      |     |      host BH    N1 BH      N2 BH      host position [kpc]       L_bol host, N1, N2 [erg/s]
#                   |  |      |     halo_number (CSV column pynbody_haloid)
#                   |  |      redshift
#                   |  snapshot (zero-padded)
#                   system number -> directory sys2_snap003072
```

`make_run.sh` does five things:

* copies the pipeline into `$RUNS_DIR/sys<N>_snap<SNAP>/`;
* fills in the `@PLACEHOLDERS@` of `input/input_params.py`;
* writes `rtgen_env.sh` from `RTGEN_ENV_SETUP`;
* warns about any placeholder left unfilled, or a missing snapshot;
* appends a line to `$RUNS_DIR/systems.txt`.

It refuses to overwrite an existing directory.

**BH order matters:** host, N1, N2. The RT grid is centred on the first one.

**Re-using system numbers:** `submit_all.sh N` needs exactly one
`sysN_snap*` directory, so archive old runs first:

```bash
mkdir -p $RUNS_DIR/old_selection && mv $RUNS_DIR/sys*_snap* $RUNS_DIR/systems.txt $RUNS_DIR/old_selection/
```

### 3b. What the filled-in `input/input_params.py` means

| Parameter | Meaning |
|---|---|
| `input_galaxy_dir`, `input_galaxy_filename` | Snapshot to read (from `SIM_DIR`, `SIM_PREFIX`, SNAP) |
| `halo_number` | tangos `halo_number`, i.e. the **1-based** row of the `AHF_halos` file |
| `halo_numbering = 'v1'` | Makes pynbody use that 1-based numbering. **Don't change it.** 0-based (`'file-order'`) is off by one: a different halo, hundreds of kpc away |
| `bh_iords` | `[host, N1, N2]`. Stage 1 stops with an error if any of them is not in the halo |
| `host_pos_kpc` | Catalogue host position, printed only as a cross-check |
| `bh_lbol_erg_s` | AGN L_bol of `[host, N1, N2]`. Read only by `postproc/band_luminosity.py`; the AGN do not emit in the RT |
| `zreal`, `z_optical`, `z_alma` | Snapshot redshift |
| `tsnap = None` | Unused. The age of the universe is read from the snapshot |

Settings you might change (edit `input/input_params.py` in the repository
before creating runs, so every system gets the same values):

| Parameter | Default | Note |
|---|---|---|
| `nx, ny, nz` | 256 | Grid cells per side |
| `sizex, sizey, sizez` | 30 kpc | Grid box, centred on the host BH |
| `nphot*` | 1e7 | Monte Carlo photon packets |
| `dust_to_gas_model` | `metallicity_dependent` | Li, Narayanan & Davé (2019) relation |
| `fsil` | 0.55 | Silicate fraction (carbon = 1 − fsil) |
| `nthreads` | 64 | |

> **Camera:** `incl` in `input_params.py` is **not used**. The viewing angle and
> image size are hard-coded in the `radmc3d image …` lines of
> `submit_master_pipeline.sh` and `submit_dust_obscuration.sh`
> (`incl 90 secondorder sizepc 20000 npixx 256 npixy 256`).
>
> * If you change them, change **both** scripts the same way, or the dusty and
>   dust-free images won't be comparable.
> * With the axis swap in `sph_interpolation.py`, `incl 90` looks along
>   simulation +z: image x = sim x, image y = sim y (see [band_luminosities.md](band_luminosities.md) §2.4).

### 3c. (Optional) Check the halo before submitting

Stage 1 already checks that all three BHs are in the halo. To catch a wrong
halo number before a queue wait:

```bash
python scripts/check_halo.py 003072 118727 "-3864.7, -1807.4, 4519.8"
# cosmo25p...003072.z1.186.AHF_halos row 118727: Mvir = 3.14e+12 Msun, Rvir = 211 kpc, host BH offset = 2.7 kpc  ->  OK
```

The host offset should be a few kpc, well inside Rvir.

---

## 4. Submit

```bash
scripts/submit_all.sh            # every sys*_snap* in $RUNS_DIR
scripts/submit_all.sh 2 4        # only systems 2 and 4
```

For each system this submits four SLURM jobs. Each one waits for the previous
one to succeed (`--dependency=afterok`):

| # | Job | Script (run from) | What it does | Resources | Typical time |
|---|---|---|---|---|---|
| 1 | `rtN_sph` | `sph_interpolation/submit_interpolation.sh` | `rewrite_sim.py`: loads only the halo; finds the 3 BHs and centres on the host with periodic unwrapping; removes BH particles from the stars; writes `time_and_z.txt` and `bh_info.txt`. `sph_interpolation.py`: deposits gas and stars on the grid | 64 CPU, 192 GB (uses about 30 GB) | 15–35 min |
| 2 | `rtN_star` | `stars_model/submit_stellar_modelling.sh` | Starburst99 spectra for the stellar sources | 16 CPU, 5 GB | 1–2 min |
| 3 | `rtN_rt` | `submit_master_pipeline.sh` | Dust opacities and densities, `radmc3d mctherm`, continuum image → `image_cont.out`; 1000 Å field → H2/HI split; CO line image → `image_line.out` | 32 CPU, 20 GB | 1.5–2.3 h |
| 4 | `rtN_dust` | `submit_dust_obscuration.sh` | `dust_density_zerodust.inp` (dust = 1e-60); re-render in `dust_obscuration_analysis/` with the **same** dust temperature, stars and camera → `image_cont_nodust.out`; plots; `band_luminosity_V_rest.txt` | 32 CPU, 20 GB | 1–1.3 h |

So a system takes about 3–4 h once its jobs start. Edit the `#SBATCH` lines at
the top of each script if your cluster needs different limits.

### Running or re-running one stage

Each stage is an ordinary sbatch script. Submit it from the directory shown in
the table, e.g. to redo only the dust-free run of system 3:

```bash
cd $RUNS_DIR/sys3_snap002816 && sbatch submit_dust_obscuration.sh
```

Stage 4 reuses `dust_density_zerodust.inp` if it already exists. Delete it if
you re-ran stage 3.

---

## 5. Monitor and check

```bash
squeue -u $USER
sacct -u $USER -S <start-date> -X -o JobID,JobName%10,Elapsed,State
scripts/check_runs.sh            # outputs present per run (4/4 = finished)
```

If one stage fails, the later ones stay pending with reason
`DependencyNeverSatisfied`. `scancel` them, fix the problem, and resubmit from
that stage.

| Stage | Log (inside the run directory) |
|---|---|
| 1 | `sph_interpolation/job_sph_interpolation.out` |
| 2 | `stars_model/job.out` |
| 3 | `job.out` (also `radmc3d.out`) |
| 4 | `job_dust_obscuration.out` |

What to check in the stage 1 log:

* each BH's offset from the grid centre should match the catalogue separations;
* the redshift should be the snapshot's;
* `Gas/Star |x|,|y|,|z| max` should be about sizex/2 (15 kpc).

RADMC-3D may warn about the partition-function temperature. That warning is
harmless: it comes from empty cells with T = 0.

A finished run contains:

* `image_cont.out` (dusty image)
* `dust_obscuration_analysis/image_cont_nodust.out` (dust-free image)
* `rtgen_dust_obscuration.png`, `rtgen_dust_obscuration_map.png`
* `band_luminosity_V_rest.txt`
* `image_line.out`: the CO line cube

---

## 6. Plots

Stage 4 makes the per-system plots. To remake them, and to build the
comparison of all systems:

```bash
cd $RUNS_DIR
python /path/to/RTGen-tripleAGN/postproc/plot_dust_obscuration.py                 # every finished sys*_snap* here
python /path/to/RTGen-tripleAGN/postproc/plot_dust_obscuration.py sys2_snap003072 # specific runs
python /path/to/RTGen-tripleAGN/postproc/plot_dust_obscuration.py old_selection/sys*
```

| Output | Content |
|---|---|
| `<run>/rtgen_dust_obscuration.png` | Top: total flux, dusty vs. dust-free. Bottom: A = −2.5 log10(F_dusty/F_dust-free). Both against observed-frame λ_obs = λ_rest (1+z) |
| `<run>/rtgen_dust_obscuration_map.png` | Map of the attenuation at observed 0.55 µm |
| `./dust_attenuation_all_systems.png` | Attenuation curves of every run given, rest-frame < 160 µm, written to the **current directory**. Overwritten each time the script runs with more than one run |

For each run it also prints:

* the bolometric fluxes and their ratio (should be about 0.7–1.1);
* the fraction of starlight (rest < 8 µm) absorbed by dust;
* A at observed 0.55 µm.

How to read the plots:

* **Observed 0.55 µm is rest-frame UV at high z.** At z = 2.5 it is ~0.16 µm
  rest, where attenuation is larger than A_V. To get A_V, compare at fixed
  rest wavelength: `band_luminosity.py --band-rest 0.5 0.6` reports A_band.
* **Negative A in the far-IR is real.** It is dust re-emitting the absorbed
  starlight; the dust-free run has no dust to emit.
* **Grey region (rest > 160 µm):** Starburst99 stops at 160 µm. RTGen
  extrapolates beyond that, so the dust-free curve there is not physical.
* **Bolometric integrals are over ν**, because the images are I_ν.

---

## 7. Band luminosities

See [band_luminosities.md](band_luminosities.md). In short:

```bash
python postproc/band_luminosity.py $RUNS_DIR/sys2_snap003072 --band-rest 0.48 0.54
python postproc/band_luminosity.py $RUNS_DIR/sys2_snap003072 --filter filters/JWST/<filter>.dat
python postproc/band_luminosity.py $RUNS_DIR/sys2_snap003072 --xray 2 10
```

---

## 8. Limitations and gotchas

* **No AGN emission in the RT.** Only stars are light sources; the BHs only set
  the grid centre. Attenuation here is of starlight. AGN band luminosities are
  an add-on estimate ([band_luminosities.md](band_luminosities.md) §3).
* **Halo numbering.** tangos `halo_number` = AHF 1-based row = pynbody
  `halo_numbers='v1'`. `rewrite_sim.py` forces `priority=['AHFCatalogue']`,
  because pynbody 2 can otherwise pick the `grp`-based catalogue.
* **Host M200 in the triple-AGN catalogue may be wrong.** In the catalogue used
  here, halo properties were merged one row off, so "host M200" belongs to a
  neighbouring halo. Use `scripts/check_halo.py` for the real Mvir.
* **Dust-free floor.** `make_zero_dust_density.py` uses 1e-60 per species.
  1e-30 is too high: it is comparable to the real dust density in many cells
  and makes the dust-free run glow in the far-IR.
* **Fixed camera** (`incl 90`, 20 kpc, 256²), in both image commands (§3b).
* **Starlight above 160 µm rest is extrapolated** (§6).
* **Disk space:** about 9 GB per run. Keep runs on scratch, and copy the
  results you want to keep (PNGs, `image_cont*.out`, band-luminosity outputs)
  somewhere permanent: scratch is purged.
* **Outputs that get overwritten:** the picks CSV (without `--out`) and
  `dust_attenuation_all_systems.png`.
