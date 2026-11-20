#!/usr/bin/env bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=20GB
#SBATCH --time=2:00:00
#SBATCH --output=job_dust_obscuration.out

# Modules + Python environment, written by scripts/make_run.sh from scripts/config.sh
source rtgen_env.sh

# --------------------------------------
# 0) Generate the zero-dust density file, if not already done
# --------------------------------------
if [ ! -f dust_density_zerodust.inp ]; then
    echo "Generating dust_density_zerodust.inp ..."
    python3 make_zero_dust_density.py
else
    echo "dust_density_zerodust.inp already exists, skipping generation."
fi

# --------------------------------------
# 1) Set up the dust-free comparison directory
# --------------------------------------
mkdir -p dust_obscuration_analysis
cd dust_obscuration_analysis

# -f so re-running this script doesn't fail on already-existing symlinks
ln -sf ../amr_grid.inp .
ln -sf ../dustopac.inp .
ln -sf ../dustkappa_*.inp .
ln -sf ../stellarsrc_density.inp .
ln -sf ../stellarsrc_templates.inp .
ln -sf ../wavelength_micron.inp .
ln -sf ../camera_wavelength_micron.inp .
ln -sf ../radmc3d.inp .
ln -sf ../dust_temperature.dat .

cp ../dust_density_zerodust.inp dust_density.inp

# --------------------------------------
# 2) Run the dust-free continuum image
#    Same flags as the real dusty run in submit_master_pipeline.sh, so the
#    two images are directly comparable -- no mctherm needed since there's
#    no dust to heat.
# --------------------------------------
radmc3d image loadlambda incl 90 secondorder sizepc 20000 npixx 256 npixy 256 setthreads 32
mv image.out image_cont_nodust.out

cd ..

echo "Done. Dust-free image at dust_obscuration_analysis/image_cont_nodust.out"
# --------------------------------------
# 3) Dusty vs dust-free flux / attenuation plots for this run
# --------------------------------------
python3 postproc/plot_dust_obscuration.py "$PWD"

# --------------------------------------
# 4) Band luminosities of the host galaxies (+ AGN estimate), default band
#    (edit or re-run by hand; see docs/band_luminosities.md)
# --------------------------------------
python3 postproc/band_luminosity.py "$PWD" --band-rest 0.5 0.6 > band_luminosity_V_rest.txt || true
