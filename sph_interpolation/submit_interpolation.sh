#!/usr/bin/env bash
#SBATCH --cpus-per-task=64
#SBATCH --mem=192GB
#SBATCH --time=03:00:00
#SBATCH --output=job_sph_interpolation.out

export NUMBA_NUM_THREADS=64

cp ../input/input_params.py .
cp ../input/constants.py .

# Modules + Python environment, written by scripts/make_run.sh from scripts/config.sh
source ../rtgen_env.sh
python3 -u rewrite_sim.py
python3 -u sph_interpolation.py

rm input_params.py
rm constants.py
