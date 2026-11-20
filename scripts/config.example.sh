# ---------------------------------------------------------------------------
# RTGen-tripleAGN configuration.
# Copy to scripts/config.sh (git-ignored) and edit:   cp config.example.sh config.sh
# Every helper script in scripts/ sources scripts/config.sh.
# ---------------------------------------------------------------------------

# Where the run directories go (~9 GB each: use scratch, not home)
export RUNS_DIR=${RUNS_DIR:-/scratch/$USER/RTGen_runs}   # env value wins if already set

# Romulus25 snapshots (directory and file-name prefix; the snapshot number is appended)
export SIM_DIR=/path/to/Romulus/cosmo25/
export SIM_PREFIX=cosmo25p.768sg1bwK1BHe75

# Directory with TripleAGN-Catalog-*.pkl (read by selection/find_dust_analogues.py)
export TRIPLEAGN_CATALOG_DIR=/path/to/tripleAGNs/physical_systems/

# Simulation constants used by scripts/check_halo.py
export SIM_H=0.6777          # little h
export SIM_BOX_KPC=25000.    # comoving box size [kpc] (Romulus25: 25 Mpc)

# Commands every SLURM job runs before Python: load modules and activate the
# Python environment (numpy < 2, numba, scipy, pandas, matplotlib, pynbody >= 2).
# Written verbatim into <run>/rtgen_env.sh by make_run.sh.
RTGEN_ENV_SETUP='
module --force purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5
source $HOME/.venvs/rtgen/bin/activate
'

# Optional: SLURM account, if your cluster requires one (sbatch reads this variable)
# export SBATCH_ACCOUNT=def-yourpi
