#!/usr/bin/env bash
# Submit the four RTGen stages for each system as a SLURM dependency chain:
#   1 SPH interpolation -> 2 stellar modelling -> 3 dusty RT (+ H2/HI, CO)
#   -> 4 dust-free image, dust-obscuration plots and default band luminosities
#
# Usage: scripts/submit_all.sh [SYSNUM ...]     (default: every sys*_snap* in $RUNS_DIR)
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"

if [ $# -gt 0 ]; then
  dirs=(); for n in "$@"; do
    m=( "$RUNS_DIR"/sys${n}_snap* )
    [ ${#m[@]} -eq 1 ] && [ -d "${m[0]}" ] || { echo "ERROR: need exactly one $RUNS_DIR/sys${n}_snap* (found ${#m[@]})" >&2; exit 1; }
    dirs+=("${m[0]}")
  done
else
  dirs=( "$RUNS_DIR"/sys*_snap* )
fi

for d in "${dirs[@]}"; do
  n=$(basename "$d" | sed 's/^sys\([^_]*\)_.*/\1/')
  j1=$(cd "$d/sph_interpolation" && sbatch --parsable -J rt${n}_sph submit_interpolation.sh)
  j2=$(cd "$d/stars_model" && sbatch --parsable -J rt${n}_star --dependency=afterok:$j1 submit_stellar_modelling.sh)
  j3=$(cd "$d" && sbatch --parsable -J rt${n}_rt --dependency=afterok:$j2 submit_master_pipeline.sh)
  j4=$(cd "$d" && sbatch --parsable -J rt${n}_dust --dependency=afterok:$j3 submit_dust_obscuration.sh)
  echo "$(basename "$d"): sph=$j1  stars=$j2  rt=$j3  dustfree=$j4"
done
