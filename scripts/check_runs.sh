#!/usr/bin/env bash
# Show which outputs exist for every run in $RUNS_DIR
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
for d in "$RUNS_DIR"/sys*_snap*; do
  n=0; for f in image_cont.out dust_obscuration_analysis/image_cont_nodust.out \
                rtgen_dust_obscuration.png rtgen_dust_obscuration_map.png; do
    [ -e "$d/$f" ] && n=$((n+1)); done
  echo "$(basename "$d"): $n/4 outputs"
done
