#!/usr/bin/env bash
# Create one RTGen run directory for a triple-AGN system.
#
# Usage:
#   scripts/make_run.sh SYSNUM SNAP Z HALO HOST_IORD N1_IORD N2_IORD "X, Y, Z" [LBOL_HOST LBOL_N1 LBOL_N2]
#
#   SYSNUM   your label for the system (1, 2, ...); the directory is sys<SYSNUM>_snap<SNAP>
#   SNAP     snapshot number, zero-padded as in the file name (e.g. 003072)
#   Z        snapshot redshift
#   HALO     tangos halo_number = 1-based AHF row (column pynbody_haloid of the picks CSV)
#   *_IORD   BH particle ids, in the order host, N1, N2 (the grid is centred on the host)
#   "X, Y, Z" host BH position [proper kpc, box frame] (cross-check only)
#   LBOL_*   optional AGN bolometric luminosities [erg/s] for postproc/band_luminosity.py
#
# scripts/runs_from_picks.py prints these commands for every row of a picks CSV.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"

if [ $# -ne 8 ] && [ $# -ne 11 ]; then
  sed -n '2,17p' "$0"; exit 1
fi
SYS=$1 SNAP=$2 Z=$3 HALO=$4 BH1=$5 BH2=$6 BH3=$7 POS=$8
LBOL=None
[ $# -eq 11 ] && LBOL="[$9, ${10}, ${11}]"

D="$RUNS_DIR/sys${SYS}_snap${SNAP}"
if [ -e "$D" ]; then echo "ERROR: $D already exists (move it away first)" >&2; exit 1; fi
mkdir -p "$RUNS_DIR"

# Copy the pipeline (everything except the repo-only folders)
rsync -a --exclude .git --exclude .gitignore --exclude docs --exclude scripts --exclude selection \
      --exclude '__pycache__' "$RTGEN_REPO/" "$D/"

sed -i -e "s|@SYSNUM@|$SYS|g" -e "s|@SNAP@|$SNAP|g" -e "s|@Z@|$Z|g" -e "s|@HALO@|$HALO|g" \
       -e "s|@BH1@|$BH1|g" -e "s|@BH2@|$BH2|g" -e "s|@BH3@|$BH3|g" -e "s|@POS@|$POS|g" \
       -e "s|@LBOL@|$LBOL|g" -e "s|@SIMDIR@|$SIM_DIR|g" -e "s|@SIMPREFIX@|$SIM_PREFIX|g" \
       "$D/input/input_params.py"

printf '# Written by scripts/make_run.sh; sourced by every submit script of this run\n%s\n' \
       "$RTGEN_ENV_SETUP" > "$D/rtgen_env.sh"

if grep -n "@[A-Z0-9]*@" "$D/input/input_params.py"; then
  echo "WARNING: placeholder left in $D/input/input_params.py" >&2
fi
f="$SIM_DIR/$SIM_PREFIX.$SNAP"
[ -e "$f" ] || echo "WARNING: snapshot $f not found" >&2

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$SYS" "$SNAP" "$Z" "$HALO" "$BH1" "$BH2" "$BH3" "$POS" >> "$RUNS_DIR/systems.txt"
echo "Created $D"
