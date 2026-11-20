# Sourced by the other scripts: locate the repo and load config.sh
RTGEN_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ ! -f "$RTGEN_REPO/scripts/config.sh" ]; then
  echo "ERROR: $RTGEN_REPO/scripts/config.sh not found. Run: cp scripts/config.example.sh scripts/config.sh  and edit it." >&2
  exit 1
fi
source "$RTGEN_REPO/scripts/config.sh"
