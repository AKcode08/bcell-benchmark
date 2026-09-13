#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Sync the BlueBEAR project tree into this local clone.
#
#   MERGE SEMANTICS: additive only. Nothing is deleted on either side, and
#   nothing that already exists locally is overwritten (--ignore-existing).
#   Files that exist on both sides but differ are REPORTED, never auto-resolved.
#
#   Dry-run is the default. Nothing is transferred until you pass --apply.
#
#   Usage:
#     ./sync_bluebear.sh                 # phases 1-3, dry run (safe)
#     ./sync_bluebear.sh --apply         # actually transfer
#     ./sync_bluebear.sh --apply --with-weights   # include model checkpoints (~13 GB)
#     ./sync_bluebear.sh --push          # additionally send local-only files UP
# ---------------------------------------------------------------------------
set -euo pipefail

REMOTE_USER="${REMOTE_USER:-axk1656}"
REMOTE_HOST="${REMOTE_HOST:-bluebear.bham.ac.uk}"
REMOTE_PATH="${REMOTE_PATH:-/rds/projects/e/elhamsak-epitope-dev/bcell_benchmark}"
LOCAL_PATH="${LOCAL_PATH:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

REMOTE="${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_PATH}/"

DRY="--dry-run"; WEIGHTS=0; PUSH=0
for arg in "$@"; do
  case "$arg" in
    --apply)        DRY="" ;;
    --with-weights) WEIGHTS=1 ;;
    --push)         PUSH=1 ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

# Never sync these in either direction.
EXCLUDES=(
  --exclude 'venv/'
  --exclude '.DS_Store'
  --exclude '__pycache__/'
  --exclude '.ipynb_checkpoints/'
  --exclude '.git/'
  --exclude '*.egg-info/'
)

# Model weights are large and re-downloadable from upstream; opt in explicitly.
WEIGHT_EXCLUDES=(
  --exclude '*.pt'  --exclude '*.pth'  --exclude '*.ckpt'
  --exclude '*.bin' --exclude '*.safetensors' --exclude '*.h5'
  --exclude 'hf_cache/' --exclude 'data/embeddings/'
)
if [[ $WEIGHTS -eq 0 ]]; then
  EXCLUDES+=("${WEIGHT_EXCLUDES[@]}")
fi

banner() { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }
[[ -n "$DRY" ]] && printf '\n\033[33mDRY RUN — no files will be transferred. Re-run with --apply.\033[0m\n'

banner "Phase 0: remote reachable?"
ssh -o BatchMode=yes -o ConnectTimeout=10 "${REMOTE_USER}@${REMOTE_HOST}" \
    "test -d '${REMOTE_PATH}' && du -sh '${REMOTE_PATH}' 2>/dev/null" \
  || { echo "Cannot reach ${REMOTE_HOST} or path missing. Check VPN / SSH key."; exit 1; }

banner "Phase 1: what exists on BlueBEAR but not here (would be pulled DOWN)"
rsync -azi --dry-run --ignore-existing "${EXCLUDES[@]}" "$REMOTE" "$LOCAL_PATH/" \
  | grep -E '^[<>]' | awk '{print "  + " $2}' | head -60
echo "  ..."
rsync -az --dry-run --stats --ignore-existing "${EXCLUDES[@]}" "$REMOTE" "$LOCAL_PATH/" \
  | grep -E 'Number of regular files transferred|Total transferred file size'

banner "Phase 2: CONFLICTS — exist on both sides but differ (NOT touched by this script)"
rsync -azi --dry-run --existing "${EXCLUDES[@]}" "$REMOTE" "$LOCAL_PATH/" \
  | grep -E '^>f.*[cst]' | awk '{print "  ! " $2}' | head -40 \
  || echo "  (none)"
echo "  ^ resolve these by hand; nothing above was modified."

banner "Phase 3: pulling BlueBEAR -> local (additive only)"
rsync -az --info=progress2 --human-readable $DRY --ignore-existing \
      "${EXCLUDES[@]}" "$REMOTE" "$LOCAL_PATH/"

if [[ $PUSH -eq 1 ]]; then
  banner "Phase 4: pushing local-only files -> BlueBEAR (additive only)"
  rsync -az --info=progress2 --human-readable $DRY --ignore-existing \
        "${EXCLUDES[@]}" "$LOCAL_PATH/" "$REMOTE"
fi

banner "Done"
[[ -n "$DRY" ]] && echo "That was a dry run. Re-run with --apply to transfer." || du -sh "$LOCAL_PATH"
