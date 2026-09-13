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
#   TWO SOURCE MODES:
#     SSH   (default) -- pulls over ssh from a BlueBEAR login node.
#     LOCAL (SOURCE_DIR=...) -- reads the RDS project directory mounted on this
#           machine over SMB. No ssh, no keys, no login prompt. Mount it in
#           Finder with Cmd-K; the exact smb:// path for your project is in BEAR
#           Admin (BEAR Identity Database) or the email sent when the project
#           storage was created. Log in as ADF\\<username> or
#           adf.bham.ac.uk/<username>. Then:
#             SOURCE_DIR="/Volumes/<project>/bcell_benchmark" ./sync_bluebear.sh
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

# SOURCE_DIR set => read from a locally mounted RDS share instead of over ssh.
SOURCE_DIR="${SOURCE_DIR:-}"
if [[ -n "$SOURCE_DIR" ]]; then
  LOCAL_MODE=1
  REMOTE="${SOURCE_DIR%/}/"
else
  LOCAL_MODE=0
  REMOTE="${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_PATH}/"
fi

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

# One authenticated connection reused by every phase, so password/2FA is
# entered once rather than per rsync invocation.
SSH_CTL="${TMPDIR:-/tmp}/bbsync-%r@%h-%p"
SSH_OPTS=(-o ControlMaster=auto -o "ControlPath=${SSH_CTL}" -o ControlPersist=10m -o ConnectTimeout=20)
SSH_CMD="ssh -o ControlMaster=auto -o ControlPath=${SSH_CTL} -o ControlPersist=10m -o ConnectTimeout=20"

banner() { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }
[[ -n "$DRY" ]] && printf '\n\033[33mDRY RUN — no files will be transferred. Re-run with --apply.\033[0m\n'

if [[ $LOCAL_MODE -eq 1 ]]; then
  banner "Phase 0: mounted source"
  if [[ ! -d "$SOURCE_DIR" ]]; then
    echo "Not a directory: $SOURCE_DIR"
    echo "Mount the RDS share first (Finder > Go > Connect to Server, Cmd-K)."
    echo "The smb:// path for your project is in BEAR Admin, or the project-creation email."
    exit 1
  fi
  echo "  source: $SOURCE_DIR"
  du -sh "$SOURCE_DIR" 2>/dev/null | sed 's/^/  /' || true
  SSH_CMD=""            # rsync runs as a plain local copy
else
banner "Phase 0: authenticate (you may be prompted once)"
if ! ssh "${SSH_OPTS[@]}" "${REMOTE_USER}@${REMOTE_HOST}" \
        "test -d '${REMOTE_PATH}' && du -sh '${REMOTE_PATH}' 2>/dev/null"; then
  cat <<DIAG

Could not authenticate to ${REMOTE_USER}@${REMOTE_HOST}, or ${REMOTE_PATH} is missing.

  "Permission denied (publickey,...)" means the host answered and rejected your
  credentials -- the network is fine, the login is not. Check, in order:

    1. ssh ${REMOTE_USER}@${REMOTE_HOST}          # does a plain login work?
    2. grep -iA5 bluebear ~/.ssh/config           # is there a Host alias with a
                                                  # different user or IdentityFile?
    3. ssh-add -l                                 # is your key loaded in the agent?
    4. Is the username right? This script assumes REMOTE_USER=${REMOTE_USER}.

  Override any of them without editing this file:
    REMOTE_USER=xxx REMOTE_HOST=yyy ./sync_bluebear.sh

DIAG
  exit 1
fi
fi

banner "Phase 1: what exists on BlueBEAR but not here (would be pulled DOWN)"
rsync -azi ${SSH_CMD:+-e "$SSH_CMD"} --dry-run --ignore-existing "${EXCLUDES[@]}" "$REMOTE" "$LOCAL_PATH/" \
  | grep -E '^[<>]' | awk '{print "  + " $2}' | head -60
echo "  ..."
rsync -az ${SSH_CMD:+-e "$SSH_CMD"} --dry-run --stats --ignore-existing "${EXCLUDES[@]}" "$REMOTE" "$LOCAL_PATH/" \
  | grep -E 'Number of regular files transferred|Total transferred file size'

banner "Phase 2: CONFLICTS — exist on both sides but differ (NOT touched by this script)"
rsync -azi ${SSH_CMD:+-e "$SSH_CMD"} --dry-run --existing "${EXCLUDES[@]}" "$REMOTE" "$LOCAL_PATH/" \
  | grep -E '^>f.*[cst]' | awk '{print "  ! " $2}' | head -40 \
  || echo "  (none)"
echo "  ^ resolve these by hand; nothing above was modified."

banner "Phase 3: pulling BlueBEAR -> local (additive only)"
rsync -az ${SSH_CMD:+-e "$SSH_CMD"} --info=progress2 --human-readable $DRY --ignore-existing \
      "${EXCLUDES[@]}" "$REMOTE" "$LOCAL_PATH/"

if [[ $PUSH -eq 1 ]]; then
  banner "Phase 4: pushing local-only files -> BlueBEAR (additive only)"
  rsync -az ${SSH_CMD:+-e "$SSH_CMD"} --info=progress2 --human-readable $DRY --ignore-existing \
        "${EXCLUDES[@]}" "$LOCAL_PATH/" "$REMOTE"
fi

[[ $LOCAL_MODE -eq 0 ]] && { ssh "${SSH_OPTS[@]}" -O exit "${REMOTE_USER}@${REMOTE_HOST}" 2>/dev/null || true; }

banner "Done"
[[ -n "$DRY" ]] && echo "That was a dry run. Re-run with --apply to transfer." || du -sh "$LOCAL_PATH"
