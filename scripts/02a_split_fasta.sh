#!/usr/bin/env bash
set -euo pipefail

ORG="${1:-}"
if [[ -z "$ORG" ]]; then
  echo "Usage: $0 <ORG_CODE>"
  exit 1
fi

IN="inputs/${ORG}/linear/${ORG}_sequences.fasta"
OUTDIR="inputs/${ORG}/linear"

if [[ ! -f "$IN" ]]; then
  echo "ERROR: input FASTA not found: $IN"
  exit 1
fi

mkdir -p "$OUTDIR"

MAX_PER_FILE=50

awk -v max="$MAX_PER_FILE" -v outdir="$OUTDIR" -v org="$ORG" '
BEGIN {
    file_index = 1
    seq_count = 0
    total = 0
    outfile = sprintf("%s/%s_part%d.fasta", outdir, org, file_index)
}

/^>/ {
    seq_count++
    total++

    if (seq_count > max) {
        file_index++
        outfile = sprintf("%s/%s_part%d.fasta", outdir, org, file_index)
        seq_count = 1
    }
}

{
    print >> outfile
}

END {
    print "Total sequences:", total > "/dev/stderr"
    print "Files created:", file_index > "/dev/stderr"
}
' "$IN"

echo ""
echo "Summary:"
for f in "$OUTDIR"/${ORG}_part*.fasta; do
    echo "$(basename "$f"): $(grep -c '^>' "$f") sequences"
done