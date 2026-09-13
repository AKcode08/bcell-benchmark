#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# PURPOSE
#   Combine all UniProt FASTA files for one organism into a
#   single multi-FASTA file for downstream epitope predictors,
#   while enforcing headers as UniProt accession ONLY.
#
# EXPECTED INPUT
#   fasta/<ORG>/*.fasta
#
# OUTPUT
#   inputs/<ORG>/<ORG>_sequences.fasta
#
# USAGE
#   ./scripts/build_multifasta_from_fasta_dir.sh PA
#   ./scripts/build_multifasta_from_fasta_dir.sh SA
# ------------------------------------------------------------

ORG="${1:-}"
if [[ -z "${ORG}" ]]; then
  echo "Usage: $0 <ORG_CODE>"
  exit 1
fi

INDIR="fasta/${ORG}/conformational"
OUTDIR="inputs/${ORG}/conformational"
OUT="${OUTDIR}/${ORG}_sequences.fasta"

if [[ ! -d "${INDIR}" ]]; then
  echo "ERROR: FASTA directory not found: ${INDIR}"
  exit 1
fi

mkdir -p "${OUTDIR}"

N_FASTA="$(ls "${INDIR}"/*.fasta 2>/dev/null | wc -l | tr -d ' ')"
if [[ "${N_FASTA}" -eq 0 ]]; then
  echo "ERROR: No FASTA files found in ${INDIR}"
  exit 1
fi

# Build multi-FASTA with headers sanitized to UniProt accession only.
# UniProt headers typically look like:
#   >sp|P12345|NAME ...   or  >tr|Q9XYZ1|NAME ...
# If not, we fall back to the first whitespace-delimited token (minus '>').

awk '
  BEGIN { OFS=""; }

  function print_record(h, s) {
    if (h == "") return;

    # Remove leading >
    header = substr(h, 2);

    # Try UniProt format: sp|P12345|...
    split(header, parts, "|");

    if (length(parts) >= 2 && (parts[1] == "sp" || parts[1] == "tr")) {
      acc = parts[2];
    } else {
      # fallback: first whitespace token
      split(header, w, /[ \t]/);
      acc = w[1];
    }

    if (!(acc in seen)) {
      seen[acc] = 1;
      print ">", acc;
      print s;
      kept++;
    } else {
      dups++;
    }
  }

  /^>/ {
    print_record(h, seq);
    h = $0;
    seq = "";
    next;
  }

  {
    gsub(/[ \r]/, "", $0);
    if ($0 != "") seq = seq $0;
  }

  END {
    print_record(h, seq);
    print "Kept sequences:", kept > "/dev/stderr";
    print "Duplicate accessions skipped:", dups > "/dev/stderr";
  }
' "${INDIR}"/*.fasta > "${OUT}"

echo "Wrote: ${OUT}"
echo "Number of sequences:"
grep -c "^>" "${OUT}"

echo "Example headers (first 10):"
grep "^>" "${OUT}" | head -n 10