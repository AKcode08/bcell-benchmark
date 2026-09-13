#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# PURPOSE
#   Download UniProt FASTA sequences referenced in an IEDB CSV
#   Export for ONE organism.

# EXPECTED INPUT STRUCTURE
#   iedb/<ORG>/<ORG>_raw.csv

# USAGE (Example)
#   ./scripts/01_fasta_download.sh ESKAPE_linear
#   ./scripts/01_fasta_download.sh ESKAPE_conformational

# OUTPUT
#   fasta/<ORG>/<ACC>.fasta
# ------------------------------------------------------------

ORG="${1:-}"
RAW="iedb/${ORG}/conformational/conformational_raw.csv"
OUTDIR="fasta/${ORG}/conformational/"

if [[ -z "${ORG}" ]]; then
  echo "Usage: $0 <ORG_CODE>"
  exit 1
fi

if [[ ! -f "${RAW}" ]]; then
  echo "ERROR: IEDB CSV not found: ${RAW}"
  exit 1
fi

mkdir -p "${OUTDIR}"

# Extract UniProt accessions from UniProt IRIs
ACCESSIONS=$(grep -Eo 'https?://(www\.)?uniprot\.org/(uniprot|uniprotkb)/[A-Z0-9]+' "${RAW}" \
  | sed -E 's|.*/||' \
  | sort -u)

if [[ -z "${ACCESSIONS}" ]]; then
  echo "ERROR: No UniProt accessions found"
  exit 1
fi

echo "Found $(echo "${ACCESSIONS}" | wc -w | tr -d ' ') UniProt accessions for ${ORG}"

FAILED=0

for acc in ${ACCESSIONS}; do
  out="${OUTDIR}/${acc}.fasta"

  if [[ -s "${out}" ]]; then
    echo "Skip (exists): ${acc}"
    continue
  fi

  echo "Downloading ${acc}"

  if ! curl -fsSL --retry 3 --retry-delay 1 \
      "https://rest.uniprot.org/uniprotkb/${acc}.fasta" \
      -o "${out}"; then
    echo "FAILED: ${acc}"
    FAILED=$((FAILED + 1))
    rm -f "${out}" || true
    continue
  fi

  # Sanity check: exactly one FASTA header
  nhead="$(grep -c '^>' "${out}" || true)"
  if [[ "${nhead}" -ne 1 ]]; then
    echo "WARNING: ${acc} FASTA has ${nhead} headers"
  fi

  sleep 0.2
done

echo "Done. FASTAs written to ${OUTDIR}"
grep -c "^>" "${OUTDIR}"/*.fasta || true

if [[ "${FAILED}" -gt 0 ]]; then
  echo "WARNING: ${FAILED} accessions failed"
  exit 2
fi