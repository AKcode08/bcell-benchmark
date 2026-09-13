#!/usr/bin/env python3

import re
import sys
import csv
from pathlib import Path

UNIPROT_PIPE_RE = re.compile(r"(?:^|[>,\s])(sp|tr)\|([A-Za-z0-9]+)\|")
UNIPROT_BARE_RE = re.compile(r"(?:^|[>,\s])([A-Za-z0-9]{6,10})(?:\s|$)")

MODEL = "BepiPred"
EPITOPE_TYPE = "conformational"
THRESHOLD = 0.151

VIRUS_ORGS = {"COV", "ENZA"}


def usage():
    print(
        "Usage: python scripts/index_bepipred3_output_conformational.py <ORG_CODE>",
        file=sys.stderr
    )
    sys.exit(1)


def extract_accession(seq_id):
    m = UNIPROT_PIPE_RE.search(seq_id)
    if m:
        return m.group(2).upper()

    m2 = UNIPROT_BARE_RE.search(seq_id.strip())
    if m2:
        return m2.group(1).upper()

    return None


def pt_sort_key(path):
    m = re.search(r"PT(\d+)raw_output\.csv", path.name)
    return int(m.group(1)) if m else 10**9


def combine_split_outputs(parts, raw):
    parts = sorted(parts, key=pt_sort_key)

    with raw.open("w", encoding="utf-8", newline="") as out_f:
        for i, part in enumerate(parts):
            lines = part.read_text(encoding="utf-8", errors="replace").splitlines(True)

            if not lines:
                raise SystemExit(f"ERROR: empty file: {part}")

            if i == 0:
                out_f.writelines(lines)
            else:
                out_f.writelines(lines[1:])

    print(f"Combined {len(parts)} split outputs into: {raw}")


def main():
    if len(sys.argv) != 2:
        usage()

    org = sys.argv[1].strip()
    group = "virus" if org in VIRUS_ORGS else "bacteria"

    outdir = Path(f"outputs/bepipred3/{org}/conformational")
    raw = outdir / "raw_output.csv"
    out = outdir / "bepipred3_indexed.csv"

    split_parts = sorted(outdir.glob("PT*raw_output.csv"), key=pt_sort_key)

    if split_parts:
        combine_split_outputs(split_parts, raw)

    if not raw.exists():
        raise SystemExit(f"ERROR: input not found: {raw}")

    out.parent.mkdir(parents=True, exist_ok=True)

    with raw.open("r", encoding="utf-8", errors="replace") as f, out.open(
        "w", newline="", encoding="utf-8"
    ) as o:

        writer = csv.writer(o)

        writer.writerow([
            "model",
            "organism",
            "group",
            "type",
            "uniprot",
            "chain",
            "res_index",
            "residue",
            "raw_score",
            "calibrated_score",
            "score",
            "threshold",
            "y_pred"
        ])

        current_uniprot = None
        res_index = 0
        rows_written = 0
        rows_skipped = 0

        header = f.readline()

        if not header:
            raise SystemExit(f"ERROR: empty file: {raw}")

        for line in f:
            line = line.strip()

            if not line:
                continue

            parts = line.rsplit(",", 3)

            if len(parts) != 4:
                rows_skipped += 1
                continue

            seq_id, aa, raw_score, cal_score = parts
            uniprot = extract_accession(seq_id)

            if not uniprot:
                rows_skipped += 1
                continue

            if uniprot != current_uniprot:
                current_uniprot = uniprot
                res_index = 1
            else:
                res_index += 1

            try:
                raw_score = float(raw_score)
                cal_score = float(cal_score)
                score = raw_score
                y_pred = int(score >= THRESHOLD)

                writer.writerow([
                    MODEL,
                    org,
                    group,
                    EPITOPE_TYPE,
                    uniprot,
                    "",
                    res_index,
                    aa,
                    raw_score,
                    cal_score,
                    score,
                    THRESHOLD,
                    y_pred
                ])

                rows_written += 1

            except ValueError:
                rows_skipped += 1
                continue

    print(f"Input:  {raw}")
    print(f"Output: {out}")
    print(f"Rows written: {rows_written}")
    print(f"Rows skipped: {rows_skipped}")


if __name__ == "__main__":
    main()