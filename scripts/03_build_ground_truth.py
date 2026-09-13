#!/usr/bin/env python3

import csv
import re
import sys
from pathlib import Path

UNIPROT_RE = re.compile(r"(?:uniprot(?:kb)?/)([A-Z0-9]+)", re.IGNORECASE)
CHEBI_RE = re.compile(r"\bCHEBI[_:/-]?\d+\b", re.IGNORECASE)

COL_CANDIDATES = {
    "iedb_iri": ["epitope id - iedb iri", "epitope iri", "iedb iri"],
    "object_type": ["epitope - object type", "object type", "epitope object type", "epitope type"],
    "sequence": ["epitope - name", "epitope sequence", "epitope", "sequence"],
    "start": ["epitope - starting position", "start position", "start", "epitope start"],
    "end": ["epitope - ending position", "end position", "end", "epitope end"],
    "organism": ["epitope - species", "species", "epitope species"],
}

def norm(s: str) -> str:
    return re.sub(r"\s+", " ", s.strip().lower())

def find_col(header: list[str], candidates: list[str]) -> int | None:
    h = [norm(x) for x in header]
    for c in candidates:
        c_norm = norm(c)
        if c_norm in h:
            return h.index(c_norm)
    return None

def extract_uniprot_from_row(row: list[str]) -> str | None:
    for cell in row:
        if not cell:
            continue
        m = UNIPROT_RE.search(cell)
        if m:
            return m.group(1).upper()
    return None

def row_has_chebi(row: list[str]) -> bool:
    for cell in row:
        if cell and CHEBI_RE.search(cell):
            return True
    return False

def sniff_dialect(path: Path) -> csv.Dialect:
    sample = path.read_text(encoding="utf-8", errors="replace")[:20000]
    try:
        return csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";"])
    except csv.Error:
        class Tsv(csv.Dialect):
            delimiter = "\t"
            quotechar = '"'
            doublequote = True
            skipinitialspace = False
            lineterminator = "\n"
            quoting = csv.QUOTE_MINIMAL
        return Tsv()

def usage() -> None:
    print("Usage: ./scripts/make_ground_truth_from_iedb.py <ORG_CODE>")
    print("Example: ./scripts/make_ground_truth_from_iedb.py SA")
    sys.exit(1)

def main() -> None:
    if len(sys.argv) != 2:
        usage()

    org = sys.argv[1].strip()
    raw = Path(f"iedb/{org}/linear/linear_raw.csv")
    out = Path(f"iedb/{org}/linear/linear_ground_truth.csv")

    if not raw.exists():
        print(f"ERROR: input not found: {raw}", file=sys.stderr)
        sys.exit(2)

    dialect = sniff_dialect(raw)

    drops = {
        "missing_uniprot": 0,
        "missing_start_end": 0,
        "non_protein_antigen_chebi": 0,
        "non_linear_peptide": 0,
        "bad_coords": 0,
        "missing_iedb_iri": 0,
        "other_missing_required": 0,
    }

    with raw.open(newline="", encoding="utf-8") as f:
        reader = csv.reader(f, dialect=dialect)
        try:
            header = next(reader)
        except StopIteration:
            print(f"ERROR: empty file: {raw}", file=sys.stderr)
            sys.exit(3)

        idx_iri = find_col(header, COL_CANDIDATES["iedb_iri"])
        idx_type = find_col(header, COL_CANDIDATES["object_type"])
        idx_seq = find_col(header, COL_CANDIDATES["sequence"])
        idx_start = find_col(header, COL_CANDIDATES["start"])
        idx_end = find_col(header, COL_CANDIDATES["end"])
        idx_org = find_col(header, COL_CANDIDATES["organism"])

        missing = [k for k, v in [("sequence", idx_seq), ("start", idx_start), ("end", idx_end)] if v is None]
        if missing:
            print("ERROR: Could not locate required columns:", ", ".join(missing), file=sys.stderr)
            print("Detected header:", header, file=sys.stderr)
            sys.exit(4)

        out.parent.mkdir(parents=True, exist_ok=True)

        rows_written = 0
        with out.open("w", newline="", encoding="utf-8") as fo:
            w = csv.writer(fo)
            w.writerow([
                "iedb_id",
                "uniprot_accession",
                "organism",
                "epitope_sequence",
                "start_pos",
                "end_pos",
                "epitope_len"
            ])

            for row in reader:
                if not row:
                    continue

                obj_type = row[idx_type].strip() if idx_type is not None and idx_type < len(row) else ""
                if obj_type and obj_type.lower() != "linear peptide":
                    drops["non_linear_peptide"] += 1
                    continue

                epitope_seq = row[idx_seq].strip()
                start = row[idx_start].strip()
                end = row[idx_end].strip()
                organism = row[idx_org].strip() if idx_org is not None and idx_org < len(row) else ""

                uniprot = extract_uniprot_from_row(row)

                if not uniprot:
                    drops["missing_uniprot"] += 1
                    if row_has_chebi(row):
                        drops["non_protein_antigen_chebi"] += 1

                if not (start and end):
                    drops["missing_start_end"] += 1

                if not (epitope_seq and start and end and uniprot):
                    drops["other_missing_required"] += 1
                    continue

                try:
                    start_i = int(float(start))
                    end_i = int(float(end))
                except ValueError:
                    drops["bad_coords"] += 1
                    continue

                if idx_iri is not None and idx_iri < len(row) and row[idx_iri].strip():
                    iri = row[idx_iri].strip()
                    iedb_id = iri.rstrip("/").split("/")[-1]
                else:
                    drops["missing_iedb_iri"] += 1
                    continue

                w.writerow([
                    iedb_id,
                    uniprot,
                    organism,
                    epitope_seq,
                    start_i,
                    end_i,
                    len(epitope_seq)
                ])
                rows_written += 1

    print(f"Detected delimiter: {repr(dialect.delimiter)}")
    print(f"Input:  {raw}")
    print(f"Output: {out}")
    print(f"Rows written: {rows_written}")
    print("Drop counts:")
    for k, v in drops.items():
        print(f"  {k}: {v}")

if __name__ == "__main__":
    main()