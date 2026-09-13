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
    "residues": ["epitope - name", "epitope residues", "residues"],
    "organism": ["epitope - species", "species", "epitope species"],
}

def norm(s: str) -> str:
    return re.sub(r"\s+", " ", s.strip().lower())

def find_col(header: list[str], candidates: list[str]) -> int | None:
    h = [norm(x) for x in header]
    for c in candidates:
        if norm(c) in h:
            return h.index(norm(c))
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

def parse_residue_set(text: str) -> set[int]:
    """
    Parse:
      D109
      K101, E102, R104
    → {101, 102, 104}
    """
    if not text:
        return set()

    nums = re.findall(r"\d+", text)
    return {int(x) for x in nums}

def usage() -> None:
    print("Usage: ./scripts/03a_build_ground_truth_disc.py <ORG_CODE>")
    print("Example: ./scripts/03a_build_ground_truth_disc.py SA")
    sys.exit(1)

def main() -> None:
    if len(sys.argv) != 2:
        usage()

    org = sys.argv[1].strip()

    # PATHS (same structure as your linear script)
    raw = Path(f"iedb/{org}/conformational/conformational_raw.csv")
    out = Path(f"iedb/{org}/conformational/conformational_ground_truth.csv")

    if not raw.exists():
        print(f"ERROR: input not found: {raw}", file=sys.stderr)
        sys.exit(2)

    dialect = sniff_dialect(raw)

    drops = {
        "missing_uniprot": 0,
        "missing_residues": 0,
        "non_protein_antigen_chebi": 0,
        "non_conformational": 0,
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
        idx_res = find_col(header, COL_CANDIDATES["residues"])
        idx_org = find_col(header, COL_CANDIDATES["organism"])

        if idx_res is None:
            print("ERROR: Could not locate residue column", file=sys.stderr)
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
                "residue_string",
                "residue_set",
                "n_residues"
            ])

            for row in reader:
                if not row:
                    continue

                obj_type = row[idx_type].strip() if idx_type is not None else ""

                # keep only discontinuous epitopes
                if obj_type and "discontinuous" not in obj_type.lower():
                    drops["non_conformational"] += 1
                    continue

                residue_text = row[idx_res].strip()
                organism = row[idx_org].strip() if idx_org is not None else ""

                uniprot = extract_uniprot_from_row(row)

                if not uniprot:
                    drops["missing_uniprot"] += 1
                    if row_has_chebi(row):
                        drops["non_protein_antigen_chebi"] += 1

                residue_set = parse_residue_set(residue_text)

                if not residue_set:
                    drops["missing_residues"] += 1

                if not (uniprot and residue_set):
                    drops["other_missing_required"] += 1
                    continue

                if idx_iri is not None and row[idx_iri].strip():
                    iri = row[idx_iri].strip()
                    iedb_id = iri.rstrip("/").split("/")[-1]
                else:
                    drops["missing_iedb_iri"] += 1
                    continue

                w.writerow([
                    iedb_id,
                    uniprot,
                    organism,
                    residue_text,
                    sorted(residue_set),
                    len(residue_set)
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