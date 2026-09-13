#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

BASE_DIR = Path("outputs/discotope3")

MODEL = "DiscoTope"
THRESHOLD = 0.90

VIRUS_ORGS = {"COV", "ENZA"}


def get_group(org: str) -> str:
    return "virus" if org in VIRUS_ORGS else "bacteria"


def extract_uniprot(pdb_field):
    # Example: A0A0H3GIG3_A -> A0A0H3GIG3
    return str(pdb_field).split("_")[0]


def process_folder(input_dir: Path, org: str, dataset_type: str):
    output_file = input_dir / "discotope3_indexed.csv"

    all_files = [
        f for f in input_dir.glob("*.csv")
        if f.name != "discotope3_indexed.csv"
    ]

    if not all_files:
        print(f"Skipping {input_dir}: no DiscoTope CSV files found.")
        return

    group = get_group(org)
    dfs = []

    for f in all_files:
        print(f"  Processing: {f.name}")

        df = pd.read_csv(f)

        required_cols = {
            "pdb",
            "chain",
            "res_id",
            "residue",
            "DiscoTope-3.0_score",
            "calibrated_score",
            "epitope",
            "rsa",
            "pLDDTs",
            "length",
        }

        if not required_cols.issubset(df.columns):
            missing = sorted(required_cols - set(df.columns))
            raise ValueError(f"{f} missing required columns: {missing}")

        score = df["calibrated_score"].astype(float)
        y_pred = (score >= THRESHOLD).astype(int)

        tmp = pd.DataFrame({
            # standardized core columns (same as BepiPred)
            "model": MODEL,
            "organism": org,
            "group": group,
            "type": dataset_type,
            "uniprot": df["pdb"].map(extract_uniprot),
            "chain": df["chain"],
            "res_index": df["res_id"].astype(int),
            "residue": df["residue"],
            "raw_score": df["DiscoTope-3.0_score"].astype(float),
            "calibrated_score": df["calibrated_score"].astype(float),
            "score": score,
            "threshold": THRESHOLD,
            "y_pred": y_pred,

            # model-specific extra columns retained
            "epitope": df["epitope"],
            "rsa": df["rsa"],
            "pLDDTs": df["pLDDTs"],
            "length": df["length"],
            "pdb": df["pdb"],
        })

        dfs.append(tmp)

    combined = pd.concat(dfs, ignore_index=True)

    combined = combined.sort_values(
        ["uniprot", "res_index"]
    ).reset_index(drop=True)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_file, index=False)

    print(f"  Wrote: {output_file}")
    print(f"  Total rows: {len(combined)}")
    print(f"  Proteins: {combined['uniprot'].nunique()}")


def main():
    if not BASE_DIR.exists():
        raise SystemExit(f"Base directory not found: {BASE_DIR}")

    for org_dir in sorted([d for d in BASE_DIR.iterdir() if d.is_dir()]):

        org = org_dir.name

        print(f"\nOrganism: {org}")

        for dataset_type in ["linear", "conformational"]:

            input_dir = org_dir / dataset_type

            if not input_dir.exists():
                print(f"Skipping {input_dir}: folder not found.")
                continue

            print(f" Dataset: {dataset_type}")

            process_folder(input_dir, org, dataset_type)


if __name__ == "__main__":
    main()