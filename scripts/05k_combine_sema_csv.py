#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

BASE_DIR = Path("outputs/SEMAi")

MODEL = "SEMA"
THRESHOLD = 0.51

VIRUS_ORGS = {"COV", "ENZA"}


def get_group(org: str) -> str:
    return "virus" if org in VIRUS_ORGS else "bacteria"


def process_folder(input_dir: Path, org: str, dataset_type: str):
    output_file = input_dir / "sema_indexed.csv"

    # avoid re-reading combined output if rerun
    all_files = [
        f for f in input_dir.glob("*.csv")
        if f.name != "sema_indexed.csv"
    ]

    if not all_files:
        print(f"Skipping {input_dir}: no SEMA CSV files found.")
        return

    group = get_group(org)
    dfs = []

    for f in sorted(all_files):
        print(f"  Processing: {f.name}")

        df = pd.read_csv(f)

        # SEMA per-antigen files contain exactly: residue, aa, score
        required_cols = {"residue", "aa", "score"}
        if not required_cols.issubset(df.columns):
            missing = sorted(required_cols - set(df.columns))
            raise ValueError(f"{f} missing required columns: {missing}")

        # uniprot is derived from the filename stem (e.g. G3XD49.csv -> G3XD49)
        uniprot = f.stem

        score = df["score"].astype(float)

        # epitope decision (inclusive threshold)
        y_pred = (score >= THRESHOLD).astype(int)
        assignment = y_pred.map({1: "E", 0: "e"})

        tmp = pd.DataFrame({
            # standardized core columns (same schema as BepiPred / DiscoTope / ElliPro)
            "model": MODEL,
            "organism": org,
            "group": group,
            "type": dataset_type,
            "uniprot": uniprot,
            "chain": "A",
            "res_index": df["residue"].astype(int),
            "residue": df["aa"],
            "raw_score": score,
            "calibrated_score": score,
            "score": score,
            "threshold": THRESHOLD,
            "y_pred": y_pred,

            # SEMA-specific retained columns
            "sema_score": score,
            "assignment": assignment,
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
