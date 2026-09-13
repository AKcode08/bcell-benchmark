#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

BASE_DIR = Path("outputs/ellipro")

MODEL = "ElliPro"
THRESHOLD = 0.50

VIRUS_ORGS = {"COV", "ENZA"}


def get_group(org: str) -> str:
    return "virus" if org in VIRUS_ORGS else "bacteria"


def process_folder(input_dir: Path, org: str, dataset_type: str):
    output_file = input_dir / "ellipro_indexed.csv"

    # avoid re-reading combined output if rerun
    all_files = [
        f for f in input_dir.glob("*.csv")
        if f.name != "ellipro_indexed.csv"
    ]

    if not all_files:
        print(f"Skipping {input_dir}: no ElliPro CSV files found.")
        return

    group = get_group(org)
    dfs = []

    for f in all_files:
        print(f"  Processing: {f.name}")

        df = pd.read_csv(f)

        required_cols = {
            "organism",
            "type",
            "uniprot",
            "chain",
            "res_index",
            "aa",
            "ellipro_score",
            "assignment",
            "y_pred",
        }

        if not required_cols.issubset(df.columns):
            missing = sorted(required_cols - set(df.columns))
            raise ValueError(f"{f} missing required columns: {missing}")

        score = df["ellipro_score"].astype(float)

        tmp = pd.DataFrame({
            # standardized core columns (same schema as BepiPred / DiscoTope)
            "model": MODEL,
            "organism": org,
            "group": group,
            "type": dataset_type,
            "uniprot": df["uniprot"],
            "chain": df["chain"],
            "res_index": df["res_index"].astype(int),
            "residue": df["aa"],
            "raw_score": score,
            "calibrated_score": score,
            "score": score,
            "threshold": THRESHOLD,
            "y_pred": df["y_pred"].astype(int),

            # ElliPro-specific retained columns
            "ellipro_score": score,
            "assignment": df["assignment"],
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