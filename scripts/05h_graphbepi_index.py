# scripts/index_graphbepi.py

from pathlib import Path
import pandas as pd
import numpy as np
import argparse

# =====================================================
# ARGUMENTS
# =====================================================

parser = argparse.ArgumentParser()

parser.add_argument(
    "--organism",
    required=True,
    help="Organism folder name (e.g. ESKAPE, GAST, INTRA, RESP)"
)

args = parser.parse_args()

ORGANISM = args.organism

# =====================================================
# PATHS
# =====================================================

BASE_DIR = Path("outputs/graphbepi") / ORGANISM

GROUP_MAP = {
    "ESKAPE": "bacteria",
    "GAST": "bacteria",
    "INTRA": "bacteria",
    "RESP": "bacteria",
    "COV": "virus",
    "ENZA": "virus"
}

group = GROUP_MAP.get(ORGANISM, "unknown")

# =====================================================
# PROCESS TYPES
# =====================================================

for epitope_type in ["linear", "conformational"]:

    type_dir = BASE_DIR / epitope_type

    if not type_dir.exists():
        print(f"Skipping missing folder: {type_dir}")
        continue

    print("\n===================================")
    print(f"Processing {ORGANISM} | {epitope_type}")
    print("===================================")

    all_dfs = []

    # -------------------------------------------------
    # READ ALL GRAPHBEPI OUTPUTS
    # -------------------------------------------------

    for csv_file in type_dir.glob("*_graphbepi.csv"):

        uniprot = csv_file.stem.replace("_graphbepi", "")

        try:
            df = pd.read_csv(csv_file)

        except Exception as e:
            print(f"Failed reading {csv_file.name}: {e}")
            continue

        print(f"Reading: {csv_file.name}")

        # -------------------------------------------------
        # DETECT COLUMNS
        # -------------------------------------------------

        score_col = None

        for col in df.columns:
            if "score" in col.lower():
                score_col = col
                break

        if score_col is None:
            print(f"No score column found in {csv_file.name}")
            continue

        residue_col = "resn"
        pred_col = "is epitope"

        # -------------------------------------------------
        # BINARY PREDICTIONS
        # -------------------------------------------------

        y_pred = (
            df[pred_col]
            .astype(str)
            .str.lower()
            .isin(["true", "1", "yes", "epitope"])
        ).astype(int)

        # -------------------------------------------------
        # STANDARDIZE
        # -------------------------------------------------

        out = pd.DataFrame({

            "model": "GraphBepi",

            "organism": ORGANISM,

            "group": group,

            "type": epitope_type,

            "uniprot": uniprot,

            "chain": "A",

            "res_index": range(1, len(df) + 1),

            "residue": df[residue_col],

            "raw_score": df[score_col],

            "calibrated_score": df[score_col],

            "score": df[score_col],

            "threshold": np.nan,

            "y_pred": y_pred,

            # -------------------------------------------------
            # ELLIPRO-STYLE STANDARDIZATION
            # -------------------------------------------------

            "graphbepi_score": df[score_col],

            "assignment": np.where(y_pred == 1, "E", "n")
        })

        all_dfs.append(out)

    # =====================================================
    # COMBINE + SAVE
    # =====================================================

    if len(all_dfs) == 0:

        print(f"No valid files found for {ORGANISM} {epitope_type}")
        continue

    combined = pd.concat(all_dfs, ignore_index=True)

    outfile = type_dir / "graphbepi_indexed.csv"

    combined.to_csv(outfile, index=False)

    print("\n===================================")
    print("Saved:")
    print(outfile)

    print("\nShape:")
    print(combined.shape)

    print("\nProteins:")
    print(combined['uniprot'].nunique())

    print("===================================")