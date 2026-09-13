#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd

MODEL = "BepiPred"
EPITOPE_TYPE = "linear"

THRESH = 0.151
K_ABS = 5
P_REL = 0.25

VIRUS_ORGS = {"COV", "ENZA"}


def usage() -> None:
    print(
        "Usage: python scripts/05_bepipred_benchmark.py <ORG_CODE>",
        file=sys.stderr,
    )
    sys.exit(1)


def get_group(org: str) -> str:
    return "virus" if org in VIRUS_ORGS else "bacteria"


def main() -> None:
    if len(sys.argv) != 2:
        usage()

    org = sys.argv[1].strip()
    group = get_group(org)

    pred_path = Path(f"outputs/bepipred3/{org}/linear/bepipred3_indexed.csv")
    gt_path = Path(f"iedb/{org}/linear/linear_ground_truth.csv")
    out_path = Path(f"tables/{org}/linear/{org}_bepipred3_results.csv")

    if not pred_path.exists():
        raise SystemExit(f"ERROR: predictions not found: {pred_path}")

    if not gt_path.exists():
        raise SystemExit(f"ERROR: ground truth not found: {gt_path}")

    pred = pd.read_csv(pred_path)
    gt = pd.read_csv(gt_path)

    required_pred_cols = {"uniprot", "res_index", "calibrated_score"}
    required_gt_cols = {
        "iedb_id",
        "uniprot_accession",
        "organism",
        "start_pos",
        "end_pos",
        "epitope_len",
    }

    if not required_pred_cols.issubset(pred.columns):
        missing = sorted(required_pred_cols - set(pred.columns))
        raise SystemExit(f"ERROR: pred file missing columns: {missing}")

    if not required_gt_cols.issubset(gt.columns):
        missing = sorted(required_gt_cols - set(gt.columns))
        raise SystemExit(f"ERROR: ground truth file missing columns: {missing}")

    pred["res_index"] = pred["res_index"].astype(int)

    pred_max = pred.groupby("uniprot")["res_index"].max().to_dict()

    results = []
    missing_pred_protein = 0
    out_of_range_ct = 0

    for _, r in gt.iterrows():
        organism_name = str(r["organism"]).strip()
        uniprot = str(r["uniprot_accession"]).strip()
        start = int(r["start_pos"])
        end = int(r["end_pos"])
        ep_len = int(r["epitope_len"])

        max_i = pred_max.get(uniprot)

        if max_i is None:
            missing_pred_protein += 1

            results.append({
                "model": MODEL,
                "organism": org,
                "organism_name": organism_name,
                "group": group,
                "type": EPITOPE_TYPE,
                "iedb_id": r["iedb_id"],
                "uniprot": uniprot,
                "start": start,
                "end": end,
                "epitope_len": ep_len,
                "predicted_hits": None,
                "recall": None,
                "hit_abs": None,
                "hit_rel": None,
                "hit": None,
                "threshold": THRESH,
                "K_abs": K_ABS,
                "P_rel": P_REL,
                "out_of_range": True,
                "reason": "missing_prediction_for_protein",
            })
            continue

        if end > int(max_i):
            out_of_range_ct += 1

            results.append({
                "model": MODEL,
                "organism": org,
                "organism_name": organism_name,
                "group": group,
                "type": EPITOPE_TYPE,
                "iedb_id": r["iedb_id"],
                "uniprot": uniprot,
                "start": start,
                "end": end,
                "epitope_len": ep_len,
                "predicted_hits": None,
                "recall": None,
                "hit_abs": None,
                "hit_rel": None,
                "hit": None,
                "threshold": THRESH,
                "K_abs": K_ABS,
                "P_rel": P_REL,
                "out_of_range": True,
                "reason": f"end_pos_exceeds_pred_len({end}>{int(max_i)})",
            })
            continue

        sub = pred[
            (pred["uniprot"] == uniprot)
            & (pred["res_index"] >= start)
            & (pred["res_index"] <= end)
        ]

        predicted_hits = int((sub["calibrated_score"] >= THRESH).sum())
        recall = predicted_hits / ep_len if ep_len > 0 else 0.0

        hit_abs = predicted_hits >= K_ABS
        hit_rel = recall >= P_REL
        hit = hit_abs or hit_rel

        results.append({
            "model": MODEL,
            "organism": org,
            "organism_name": organism_name,
            "group": group,
            "type": EPITOPE_TYPE,
            "iedb_id": r["iedb_id"],
            "uniprot": uniprot,
            "start": start,
            "end": end,
            "epitope_len": ep_len,
            "predicted_hits": predicted_hits,
            "recall": round(recall, 3),
            "hit_abs": hit_abs,
            "hit_rel": hit_rel,
            "hit": hit,
            "threshold": THRESH,
            "K_abs": K_ABS,
            "P_rel": P_REL,
            "out_of_range": False,
            "reason": "",
        })

    res_df = pd.DataFrame(results)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    res_df.to_csv(out_path, index=False)

    eval_df = res_df[res_df["out_of_range"] == False].copy()

    print(eval_df.head(10))

    print("\nSUMMARY (EVALUABLE ONLY)")
    print("Model:", MODEL)
    print("ORG:", org)
    print("Group:", group)
    print("Type:", EPITOPE_TYPE)
    print("Threshold:", THRESH, "K_abs:", K_ABS, "P_rel:", P_REL)
    print("Total epitopes:", len(res_df))
    print("Evaluable epitopes:", len(eval_df))
    print("Out-of-range / missing:", int(res_df["out_of_range"].sum()))
    print("  - missing prediction for protein:", missing_pred_protein)
    print("  - end_pos exceeds predicted length:", out_of_range_ct)

    if len(eval_df) > 0:
        print("Hits combined:", int(eval_df["hit"].sum()))
        print("Mean recall:", round(float(eval_df["recall"].mean()), 3))

        by_protein = (
            eval_df
            .groupby(["model", "organism", "group", "type", "uniprot"])
            .agg(
                n_epitopes=("iedb_id", "count"),
                n_hits=("hit", lambda x: int(x.sum())),
                mean_recall=("recall", "mean"),
            )
            .sort_values("n_epitopes", ascending=False)
        )

        by_protein_out = Path(f"tables/{org}/linear/{org}_bepipred3_by_protein.csv")
        by_protein.to_csv(by_protein_out)

        print("\nWrote by-protein table:", by_protein_out)

    print("\nWrote:", out_path)


if __name__ == "__main__":
    main()