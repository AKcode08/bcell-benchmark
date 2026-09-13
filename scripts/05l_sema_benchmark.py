#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd

MODEL = "SEMA"
EPITOPE_TYPE = "linear"

THRESH = 0.51
K_ABS = 5
P_REL = 0.25

VIRUS_ORGS = {"COV", "ENZA"}


def usage():
    print(
        "Usage: python scripts/05f_sema_benchmark.py <ORG1> [<ORG2> ...]",
        file=sys.stderr
    )
    sys.exit(1)


def get_group(org: str) -> str:
    return "virus" if org in VIRUS_ORGS else "bacteria"


def benchmark_one_org(org):

    group = get_group(org)

    pred_path = Path(f"outputs/SEMAi/{org}/linear/sema_indexed.csv")
    gt_path = Path(f"iedb/{org}/linear/linear_ground_truth.csv")

    out_path = Path(f"tables/{org}/linear/{org}_sema_results.csv")

    if not pred_path.exists():
        print(f"Skipping {org}: no prediction file")
        return None

    if not gt_path.exists():
        print(f"Skipping {org}: no ground truth file")
        return None

    pred = pd.read_csv(pred_path)
    gt = pd.read_csv(gt_path)

    required_pred_cols = {"uniprot", "res_index", "raw_score"}
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
        raise SystemExit(f"{org}: pred file missing columns: {missing}")

    if not required_gt_cols.issubset(gt.columns):
        missing = sorted(required_gt_cols - set(gt.columns))
        raise SystemExit(f"{org}: gt file missing columns: {missing}")

    pred["res_index"] = pred["res_index"].astype(int)

    pred_max = pred.groupby("uniprot")["res_index"].max().to_dict()

    rows = []

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

            rows.append({
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

            rows.append({
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
            (pred["uniprot"] == uniprot) &
            (pred["res_index"] >= start) &
            (pred["res_index"] <= end)
        ]

        predicted_hits = int((sub["raw_score"] >= THRESH).sum())

        recall = predicted_hits / ep_len if ep_len > 0 else 0.0

        hit_abs = predicted_hits >= K_ABS
        hit_rel = recall >= P_REL
        hit = hit_abs or hit_rel

        rows.append({
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

    res_df = pd.DataFrame(rows)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    res_df.to_csv(out_path, index=False)

    eval_df = res_df[res_df["out_of_range"] == False].copy()

    print("\nSUMMARY")
    print("Model:", MODEL)
    print("ORG:", org)
    print("Group:", group)
    print("Type:", EPITOPE_TYPE)
    print("Threshold:", THRESH)
    print("Total epitopes:", len(res_df))
    print("Evaluable:", len(eval_df))
    print("Excluded:", int(res_df["out_of_range"].sum()))

    if len(eval_df) > 0:

        print("Hits:", int(eval_df["hit"].sum()))
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

        by_protein_out = Path(
            f"tables/{org}/linear/{org}_sema_by_protein.csv"
        )

        by_protein.to_csv(by_protein_out)

    print("Saved:", out_path)

    return res_df


def main():

    if len(sys.argv) < 2:
        usage()

    orgs = sys.argv[1:]
    all_results = []

    for org in orgs:

        df = benchmark_one_org(org)

        if df is not None:
            all_results.append(df)

    if all_results:

        combined = pd.concat(all_results, ignore_index=True)

        out = Path("tables/sema_linear_all_orgs_results.csv")
        out.parent.mkdir(parents=True, exist_ok=True)

        combined.to_csv(out, index=False)

        print("\nCombined results saved:", out)


if __name__ == "__main__":
    main()
