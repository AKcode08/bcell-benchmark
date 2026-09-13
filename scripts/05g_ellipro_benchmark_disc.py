#!/usr/bin/env python3

import sys
import ast
from pathlib import Path
import pandas as pd

MODEL = "ElliPro"
EPITOPE_TYPE = "conformational"

THRESH = 0.50
K_ABS = 2
P_REL = 0.25

VIRUS_ORGS = {"COV", "ENZA"}


def usage():
    print(
        "Usage: python scripts/06d_ellipro_conformational_benchmark.py <ORG1> [<ORG2> ...]",
        file=sys.stderr
    )
    sys.exit(1)


def get_group(org: str) -> str:
    return "virus" if org in VIRUS_ORGS else "bacteria"


def parse_residue_set(x):
    if pd.isna(x):
        return set()

    try:
        parsed = ast.literal_eval(str(x))

        if isinstance(parsed, (list, tuple, set)):
            return {int(v) for v in parsed}

    except Exception:
        pass

    return set()


def benchmark_one_org(org):

    group = get_group(org)

    pred_path = Path(
        f"outputs/ellipro/{org}/conformational/ellipro_indexed.csv"
    )

    gt_path = Path(
        f"iedb/{org}/conformational/conformational_ground_truth.csv"
    )

    out_path = Path(
        f"tables/{org}/conformational/{org}_ellipro_results.csv"
    )

    if not pred_path.exists():
        print(f"Skipping {org}: missing predictions")
        return None

    if not gt_path.exists():
        print(f"Skipping {org}: missing ground truth")
        return None

    pred = pd.read_csv(pred_path)
    gt = pd.read_csv(gt_path)

    required_pred_cols = {"uniprot", "res_index", "raw_score"}

    required_gt_cols = {
        "iedb_id",
        "uniprot_accession",
        "organism",
        "residue_string",
        "residue_set",
        "n_residues",
    }

    if not required_pred_cols.issubset(pred.columns):
        missing = sorted(required_pred_cols - set(pred.columns))
        raise SystemExit(f"{org}: pred file missing columns: {missing}")

    if not required_gt_cols.issubset(gt.columns):
        missing = sorted(required_gt_cols - set(gt.columns))
        raise SystemExit(f"{org}: gt file missing columns: {missing}")

    pred["res_index"] = pred["res_index"].astype(int)

    pred_max = pred.groupby("uniprot")["res_index"].max().to_dict()

    pred_pos = (
        pred.loc[pred["raw_score"] >= THRESH, ["uniprot", "res_index"]]
        .groupby("uniprot")["res_index"]
        .apply(lambda s: set(map(int, s.tolist())))
        .to_dict()
    )

    rows = []

    missing_pred_protein = 0
    all_residues_out_of_range_ct = 0
    partial_out_of_range_ct = 0
    empty_residue_set_ct = 0

    for _, r in gt.iterrows():

        organism_name = str(r["organism"]).strip()
        uniprot = str(r["uniprot_accession"]).strip()
        iedb_id = str(r["iedb_id"]).strip()

        residue_string = str(r["residue_string"]).strip()
        true_set = parse_residue_set(r["residue_set"])
        n_residues = int(r["n_residues"])

        if not true_set:

            empty_residue_set_ct += 1

            rows.append({
                "model": MODEL,
                "organism": org,
                "organism_name": organism_name,
                "group": group,
                "type": EPITOPE_TYPE,
                "iedb_id": iedb_id,
                "uniprot": uniprot,
                "residue_string": residue_string,
                "n_residues_total": n_residues,
                "n_residues_in_range": None,
                "fraction_in_range": None,
                "predicted_hits": None,
                "recall": None,
                "precision": None,
                "hit_abs": None,
                "hit_rel": None,
                "hit": None,
                "threshold": THRESH,
                "K_abs": K_ABS,
                "P_rel": P_REL,
                "partially_out_of_range": None,
                "out_of_range": True,
                "reason": "empty_true_residue_set",
            })

            continue

        max_i = pred_max.get(uniprot)

        if max_i is None:

            missing_pred_protein += 1

            rows.append({
                "model": MODEL,
                "organism": org,
                "organism_name": organism_name,
                "group": group,
                "type": EPITOPE_TYPE,
                "iedb_id": iedb_id,
                "uniprot": uniprot,
                "residue_string": residue_string,
                "n_residues_total": n_residues,
                "n_residues_in_range": None,
                "fraction_in_range": None,
                "predicted_hits": None,
                "recall": None,
                "precision": None,
                "hit_abs": None,
                "hit_rel": None,
                "hit": None,
                "threshold": THRESH,
                "K_abs": K_ABS,
                "P_rel": P_REL,
                "partially_out_of_range": None,
                "out_of_range": True,
                "reason": "missing_prediction_for_protein",
            })

            continue

        valid_true = {x for x in true_set if x <= int(max_i)}
        n_valid = len(valid_true)

        if n_valid == 0:

            all_residues_out_of_range_ct += 1

            rows.append({
                "model": MODEL,
                "organism": org,
                "organism_name": organism_name,
                "group": group,
                "type": EPITOPE_TYPE,
                "iedb_id": iedb_id,
                "uniprot": uniprot,
                "residue_string": residue_string,
                "n_residues_total": n_residues,
                "n_residues_in_range": 0,
                "fraction_in_range": 0.0,
                "predicted_hits": None,
                "recall": None,
                "precision": None,
                "hit_abs": None,
                "hit_rel": None,
                "hit": None,
                "threshold": THRESH,
                "K_abs": K_ABS,
                "P_rel": P_REL,
                "partially_out_of_range": False,
                "out_of_range": True,
                "reason": f"all_true_residues_out_of_range(max_true={max(true_set)},pred_len={int(max_i)})",
            })

            continue

        partially_out_of_range = n_valid < len(true_set)

        if partially_out_of_range:
            partial_out_of_range_ct += 1

        pred_set = pred_pos.get(uniprot, set())

        overlap = valid_true & pred_set

        predicted_hits = len(overlap)

        recall = predicted_hits / n_valid if n_valid > 0 else 0.0
        precision = predicted_hits / len(pred_set) if len(pred_set) > 0 else 0.0

        hit_abs = predicted_hits >= K_ABS
        hit_rel = recall >= P_REL
        hit = hit_abs or hit_rel

        rows.append({
            "model": MODEL,
            "organism": org,
            "organism_name": organism_name,
            "group": group,
            "type": EPITOPE_TYPE,
            "iedb_id": iedb_id,
            "uniprot": uniprot,
            "residue_string": residue_string,
            "n_residues_total": n_residues,
            "n_residues_in_range": n_valid,
            "fraction_in_range": round(n_valid / len(true_set), 3),
            "predicted_hits": predicted_hits,
            "recall": round(recall, 3),
            "precision": round(precision, 3),
            "hit_abs": hit_abs,
            "hit_rel": hit_rel,
            "hit": hit,
            "threshold": THRESH,
            "K_abs": K_ABS,
            "P_rel": P_REL,
            "partially_out_of_range": partially_out_of_range,
            "out_of_range": False,
            "reason": "partial_in_range_evaluation" if partially_out_of_range else "",
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
        print("Mean precision:", round(float(eval_df["precision"].mean()), 3))

        by_protein = (
            eval_df
            .groupby(["model", "organism", "group", "type", "uniprot"])
            .agg(
                n_epitopes=("iedb_id", "count"),
                n_hits=("hit", lambda x: int(x.sum())),
                mean_recall=("recall", "mean"),
                mean_precision=("precision", "mean"),
            )
            .sort_values("n_epitopes", ascending=False)
        )

        by_protein_out = Path(
            f"tables/{org}/conformational/{org}_ellipro_by_protein.csv"
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

        res = benchmark_one_org(org)

        if res is not None:
            all_results.append(res)

    if all_results:

        combined = pd.concat(all_results, ignore_index=True)

        out = Path("tables/ellipro_conformational_all_orgs.csv")
        out.parent.mkdir(parents=True, exist_ok=True)

        combined.to_csv(out, index=False)

        print("\nCombined results saved:")
        print(out)
        print("Total rows:", len(combined))


if __name__ == "__main__":
    main()