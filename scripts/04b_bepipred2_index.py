#!/usr/bin/env python3
"""
Combine BepiPred-2.0 web-server raw outputs into the shared 13-column benchmark schema.

Input tree (organism/type encoded in folder path):
    {INPUT_ROOT}/{org}/{type}/*raw*ouput*.csv   (or *raw*output*.csv)
    - a dir may hold PT1raw_ouput.csv .. PT4raw_ouput.csv  (web-server batch splits) -> concatenated
    - or a single raw_ouput.csv

Output tree:
    {OUTPUT_ROOT}/{org}/{type}/BepiPred-2.0_{org}_{type}.csv

Shared schema (first 13 cols), then predictor-specific extras appended:
    model, organism, group, type, uniprot, chain, res_index, residue,
    raw_score, calibrated_score, score, threshold, y_pred,
    + exposed_buried, rsa, helix_prob, sheet_prob, coil_prob
"""
import sys, glob, os, re
import pandas as pd

# ============================ CONFIG ============================
MODEL_NAME = "BepiPred-2.0"
THRESHOLD  = 0.5          # confirmed native cutoff
CHAIN_ID   = "A"          # <FLAG> AlphaFold monomers -> single chain 'A'. Confirm vs your other models.

# organism = folder token (ESKAPE/GAST/INTRA/RESP); group = broad category constant.
GROUP_VALUE = "bacteria"         # "bacteria" here; set "virus" for the COV/ENZA sets later

# <FLAG> BepiPred-2.0 web output has no separate calibrated score.
# Set how calibrated_score should be filled to match your other models:
#   "same"  -> copy the epitope probability
#   "blank" -> leave empty (NaN)
CALIBRATED_SCORE_MODE = "same"
# ===============================================================

RAW_GLOBS = ["*raw*ouput*.csv", "*raw*output*.csv"]   # tolerate the 'ouput' misspelling + PTn prefixes


def find_raw_files(type_dir):
    files = []
    for pat in RAW_GLOBS:
        files.extend(glob.glob(os.path.join(type_dir, pat)))
    files = sorted(set(files))  # PT1,PT2,... sort naturally; single file also caught
    return files


def load_and_combine(type_dir):
    files = find_raw_files(type_dir)
    if not files:
        return None, []
    frames = [pd.read_csv(f) for f in files]
    combined = pd.concat(frames, ignore_index=True)
    return combined, [os.path.basename(f) for f in files]


def map_schema(df, org, typ):
    ep = df["EpitopeProbability"].astype(float)
    out = pd.DataFrame()
    out["model"]            = [MODEL_NAME] * len(df)
    out["organism"]         = org
    out["group"]            = GROUP_VALUE
    out["type"]             = typ
    out["uniprot"]          = df["Entry"].values
    out["chain"]            = CHAIN_ID
    out["res_index"]        = df["Position"].astype(int).values
    out["residue"]          = df["AminoAcid"].values
    out["raw_score"]        = ep.values
    if CALIBRATED_SCORE_MODE == "same":
        out["calibrated_score"] = ep.values
    else:
        out["calibrated_score"] = pd.NA
    out["score"]            = ep.values
    out["threshold"]        = THRESHOLD
    out["y_pred"]           = (ep.values >= THRESHOLD).astype(int)
    # predictor-specific extras
    out["exposed_buried"]   = df["Exposed/Buried"].values
    out["rsa"]              = df["RelativeSurfaceAccessilibity"].values
    out["helix_prob"]       = df["HelixProbability"].values
    out["sheet_prob"]       = df["SheetProbability"].values
    out["coil_prob"]        = df["CoilProbability"].values
    return out


def main(input_root, output_root):
    n_written = 0
    for org in sorted(os.listdir(input_root)):
        org_dir = os.path.join(input_root, org)
        if not os.path.isdir(org_dir):
            continue
        for typ in sorted(os.listdir(org_dir)):
            type_dir = os.path.join(org_dir, typ)
            if not os.path.isdir(type_dir):
                continue
            combined, src = load_and_combine(type_dir)
            if combined is None:
                continue

            # integrity: same antigen must not appear in two batch splits
            dup = combined["Entry"].value_counts()
            per_ent_pos_dup = combined.duplicated(subset=["Entry", "Position"]).sum()
            if per_ent_pos_dup > 0:
                print(f"  [WARN] {org}/{typ}: {per_ent_pos_dup} duplicate (Entry,Position) rows across splits {src}")

            out = map_schema(combined, org, typ)
            out = out.sort_values(["uniprot", "res_index"]).reset_index(drop=True)

            dst_dir = os.path.join(output_root, org, typ)
            os.makedirs(dst_dir, exist_ok=True)
            dst = os.path.join(dst_dir, "bepipred2_indexed.csv")
            out.to_csv(dst, index=False)
            n_written += 1
            print(f"  {org}/{typ}: {len(src)} file(s) {src} -> {out['uniprot'].nunique()} antigens, "
                  f"{len(out)} rows -> {dst}")
    print(f"\nDone. {n_written} combined file(s) written under {output_root}/")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python combine_bepipred2.py <INPUT_ROOT> <OUTPUT_ROOT>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
