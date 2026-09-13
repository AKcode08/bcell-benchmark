#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import requests
import argparse
import time
import re

REQUEST_TIMEOUT = 30
SLEEP_BETWEEN_REQUESTS = 0.2

AF_FILE_BASE = "https://alphafold.ebi.ac.uk/files"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_PDB_FILE_URL = "https://files.rcsb.org/download"


def normalize_uniprot_id(x):
    if pd.isna(x):
        return None
    s = str(x).strip()
    if not s:
        return None
    s = re.split(r"[;,]", s)[0].strip()
    if "|" in s:
        parts = s.split("|")
        if len(parts) >= 2 and parts[1].strip():
            s = parts[1].strip()
    s = s.upper()
    return s if s else None


def load_unique_uniprots(csv_path: Path, column: str):
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()

    if column not in df.columns:
        raise ValueError(
            f"Column '{column}' not found in {csv_path}. "
            f"Available columns: {df.columns.tolist()}"
        )

    ids = (
        df[column]
        .dropna()
        .astype(str)
        .map(normalize_uniprot_id)
        .dropna()
        .unique()
        .tolist()
    )
    return sorted(ids)


def candidate_af_urls(uniprot_id: str):
    af_id = f"AF-{uniprot_id}-F1"
    return [
        (f"{AF_FILE_BASE}/{af_id}-model_v6.pdb", "v6"),
        (f"{AF_FILE_BASE}/{af_id}-model_v4.pdb", "v4"),
    ]


def try_download_alphafold(uniprot_id: str, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)

    for url, version in candidate_af_urls(uniprot_id):
        try:
            r = requests.get(url, stream=True, timeout=REQUEST_TIMEOUT)
        except requests.RequestException as e:
            return {
                "download_status": "af_request_failed",
                "source": "AlphaFold",
                "version": None,
                "url": url,
                "error": str(e),
                "structure_path": None,
                "structure_id": None,
            }

        if r.status_code == 200:
            outfile = outdir / f"{uniprot_id}.pdb"
            with open(outfile, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            return {
                "download_status": "downloaded",
                "source": "AlphaFold",
                "version": version,
                "url": url,
                "error": None,
                "structure_path": str(outfile),
                "structure_id": f"AF-{uniprot_id}-F1",
            }

        if r.status_code not in (403, 404):
            return {
                "download_status": f"af_http_{r.status_code}",
                "source": "AlphaFold",
                "version": version,
                "url": url,
                "error": r.text[:200],
                "structure_path": None,
                "structure_id": None,
            }

    return {
        "download_status": "not_found",
        "source": "AlphaFold",
        "version": None,
        "url": None,
        "error": None,
        "structure_path": None,
        "structure_id": None,
    }


def search_rcsb_by_uniprot(uniprot_id: str, rows: int = 10):
    payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id,
            },
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": False,
            "results_content_type": ["experimental"],
            "paginate": {"start": 0, "rows": rows},
            "sort": [{"sort_by": "score", "direction": "desc"}],
        },
    }

    r = requests.post(RCSB_SEARCH_URL, json=payload, timeout=REQUEST_TIMEOUT)
    r.raise_for_status()
    data = r.json()
    results = data.get("result_set", [])
    return [x["identifier"] for x in results if "identifier" in x]


def try_download_pdb(uniprot_id: str, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        pdb_ids = search_rcsb_by_uniprot(uniprot_id)
    except requests.RequestException as e:
        return {
            "download_status": "pdb_search_failed",
            "source": "PDB",
            "version": None,
            "url": None,
            "error": str(e),
            "structure_path": None,
            "structure_id": None,
        }

    if not pdb_ids:
        return {
            "download_status": "not_found",
            "source": "PDB",
            "version": None,
            "url": None,
            "error": None,
            "structure_path": None,
            "structure_id": None,
        }

    for pdb_id in pdb_ids:
        url = f"{RCSB_PDB_FILE_URL}/{pdb_id}.pdb"
        try:
            r = requests.get(url, stream=True, timeout=REQUEST_TIMEOUT)
        except requests.RequestException as e:
            return {
                "download_status": "pdb_request_failed",
                "source": "PDB",
                "version": None,
                "url": url,
                "error": str(e),
                "structure_path": None,
                "structure_id": pdb_id,
            }

        if r.status_code == 200:
            outfile = outdir / f"{uniprot_id}_{pdb_id}.pdb"
            with open(outfile, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            return {
                "download_status": "downloaded",
                "source": "PDB",
                "version": None,
                "url": url,
                "error": None,
                "structure_path": str(outfile),
                "structure_id": pdb_id,
            }

    return {
        "download_status": "not_found",
        "source": "PDB",
        "version": None,
        "url": None,
        "error": None,
        "structure_path": None,
        "structure_id": None,
    }


def try_download_structure(uniprot_id: str, outdir: Path):
    af_res = try_download_alphafold(uniprot_id, outdir)
    if af_res["download_status"] == "downloaded":
        return af_res

    pdb_res = try_download_pdb(uniprot_id, outdir)
    if pdb_res["download_status"] == "downloaded":
        return pdb_res

    return {
        "download_status": "not_found",
        "source": "none",
        "version": None,
        "url": None,
        "error": af_res.get("error") or pdb_res.get("error"),
        "structure_path": None,
        "structure_id": None,
    }


def download_set(uniprot_ids, outdir: Path, organism: str, dataset_type: str):
    rows = []

    for i, uid in enumerate(uniprot_ids, start=1):
        print(f"[{organism} | {dataset_type}] {i}/{len(uniprot_ids)}: {uid}")

        existing = list(outdir.glob(f"{uid}*.pdb"))
        if existing:
            rows.append({
                "organism": organism,
                "dataset_type": dataset_type,
                "uniprot_id": uid,
                "download_status": "already_exists",
                "source": "existing",
                "alphafold_version": None,
                "structure_id": None,
                "download_url": None,
                "structure_path": str(existing[0]),
                "error": None,
            })
            time.sleep(SLEEP_BETWEEN_REQUESTS)
            continue

        res = try_download_structure(uid, outdir)

        rows.append({
            "organism": organism,
            "dataset_type": dataset_type,
            "uniprot_id": uid,
            "download_status": res["download_status"],
            "source": res["source"],
            "alphafold_version": res["version"],
            "structure_id": res["structure_id"],
            "download_url": res["url"],
            "structure_path": res["structure_path"],
            "error": res["error"],
        })

        time.sleep(SLEEP_BETWEEN_REQUESTS)

    manifest = pd.DataFrame(rows)
    manifest.to_csv(outdir / "structure_download_manifest.csv", index=False)
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism", required=True)
    parser.add_argument("--iedb-root", default="iedb")
    parser.add_argument("--structures-root", default="structures")
    parser.add_argument("--uniprot-column", default="uniprot_accession")
    args = parser.parse_args()

    organism = args.organism

    linear_csv = Path(args.iedb_root) / organism / "linear" / "linear_ground_truth.csv"
    conformational_csv = Path(args.iedb_root) / organism / "conformational" / "conformational_ground_truth.csv"

    linear_out = Path(args.structures_root) / organism / "linear"
    conformational_out = Path(args.structures_root) / organism / "conformational"

    linear_ids = load_unique_uniprots(linear_csv, args.uniprot_column)
    conformational_ids = load_unique_uniprots(conformational_csv, args.uniprot_column)

    print(f"{organism} linear unique UniProt IDs: {len(linear_ids)}")
    print(f"{organism} conformational unique UniProt IDs: {len(conformational_ids)}")

    manifest_lin = download_set(linear_ids, linear_out, organism, "linear")
    manifest_con = download_set(conformational_ids, conformational_out, organism, "conformational")

    print("\nLinear summary:")
    print(manifest_lin["download_status"].value_counts(dropna=False))

    print("\nConformational summary:")
    print(manifest_con["download_status"].value_counts(dropna=False))

    print("\nLinear source summary:")
    print(manifest_lin["source"].value_counts(dropna=False))

    print("\nConformational source summary:")
    print(manifest_con["source"].value_counts(dropna=False))


if __name__ == "__main__":
    main()