# Benchmarking B-cell Epitope Prediction Across Bacterial Antigens

Residue- and region-level evaluation of seven B-cell epitope predictors on AlphaFold-modelled
bacterial antigens, using experimentally determined epitope annotations from the IEDB.

This repository contains the dataset construction, predictor execution and analysis code for:

> Kumar, A. *Benchmarking B-cell Epitope Prediction Across Bacterial Antigens: Improved Performance
> Metrics with Residue- and Region-Level Evaluation Across Seven Approaches and Algorithms.*
> MSc Bioinformatics thesis, University of Birmingham, September 2026.

---

## Aim

Reported B-cell epitope predictor performance is difficult to compare across studies because
methods are developed on different antigen sets, epitope definitions, homology filters and
residue subsets. Separately, the metrics in common use (ROC-AUC, PR-AUC) score residues
individually and do not establish whether a predictor's highest-scoring residues fall inside the
same experimentally annotated epitope.

This benchmark tests one hypothesis: **residue-level discrimination and epitope-region recovery
measure distinct aspects of predictor performance.** All seven predictors are evaluated on a
single antigen set under a common framework, addressing three questions:

1. Which predictors rank epitope residues most effectively?
2. Does residue-level discrimination translate into recovery of individual linear or
   conformational epitopes?
3. Does region recovery exceed what is expected from random residue selection?

Two design choices carry most of the methodological weight:

- **Matched prediction count.** For every antigen, each predictor is allowed exactly `N` positive
  predictions, where `N` is that antigen's annotated epitope-residue count. This removes
  threshold and calibration differences between methods and isolates ranking quality.
- **Matched random baseline.** Predictor scores are permuted within each antigen, the top `N`
  are retaken, and recall is recomputed (20 repeats). Recovery is reported as *observed minus
  random*, because positional and spatial tolerance inflate both.

Linear and conformational epitopes are evaluated separately throughout, on **common-coverage
sets** restricted to antigens and residue positions scored by all seven predictors.

---

## Predictors benchmarked

| Predictor | Year | Architecture (encoder → network) | Structure input | Executed via |
|---|---|---|---|---|
| ElliPro | 2008 | Ellipsoid fit → protrusion index → distance clustering (no ML) | Yes | Next-Gen IEDB Tools API |
| BepiPred-2.0 | 2017 | Sequence-derived RSA features → random forest | No | DTU Health Tech server |
| BepiPred-3.0 | 2022 | ESM-2 (650M) → feed-forward NN (+RSA, +length) | No | DTU Health Tech server |
| GraphBepi | 2023 | ESM-2 + AlphaFold graph + DSSP → BiLSTM ∥ edge-enhanced GNN → MLP | Yes | Self-hosted (BlueBEAR) |
| DiscoTope-3.0 | 2024 | ESM-IF1 inverse folding + RSA/pLDDT → XGBoost, PU learning | Yes | DTU Health Tech server |
| SEMA-3D 2.0 | 2024 | Ensemble of 5 SaProt (3Di structure-aware vocabulary) → linear head | Yes | Self-hosted (BlueBEAR) |
| RoBep | 2026 | ESM-C + RSA/dihedrals → E(n)-equivariant GNN + region constraint | Yes | Forge API + BlueBEAR |

Together these span geometric, classical machine-learning, protein-language-model and
graph-based approaches. Per-predictor parameters, checkpoints and access dates are in the
manuscript Methods and Supplementary Tables S8–S9.

---

## Dataset and workflow

```
IEDB export  ──▶  UniProt REST  ──▶  AlphaFold DB v6  ──▶  7 predictors  ──▶  common schema  ──▶  analysis
(annotations)     (canonical seq)    (AF-{acc}-F1 v6)      (sequence or       (13 columns)       (notebooks)
                                                            structure in)
```

1. **Epitope annotations.** IEDB B-cell assay records with positive experimental outcomes, no
   host or disease restriction. Linear (continuous) and conformational (discontinuous) epitopes
   are kept as separate annotation classes.
2. **Sequences.** Canonical sequences from the UniProt REST API (release 2026_01), keyed by the
   accession in the IEDB cross-reference and deduplicated so multiple annotations map to one
   sequence. Residue positions use 1-based UniProt numbering throughout.
3. **Structures.** Single-fragment `AF-{accession}-F1` model, **version 6**, from the AlphaFold
   Protein Structure Database, matched to the canonical sequence used for annotation. Per-residue
   pLDDT is read from the Cα B-factor field.
4. **Labelling.** Residues inside an annotated epitope are positive; all other resolved residues
   are negative for benchmarking purposes. (See *Limitations* — negatives are unlabelled, not
   verified non-epitopes.)
5. **Standardisation.** Every predictor's output is mapped to 1-based UniProt positions and
   converted to a continuous per-residue score where higher means greater epitope propensity.

### Antigen attrition

| Stage | Antigens |
|---|---|
| Retrievable UniProt sequence | 364 |
| Excluded after ground-truth filtering | −7 |
| Structure retrieval attempted | 357 |
| No structure available | −6 |
| Experimental PDB structure only (excluded) | −1 |
| **AlphaFold DB v6 models retained** | **350** |
| Linear-epitope common-coverage set | **319** (119,464 residues) |
| Conformational-epitope common-coverage set | **36** (13,513 residues) |

The structural dataset contains 2,920 linear and 148 conformational epitope annotations
(mean size 14.5 and 8.4 residues).

### Organism groups

Analyses use the four bacterial groups. `COV` and `ENZA` are curated in the tree but excluded
from the benchmark (multi-chain complexity); the notebooks set `ORGS = ["ESKAPE", "GAST",
"INTRA", "RESP"]`.

| Code | Group | Species |
|---|---|---|
| `ESKAPE` | Priority nosocomial | *E. faecium, S. aureus, K. pneumoniae, A. baumannii, P. aeruginosa, Enterobacter* spp. |
| `GAST` | Gastrointestinal | *S. enterica, S. flexneri, C. jejuni, V. cholerae* |
| `INTRA` | Intracellular & zoonotic | *L. monocytogenes*, *Brucella* spp., *C. burnetii*, *C. trachomatis* |
| `RESP` | Respiratory | *M. tuberculosis, S. pneumoniae, H. influenzae, B. pertussis* |
| `COV`, `ENZA` | Viral (curated, not analysed) | — |

---

## Analysis pipeline

| Level | Unit | Measure | Matching rule | Baseline |
|---|---|---|---|---|
| Residue ranking | Residue | ROC-AUC, PR-AUC (micro + macro) | Exact residue | — |
| Linear recovery | Continuous epitope | Per-epitope recall, detection rate | ±0–3 sequence positions, plus 8 Å Cβ criterion for k>0 | Random selection at same `N` |
| Conformational recovery | Discontinuous epitope | Per-epitope recall, detection rate | 0–8 Å minimum Cβ distance | Random selection at same `N` |
| Structural confidence | Residue | ROC-AUC, PR-AUC | pLDDT strata (<50, 50–<70, 70–<90, ≥90) | Band-specific prevalence |
| Predictor consensus | Residue | Precision, recall | Agreement ≥1 to ≥7 predictors | Conformational prevalence (0.063) |

Distances are measured between Cβ atoms, using Cα for glycine. Detection rate `D(τ)` is the
fraction of epitopes with recall ≥ τ, evaluated at τ = 0.25, 0.50 and 0.75.

### Shared output schema

Every predictor is reduced to the same 13 leading columns, with predictor-specific fields
appended after:

```
model, organism, group, type, uniprot, chain, res_index, residue,
raw_score, calibrated_score, score, threshold, y_pred
```

---

## Repository structure

```
bcell_benchmark/
├── iedb/{ORG}/{linear|conformational}/
│     ├── {type}_raw.csv                 IEDB export as downloaded
│     ├── {type}_meta_raw.csv            assay / reference metadata
│     ├── {type}_receptor_raw.csv        receptor records
│     └── {type}_ground_truth.csv        parsed labels used by every analysis
├── fasta/{ORG}/{type}/{ACCESSION}.fasta per-antigen UniProt sequences (443 files)
├── inputs/{ORG}/{type}/                 multi-FASTA batches for web-server submission
│                                        (*_sequences.fasta, plus part1/part2 splits)
├── structures/{ORG}/{type}/{ACC}.pdb    AlphaFold DB v6 models (425 files)
├── outputs/{predictor}/{ORG}/{type}/
│     ├── <raw per-antigen predictor output>
│     └── {predictor}_indexed.csv        standardised 13-column table
├── tables/{ORG}/{type}/                 per-group benchmark tables
│     └── *_results.csv, *_by_protein.csv
├── tables/*_all_orgs*.csv               pooled benchmark tables
├── figures/                             manuscript figures (.png + .pdf)
├── scripts/                             pipeline scripts 01–05l and notebooks 06–08
├── slurm/                               BlueBEAR job scripts
├── logs/                                job logs
└── requirements.txt
```

`{predictor}` is one of `ellipro`, `bepipred2`, `bepipred3`, `graphbepi`, `discotope3`,
`sema`, `robep`.

### Scripts

| Script | Purpose |
|---|---|
| `01_fasta_download.sh <ORG_TYPE>` | Download UniProt FASTA for the accessions in an IEDB export |
| `02_make_multifasta.sh <ORG>` | Combine per-antigen FASTA into one multi-FASTA, headers = accession only |
| `02a_split_fasta.sh <ORG>` | Split a multi-FASTA into parts for web-server size limits |
| `03_build_ground_truth.py <ORG>` | Parse linear epitope annotations → `linear_ground_truth.csv` |
| `03a_build_ground_truth_disc.py <ORG>` | Parse conformational annotations → `conformational_ground_truth.csv` |
| `03b_pdb_download.py --organism <ORG>` | Fetch AlphaFold DB v6 models (falls back to RCSB lookup) |
| `03c_zip_structures.sh` | Archive structure directories for upload to web servers |
| `04_bepipred_index_raw.py` | Index BepiPred raw output (linear) into the shared schema |
| `04a_bepipred_index_disc.py` | Index BepiPred raw output (conformational) |
| `04b_bepipred2_index.py <IN_ROOT> <OUT_ROOT>` | Combine BepiPred-2.0 web-server batch splits |
| `05*_combine_*_csv.py` | Merge per-antigen predictor outputs into `{predictor}_indexed.csv` |
| `05_bepipred_benchmark.py <ORG>`, `05a_..._disc.py <ORG>` | BepiPred benchmark tables — **one organism per call** |
| `05c/05d`, `05f/05g`, `05i/05j`, `05l` `<ORG...>` | DiscoTope, ElliPro, GraphBepi, SEMA benchmark tables — accept an organism list |
| `05h_graphbepi_index.py` | Index GraphBepi output |

### Notebooks

All notebooks run from `scripts/` and read `../outputs`, `../iedb`, `../structures`.

| Notebook | Produces |
|---|---|
| `06_residue_level.ipynb` | Residue-level ROC-AUC / PR-AUC (Table 4, Fig. 2), pLDDT stratification (Fig. 6), `conformational_consensus.csv` (Table S7) |
| `06a_linear_level.ipynb` | Linear region recovery vs random, detection rates (Fig. 3, Table S3) |
| `06b_conformational_level.ipynb` | Conformational region recovery, predictor agreement and spatial overlap (Fig. 4, Fig. 7, Tables S4, S7) |
| `07_dataset_characterization.ipynb` | Dataset composition, length and pLDDT distributions (Tables 1–2, Fig. S1) |
| `08_bias_analysis_patched.ipynb` | Antigen-property Spearman correlations and OLS models (Fig. S2, Table S5) |

---

## Environment setup

Predictor execution and analysis use different environments. Predictors that were run through
public web servers require no local installation.

### Analysis environment (notebooks 06–08)

The environment used for the reported results (manuscript Table S8):

| | Version |
|---|---|
| Python | 3.10.19 |
| NumPy | 2.2.6 |
| pandas | 2.3.3 |
| scikit-learn | 1.7.2 |
| SciPy | 1.15.3 |
| matplotlib | 3.10.8 |

```bash
conda create -n metabolab python=3.10.19
conda activate metabolab
pip install numpy==2.2.6 pandas==2.3.3 scikit-learn==1.7.2 scipy==1.15.3 \
            matplotlib==3.10.8 seaborn biopython jupyterlab
```

> **Note:** the `requirements.txt` at the repository root describes a separate Python 3.13
> `venv/` used for plotting only. It omits `scikit-learn` and `scipy`, which the notebooks
> import, so it is **not** sufficient to reproduce the analysis. Use the specification above.

### Self-hosted predictors (BlueBEAR)

| Predictor | Python | Model weights | Key libraries |
|---|---|---|---|
| GraphBepi | 3.9.23 | `BCE_633_GraphBepi/model_-1.ckpt`; ESM-2 `esm2_t36_3B_UR50D` | PyTorch 1.12.1, fair-esm 2.0.0, PyTorch-Lightning 1.6.4, NumPy 1.21.5, pandas 1.4.2 |
| SEMA-3D | 3.10.0 | `sema_3d_0.pth`–`sema_3d_4.pth`; `westlake-repl/SaProt_650M_PDB` | PyTorch 2.4.1, transformers 4.37.2, PyG 2.6.1, Biopython 1.78 |
| RoBep | 3.10.20 | `best_mcc_model.bin`; ESM-C configuration | PyTorch 2.5.0, transformers 4.46.3, esm 3.1.3, PyG 2.6.1, scikit-learn 1.5.2 |

Each predictor runs in its own isolated conda environment. Additional dependencies: **DSSP
4.4.0** (secondary structure for GraphBepi) and **Foldseek** (3Di tokens for SEMA-3D; the
commit bundled with the SEMAi repository was used). Jobs were managed with **SLURM 26.05.1** on
BlueBEAR (RHEL 8.10), CPU-only inference, no internet access on compute nodes.

> **BlueBEAR note:** a bare `conda activate` fails inside a SLURM script — the shell is
> non-interactive and conda's hook is not loaded. Either source the hook first
> (`source $(conda info --base)/etc/profile.d/conda.sh` before `conda activate`, as
> `slurm/run_graphbepi_batch.slurm` does), or call the environment's Python by absolute path.
> Compute nodes are CPU-only and have no internet access, so models and weights must be cached
> in advance.

---

## Reproducing the analysis

Run from the repository root unless stated otherwise. Steps 1–3 rebuild the dataset; steps 4–5
regenerate predictions; step 6 reproduces every figure and table in the manuscript. **If you
only want to reproduce the reported results, start at step 6** — the standardised
`{predictor}_indexed.csv` files and ground-truth tables are already in the tree.

### 1. Sequences

```bash
./scripts/01_fasta_download.sh ESKAPE_linear
./scripts/02_make_multifasta.sh ESKAPE
./scripts/02a_split_fasta.sh ESKAPE        # only if a server rejects the full file
```

### 2. Ground truth

```bash
for ORG in ESKAPE GAST INTRA RESP; do
  ./scripts/03_build_ground_truth.py "$ORG"
  ./scripts/03a_build_ground_truth_disc.py "$ORG"
done
```

### 3. Structures

```bash
for ORG in ESKAPE GAST INTRA RESP; do
  ./scripts/03b_pdb_download.py --organism "$ORG"
done
```

### 4. Predictions

Four predictors were run through public services and cannot be reproduced offline; results
reflect the server state on the access date (manuscript Table S9).

| Predictor | Service | Accessed |
|---|---|---|
| BepiPred-2.0 | DTU Health Tech BepiPred-2.0 server | June 2026 |
| BepiPred-3.0 | DTU Health Tech BepiPred-3.0 server | March 2026 |
| DiscoTope-3.0 | DTU Health Tech DiscoTope-3.0 server | March 2026 |
| ElliPro | Next-Generation IEDB Tools API | April 2026 |
| RoBep (ESM-C embeddings only) | EvolutionaryScale Forge API, `esmc-6b-2024-12` | May 2026 |

Submit `inputs/{ORG}/{type}/*_sequences.fasta` (sequence-based methods) or the archives from
`03c_zip_structures.sh` (structure-based methods), and place the returned files under
`outputs/{predictor}/{ORG}/{type}/`.

The remaining three — GraphBepi, SEMA-3D and RoBep — were self-hosted rather than run through a
web server. "Locally" in the manuscript means *on BlueBEAR*, not on a workstation: each predictor's
upstream implementation and checkpoints live under `models/` in the cluster project tree
(`/rds/projects/e/elhamsak-epitope-dev/bcell_benchmark/models/`), which is **not** part of this
clone. Clone each upstream repository there and set `PROJECT_DIR` before submitting:

```bash
# reads structures/{ORG}/{type}, writes outputs/graphbepi/{ORG}/{type}
PATHOGENS=ESKAPE,GAST,INTRA,RESP TYPES=linear,conformational \
  sbatch slurm/run_graphbepi_batch.slurm

sbatch slurm/run_graphbepi_eskape_conf.slurm   # single-group re-run
```

RoBep runs in two phases: Phase 1 requests ESM-C embeddings from the Forge API (internet is
unavailable on compute nodes, so this runs from a login node and caches to
`data/embeddings/esmc/`), Phase 2 batch-predicts offline.

### 5. Standardise outputs

```bash
python scripts/04b_bepipred2_index.py outputs/bepipred2 outputs/bepipred2
python scripts/05b_combine_discotope_csv.py
python scripts/05e_combine_ellipro_csv.py
python scripts/05k_combine_sema_csv.py
python scripts/05h_graphbepi_index.py --organism ESKAPE
```

Each writes `outputs/{predictor}/{ORG}/{type}/{predictor}_indexed.csv`. The `05*_benchmark*.py`
scripts then emit per-group tables into `tables/`:

```bash
# 05 and 05a take exactly one organism; the rest accept a list
for ORG in ESKAPE GAST INTRA RESP; do
  python scripts/05_bepipred_benchmark.py "$ORG"
  python scripts/05a_bepipred_benchmark_disc.py "$ORG"
done

python scripts/05c_discotope_benchmark.py ESKAPE GAST INTRA RESP
python scripts/05d_discotope_benchmark_disc.py ESKAPE GAST INTRA RESP
python scripts/05f_ellipro_benchmark.py ESKAPE GAST INTRA RESP
python scripts/05g_ellipro_benchmark_disc.py ESKAPE GAST INTRA RESP
python scripts/05i_graphbepi_benchmark.py ESKAPE GAST INTRA RESP
python scripts/05j_graphbepi_benchmark_disc.py ESKAPE GAST INTRA RESP
python scripts/05l_sema_benchmark.py ESKAPE GAST INTRA RESP
```

### 6. Analysis and figures

```bash
conda activate metabolab
cd scripts
jupyter lab
```

Execute in order: `06_residue_level` → `06a_linear_level` → `06b_conformational_level` →
`07_dataset_characterization` → `08_bias_analysis_patched`.

> **Note:** `06`, `06a` and `08` write figures to the working directory (`scripts/`); `06b` and
> `07` write to `../figures/`. Move the former into `figures/` after running, or launch Jupyter
> from the repository root and adjust the relative paths.

---

## Key results

| | Linear | Conformational |
|---|---|---|
| Best micro ROC-AUC | BepiPred-3.0, 0.64 | GraphBepi and BepiPred-3.0, 0.65 |
| Best region recovery over random | BepiPred-3.0, **+0.044** | RoBep **+0.094**, GraphBepi +0.072, DiscoTope-3.0 +0.053 |
| Legacy methods | ElliPro 0.51, BepiPred-2.0 0.53 | At or near chance |

- No predictor exceeded ROC-AUC 0.65 on either epitope type.
- BepiPred-3.0 matched GraphBepi on conformational residue-level ROC-AUC (0.65) yet recovered
  conformational epitopes only +0.005 above random — the central dissociation the benchmark
  was designed to detect.
- Positional and spatial tolerance raised observed recall, but raised the random baseline
  faster: at k ≥ 1 (linear) and d ≥ 6 Å (conformational) every predictor fell below its matched
  random baseline. Tolerance-based recovery is uninterpretable without a matched control.
- Agreement among ≥5 of 7 predictors raised conformational precision to 0.31 (≈5× the 0.063
  prevalence) while recall fell to 0.05. This uses the known epitope-residue count to set `N`
  and is therefore not a deployable ensemble.

---

## Limitations

- **Unlabelled negatives.** Residues without an IEDB annotation are not verified non-epitopes;
  an antigen may have been tested against few antibodies. Precision and other measures that
  penalise false positives are therefore likely underestimated.
- **Antigen-only prediction.** An epitope depends on the binding antibody. These predictors
  estimate general epitope propensity without modelling a paratope.
- **Scale.** The conformational set contains 36 antigens. Differences between structure-based
  predictors need confirmation on larger independent data.
- **Scope.** Bacterial antigens only; unbound monomeric AlphaFold models, so oligomeric
  interfaces, glycosylation and antibody-induced conformational change are not represented.

---

## Repository notes

Items a reader may otherwise trip over:

- **No predictor model code in this clone.** `scripts/` holds dataset construction, output
  standardisation and benchmarking only. The upstream implementations and checkpoints for
  GraphBepi, SEMA-3D and RoBep, plus RoBep's own two-phase driver scripts, live under `models/`
  in the BlueBEAR project tree (`/rds/projects/e/elhamsak-epitope-dev/bcell_benchmark/`); only
  results were synced here. Sync the RoBep drivers in — or document their upstream commits —
  before deposition, since RoBep is the only predictor with no local script covering its run.
- **RoBep output filenames are lowercased accessions** (`a0a0h3ggm3.json`, not `A0A0H3GGM3.json`)
  because the upstream implementation lowercases internally. The `robep_indexed.csv` step
  restores canonical UniProt casing; match case-insensitively if you read the raw files.
- **`outputs/robep/` has no `COV`/`ENZA`.** Expected — RoBep was only run on the four bacterial
  groups used in the analysis.
- **`05k_combine_sema_csv.py` points at `outputs/SEMAi`** while the directory is `outputs/sema`;
  adjust `BASE_DIR` before running. Its `THRESHOLD` constant (0.51) also differs from the SEMA-3D
  epitope threshold used elsewhere in the project (0.361, inclusive `>=`).
- **`02a_split_fasta.sh` expects a filename that is not on disk.** It reads
  `inputs/{ORG}/linear/{ORG}_sequences.fasta`, but the files written by step 1 are named
  `{ORG}_linear_sequences.fasta`. Rename or adjust `IN=` before running it.
- **`05l_sema_benchmark.py`** prints a usage string referring to `05f_sema_benchmark.py`.
- **`08_bias_analysis_patched.ipynb`** is the version corresponding to the reported results;
  `08_bias_analysis.ipynb` is retained as the pre-correction copy.
- **`venv/`** and `.DS_Store` files should be excluded before publication.

---

## Citation

If you use this code or the benchmark dataset, please cite:

```bibtex
@mastersthesis{kumar2026bcell,
  author = {Kumar, Aman},
  title  = {Benchmarking B-cell Epitope Prediction Across Bacterial Antigens:
            Improved Performance Metrics with Residue- and Region-Level Evaluation
            Across Seven Approaches and Algorithms},
  school = {University of Birmingham},
  year   = {2026},
  month  = {September},
  type   = {{MSc} thesis}
}
```

<!-- TODO: replace with the journal citation once the manuscript is published. -->

Please also cite the underlying resources — IEDB, UniProt, the AlphaFold Protein Structure
Database — and the individual predictors, listed in the manuscript References.

## Contact

**Aman Kumar** — MSc Bioinformatics, University of Birmingham
amanaastha.ak@gmail.com

<!-- TODO: add institutional email, ORCID and the public repository URL before deposition. -->

Supervisors: Dr. Andreas Bender, Dr. Mohamed El-Hadidi.

## Acknowledgements

Thanks to Dr. Vlad Cojocaru, Patricia-Raluca Trăistaru and Ana Lucanu (Babeș-Bolyai University)
and Yuan Hao (Khalifa University) for technical discussion and feedback. This work used the
University of Birmingham's Birmingham Environment for Academic Research (BEAR), including the
BlueBEAR high-performance computing facility.

## License

Code in this repository is released under the MIT License — see [`LICENSE`](LICENSE).

Derived data are subject to the terms of their sources: epitope annotations from the
[IEDB](https://www.iedb.org/), sequences from [UniProt](https://www.uniprot.org/) (CC BY 4.0),
and structures from the [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/)
(CC BY 4.0). Each benchmarked predictor carries its own licence; consult the original
repositories before redistribution.
