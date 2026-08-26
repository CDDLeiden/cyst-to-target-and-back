# From cyst to target and back

Code, data, trained models, and analysis outputs accompanying the manuscript
*“From Phenotypic Screening to Target and Compound Prioritization for Autosomal
Dominant Polycystic Kidney Disease.”*

The repository is organized in the order of the study workflow. Precomputed
notebook outputs and figures are committed so that the reported results can be
inspected without rerunning computationally expensive or network-dependent
steps.

## Quick start

### 1. Clone the complete repository

[Git LFS](https://git-lfs.com/) is required for four large CSV files.
Install Git LFS before cloning, then run:

```shell
git lfs install
git clone https://github.com/CDDLeiden/cyst-to-target-and-back.git
cd cyst-to-target-and-back
git lfs pull
git lfs fsck
```

If you already cloned without Git LFS, install it and run `git lfs pull` from
inside the clone.

### 2. Create the locked environment

Install [Pixi](https://pixi.sh/latest/installation/), then run from the
repository root:

```shell
pixi install --locked
```

Pixi installs the exact Linux or macOS packages recorded in `pixi.lock`; there
is no separate Conda activation or pip installation step. The lock currently
supports `linux-64` and Apple Silicon (`osx-arm64`).

Use `pixi shell` if you want an interactive shell inside the environment.

## GPU model training

A GPU is not needed to inspect results or rerun the phenotypic analyses. The locked default environment intentionally uses CPU
XGBoost on both supported platforms. Full QSAR retraining in
`modelTrain-adora1_nr3c2.ipynb` was performed with XGBoost 2.1.0 on NVIDIA GPUs.
To retrain, use a Linux environment with the matching GPU-enabled XGBoost/CUDA
build and change the notebook's `device` values to an available device (for
example, `cuda:0`). The exact training/validation subsets, optimized model
files, parameters, metrics, and Optuna study databases are committed, so model
claims can be inspected without GPU retraining.

## What can be reproduced?

The repository supports three different levels of verification:

1. **Inspect the reported outputs.** The notebooks retain their executed cells,
   plots, model metrics, and statistical results. Figures, trained models,
   Optuna studies, and processed tables are also committed.
2. **Rerun analyses from supplied data.** The screening, target-validation,
   compound-exploration, gene-expression, and model-analysis inputs are
   included. Use the working directories listed below.
3. **Regenerate upstream datasets.** Some steps use live PubChem, ChEMBL, or
   UniProt services and Papyrus 05.6. These can be slower and may differ as
   external services evolve. SmallWorld analogue searches and supplier stock
   searches are represented by their supplied result tables; the commercial
   queries themselves are not replayed by this repository. Experimental assay
   generation is outside the computational workflow.

See [`data/README.md`](data/README.md) for file-level descriptions and
provenance.

## Workflow and execution order

Notebook paths are relative to their own numbered directory. Start Jupyter from
that directory (not from the repository root) so those paths resolve correctly.

| Stage | Purpose | How to inspect or rerun |
|---|---|---|
| `01_ADPKD_ScreeningProcess` | Reprocess the published phenotypic screens and identify hits | From the repository root: `pixi run python 01_ADPKD_ScreeningProcess/screening_process.py`, then `pixi run python 01_ADPKD_ScreeningProcess/hit_identification.py`. These steps use live PubChem/ChEMBL queries; the exact precomputed hit table is supplied in `data/adpkd_screening/identified_hits/`. |
| `02_TargetID_and_Prioritization` | Link screened compounds to Papyrus 05.6 bioactivities and prioritize targets | `pixi run jupyter-target-identification`. For the supplied-data route, use the notebook section **Loading from saved bioactivityLinker**. The first-time route downloads Papyrus and uses UniProt. |
| `02_TargetID_and_Prioritization/gene_expression` | Plot expression of prioritized targets in mIMCD3 cells | `pixi run jupyter-gene-expression`. |
| `03_Target_Validation` | Normalize, plot, and statistically analyze the target-validation screen | Run `pixi run jupyter-target-validation`. The normalization audit can be regenerated with `pixi run audit-normalization`; it writes the [per-well raw and normalized values](data/normalization_audit/per_well_raw_and_normalized.csv) and [within-plate Mann–Whitney audit](data/normalization_audit/within_plate_mann_whitney_audit.csv). See [`03_Target_Validation/README.md`](03_Target_Validation/README.md). |
| `04_Virtual_Screening` | Train/evaluate A1AR and MR QSAR models, screen analogues, and analyze selected compounds | `pixi run jupyter-virtual-screening`, then open `modelTrain-adora1_nr3c2.ipynb`, `virtualScreening.ipynb`, or `ADPKD-ExplorationScreening-analysis.ipynb`. Full training requires a configured GPU; committed models and outputs allow inspection without retraining. See [`04_Virtual_Screening/README.md`](04_Virtual_Screening/README.md). |
| `05_A1R_Screening` | Report radioligand displacement results for selected A1AR compounds | See [`05_A1R_Screening/README.md`](05_A1R_Screening/README.md). This is an experimental result summary, not a computational pipeline. |

To register the environment as a selectable Jupyter kernel if needed:

```shell
pixi run python -m ipykernel install --user --name cystToTarget --display-name "cystToTarget"
```

## Repository layout

- [`data/`](data/) — supplied raw/processed computational inputs and audit tables
- [`figures/`](figures/) — manuscript figures and molecular drawings
- [`paper_tables/`](paper_tables/) — machine-readable manuscript tables
- [`src/papyrusBioactivityLinker/`](src/papyrusBioactivityLinker/) — local
  compound-to-Papyrus linking package
- `01_...` through `05_...` — ordered analysis stages described above

## Key supplied outputs

- Target-prioritization figures:
  [`figures/ADPKD-Prioritized_Target_Set-prioritized-geneName.svg`](figures/ADPKD-Prioritized_Target_Set-prioritized-geneName.svg)
  and
  [`figures/ADPKD-Prioritized_Target_Set-prioritized-proteinName.svg`](figures/ADPKD-Prioritized_Target_Set-prioritized-proteinName.svg)
- Target-validation statistical results and figures:
  [`03_Target_Validation/README.md`](03_Target_Validation/README.md)
- QSAR, virtual-screening, and compound-exploration results:
  [`04_Virtual_Screening/README.md`](04_Virtual_Screening/README.md)
- A1AR predicted and experimental affinity summary:
  [`05_A1R_Screening/README.md`](05_A1R_Screening/README.md)

## License

The code is released under the [MIT License](LICENSE). Please cite the
accompanying paper when using this repository or its data.
