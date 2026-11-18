# Psoriasis multi-signature benchmark

Analysis code accompanying the PeerJ Research Article:

> **Limited cross-cohort added value of published psoriasis transcriptomic signatures: a three-cohort reanalysis**

This repository implements the prespecified D2 case-only lesion/nonlesional benchmark described in the manuscript. It covers signature registration, platform gene coverage, equal-weight positive-direction z-score scoring, patient-level incremental discrimination, bootstrap confidence intervals, Spearman redundancy, participation-ratio effective dimension, and figure/supplement generation.

## Scope and interpretation

The code is a transparent comparator benchmark. It does not reconstruct any published classifier exactly. None of the ten registered signatures had both complete weights and intercepts available in the assessed records; the comparator is therefore an equal-weight, positive-direction gene-level z-score. Missing genes are not imputed.

The endpoint is D2: lesional versus nonlesional skin within patients with psoriasis. It is not a lesion-versus-healthy-control diagnostic task, a severity or treatment-response endpoint, or evidence of clinical utility.

## Workflow

| Script | Purpose |
|---|---|
| `scripts/benchmark_gate.mjs` | Register signatures and audit model recoverability and platform coverage. |
| `scripts/benchmark_score.mjs` | Score GSE41664 and GSE6710 and compare with four biological baselines. |
| `scripts/benchmark_increment.mjs` | Patient-level incremental models, leave-one-patient-out predictions, bootstrap, and Delta AUC calls. |
| `scripts/audit_gse154200.mjs` | Audit GSE154200 metadata, treatment, and pairing. |
| `scripts/benchmark_gse154200.mjs` | Run the independent GSE154200 benchmark. |
| `scripts/generate_peerj_figures_supplements.R` | Generate figures and supplemental outputs. |

## Requirements

- Node.js with ES modules support.
- R 4.4.2 or compatible.
- R packages: `readr`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `scales`, `grid`, `jsonlite`, `digest`, and `png`.

## Verify a checkout

Run the repository checks from any working directory:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\verify_repository.ps1
```

The verifier checks Node/R syntax and compares every tracked or publishable file with `SHA256_MANIFEST.csv`. It deliberately ignores local-only files such as `ZENODO_DOI.txt`. After reviewing an intentional code or metadata change, regenerate the manifest once and rerun the verifier:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\verify_repository.ps1 -UpdateManifest
powershell -NoProfile -ExecutionPolicy Bypass -File .\verify_repository.ps1
```

## Run order

```powershell
node scripts/benchmark_gate.mjs
node scripts/benchmark_score.mjs
node scripts/benchmark_increment.mjs
node scripts/audit_gse154200.mjs
node scripts/benchmark_gse154200.mjs
Rscript scripts/generate_peerj_figures_supplements.R
```

The scripts write derived outputs to `data/benchmark/` and `artifacts/`. Raw public inputs must first be placed under `sources/raw/` as documented in `data/input_manifest.csv` and `README_INPUTS.md`.

`run_all.ps1` resolves paths from its own location, so it can be invoked from outside the repository root.

## Prespecified decisions

- Endpoint: D2 lesion/nonlesional pairs only.
- Signature comparator: equal-weight positive-direction z-score unless recoverable model parameters exist.
- Biological baselines: KC differentiation, inflammatory response, IFN-gamma response, and E2F proliferation.
- Incremental model: signature plus the four baseline scores.
- Resampling: patient-level leave-one-patient-out predictions and 2,000 patient-level bootstrap replicates.
- Minimum meaningful Delta AUC: 0.03.
- Missing genes: retained as missing; no imputation.

## Citation and license

Citation metadata are in `CITATION.cff`. The code is released under the MIT license.
