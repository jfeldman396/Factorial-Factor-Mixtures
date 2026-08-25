# Reproduction Guide

For the most complete instructions, see `REPLICATION.md`.

All commands below assume the working directory is the repository root:

```sh
cd "/Users/joefeldman/Documents/Deep Factor Models/factorial-factor-mixtures"
```

## R Package Requirements

The current scripts use base R plus common CRAN packages including `MASS`, `truncnorm`, `ggplot2`, `reshape2`, and `mclust`. Some plotting scripts also use Python packages such as `pandas`, `numpy`, `matplotlib`, and `plotly`.

## Sample-Size Simulation

The active simulation launcher is:

```sh
Rscript scripts/sample_size/run_sampledZ_full_pgrid_smallp_gibbs_MAP_intercepts.R
```

The canonical plotting scripts are:

```sh
OUT_DIR=results/full/sampledZ_pgrid_balanced_crossloading_smallp_gibbs_MAP_intercepts \
Rscript scripts/sample_size/plot_dgp_loading_heatmaps.R

OUT_DIR=results/full/sampledZ_pgrid_balanced_crossloading_smallp_gibbs_MAP_intercepts \
Rscript scripts/sample_size/plot_sample_size_rmse_panels.R

OUT_DIR=results/full/sampledZ_pgrid_balanced_crossloading_smallp_gibbs_MAP_intercepts \
Rscript scripts/sample_size/plot_sample_size_timing_lines.R
```

Use the unbalanced output directory in the same commands to refresh the
moderately IFEval-like block-size results:

```sh
OUT_DIR=results/full/sampledZ_pgrid_unbalanced_crossloading_smallp_gibbs_MAP_intercepts
```

The main launcher uses sampled augmented `Z` pretraining, MAP refinement,
`p in {250,500,1000,2000}`, and runs the Gibbs comparator only for
`p in {250,500}`. Important settings are also recorded in
`configs/sample_size_intercepts_centered.env`.

## IFEval Analysis

### Recreate The Analysis Matrix

The fitted IFEval analysis uses `data/ifeval/openeval_ifeval_only_binary_matrix.csv`.
This is a binary model-by-item matrix.  The entries are formed from OpenEval
responses by extracting the numeric `ifeval_strict_accuracy` values from each
nested score object, averaging within model-item pairs, and coding the pair as
correct when the average score is at least `0.5`.

To rebuild the source OpenEval matrix from Hugging Face, install the data
dependencies and run:

```sh
python3 -m pip install pandas pyarrow huggingface_hub
```

```sh
python3 scripts/data/format_openeval_binary_matrix.py \
  --benchmarks ifeval \
  --out-dir data/openeval_ifeval_formatted_uncapped \
  --min-item-response-prop 0.25 \
  --min-model-response-prop 0.25
```

If the OpenEval snapshot is already cached locally, the same command can be
run without network access by adding:

```sh
--local-snapshot-dir "$HOME/.cache/huggingface/hub/datasets--human-centered-eval--OpenEval/snapshots/<snapshot-id>"
```

The formatter writes both response data and item provenance:

```text
data/openeval_ifeval_formatted_uncapped/openeval_binary_matrix_raw.csv
data/openeval_ifeval_formatted_uncapped/openeval_response_long_scores.csv
data/openeval_ifeval_formatted_uncapped/openeval_item_metadata.csv
data/openeval_ifeval_formatted_uncapped/openeval_item_instruction_metadata_long.csv
```

`openeval_item_metadata.csv` includes the raw OpenEval item payload plus
parsed columns for `prompt`, `instruction_ids`, `instruction_families`,
`n_instructions`, and `instruction_kwargs`.  The long instruction metadata has
one row per item-instruction pair, which is the table to use when relating the
25 retained IFEval instruction ids to the lower-dimensional factor solution.

Then create the model-facing IFEval matrix by starting from the full
model-by-item matrix, selecting IFEval columns, removing low-coverage models
and items, requiring complete item columns for the retained models, and
retaining only nonconstant items:

```sh
OPENEVAL_FULL_MATRIX=data/openeval_ifeval_formatted_uncapped/openeval_binary_matrix_raw.csv \
OPENEVAL_ITEM_METADATA=data/openeval_ifeval_formatted_uncapped/openeval_item_metadata.csv \
OPENEVAL_ITEM_INSTRUCTION_METADATA=data/openeval_ifeval_formatted_uncapped/openeval_item_instruction_metadata_long.csv \
OUT_DIR=data/ifeval \
Rscript scripts/ifeval/create_ifeval_analysis_matrix.R
```

The currently committed uncapped build starts from `124` models and `541`
IFEval items.  It drops `2` low-coverage models, drops `0` item columns for
coverage or missingness after that model filter, and drops `7` all-correct
constant item columns.  The final analysis matrix has `122` models and `534`
retained IFEval items.  The exact counts are stored in
`data/ifeval/ifeval_analysis_matrix_build_summary.csv`.  The retained
`score >= 0.5` analysis metadata contains `25` unique instruction ids, `9`
instruction families, and `820` retained item-instruction rows.

To build threshold-specific matrices for IFEval strict-accuracy scores
`score >= 0.5`, `score >= 2/3`, and `score >= 1.0`, run:

```sh
scripts/ifeval/build_ifeval_threshold_matrices.sh
```

This creates:

```text
data/ifeval_threshold_0p5
data/ifeval_threshold_0p67
data/ifeval_threshold_1
```

The `0p67` label denotes the exact two-thirds cutoff.  A literal threshold of
`0.67` would exclude scores equal to `2/3` and would therefore produce the same
matrix as `score >= 1.0` for this snapshot.

The regenerated threshold-specific analysis matrices have the following
dimensions:

| Folder | Strict-accuracy rule | Models | Retained prompts | Mean binary score |
| --- | --- | ---: | ---: | ---: |
| `data/ifeval_threshold_0p5/` | `score >= 0.5` | 122 | 534 | 0.757 |
| `data/ifeval_threshold_0p67/` | `score >= 2/3` | 122 | 538 | 0.653 |
| `data/ifeval_threshold_1/` | `score >= 1` | 122 | 534 | 0.619 |

To run the full threshold sensitivity analysis:

```sh
scripts/ifeval/run_ifeval_threshold_analyses.sh
```

The resulting fits are written under
`results/full/ifeval_threshold_sensitivity`.

Run the full IFEval pipeline with:

```sh
zsh scripts/ifeval/run_full_analysis.sh
```

The major steps are:

1. Held-out predictive likelihood tuning over rank, component-wise `G`, and sparse loading penalty.
2. Selected mixture model refit with MAP refinement.
3. Loading/cross-loading summaries.
4. Factor-score and mixture-profile visualization.
5. Writeup rendering.

Selected generated tables and plots are stored under `results/selected_tables/ifeval` and `results/selected_plots/ifeval`.

The current component-wise interpretation reported in the writeup uses
`H = 4`, `G = (3,3,1,3)`, and sparse-loading MAP refinement with
`lambda_l1_penalty = 4`. To reproduce that selected fit directly:

```sh
MATRIX_PATH=data/ifeval/openeval_ifeval_only_binary_matrix.csv \
ITEM_METADATA_PATH=data/ifeval/openeval_item_metadata.csv \
OUT_DIR=results/full/ifeval/reproduced_componentwise_H4_G3313_lambda4 \
H_FIXED=4 \
G_FIXED=3,3,1,3 \
WORKERS=8 \
PRETRAIN_AUG_ITER=20 \
REFINE_ITER=20 \
MIXTURE_MAX_ITER=200 \
REQUIRE_MIXTURE_CONVERGENCE=TRUE \
REFINEMENT_LAMBDA_L1_PENALTY=4 \
Rscript scripts/ifeval/fit_interpret_ifeval_mixture.R
```

The compact writeup for this fit is `writeup/ifeval_componentwise_G3313.pdf`.

## Notes For Audit

The scripts are intentionally close to the working analysis versions. The next
code-review pass should split data generation, fitting, alignment, metrics, and
plotting into smaller functions with unit tests. See `FUNCTION_MAP.md` for the
current function inventory and `CODE_AUDIT.md` for the latest static audit.
