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

Run the full IFEval pipeline with:

```sh
bash scripts/ifeval/run_full_analysis.sh
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
