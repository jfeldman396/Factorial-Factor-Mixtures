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
Rscript scripts/sample_size/run_moderate_crossloading_sample_size_MAP_intercepts.R
```

The plotting/checkpoint script is:

```sh
Rscript scripts/sample_size/plot_moderate_crossloading_checkpoint.R
```

The shell refresher can be used while a long run is active:

```sh
OUT_DIR=results/moderate_crossloading_joint_mfa_sample_size_p500_25reps_H3H4_G2G3_MAP_intercepts_centered \
  bash scripts/sample_size/refresh_moderate_crossloading_plots.sh
```

Important settings are also recorded in `configs/sample_size_intercepts_centered.env`.

## IFEval Analysis

Run the full IFEval pipeline with:

```sh
bash scripts/ifeval/run_full_analysis.sh
```

The major steps are:

1. Singular-value shelf diagnostic.
2. Held-out predictive likelihood tuning over rank, component-wise `G`, and sparse loading penalty.
3. Selected mixture model refit with MAP refinement.
4. Loading/cross-loading summaries.
5. Factor-score and mixture-profile visualization.
6. Writeup rendering.

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
