# Reproduction Guide

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

1. Rank selection and predictive comparison across candidate ranks.
2. Ordinary binary probit factor model fit for the selected rank.
3. Mixture binary probit factor fit with tuned sparse loading penalty.
4. Loading/cross-loading summaries.
5. 3D factor-score visualization.
6. LaTeX writeup rendering.

Selected generated tables and plots are stored under `results/selected_tables/ifeval` and `results/selected_plots/ifeval`.

## Notes For Audit

The scripts are intentionally close to the working analysis versions. The next code-review pass should split data generation, fitting, alignment, metrics, and plotting into smaller functions with unit tests. See `FUNCTION_MAP.md` for the current function inventory.
