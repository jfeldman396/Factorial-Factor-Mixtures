# Replication Instructions

This guide explains how to replicate the two analysis tracks in this repository:

1. the IFEval empirical analysis;
2. the sample-size simulation comparing the proposed factorial factor-mixture method with Viroli-style Gibbs samplers.

Run these commands in a local macOS Terminal, not in the GitHub web interface.
All commands assume the repository root is the working directory:

```sh
cd "/Users/joefeldman/Documents/Deep Factor Models/factorial-factor-mixtures"
```

## Software

The code is written primarily in R. The main scripts use:

- base R and `parallel`;
- `MASS`;
- `truncnorm`;
- `glmnet`;
- `coda`;
- `ggplot2`;
- `reshape2`;
- `mclust`.

Optional 3D HTML visualizations use Python packages:

- `numpy`;
- `pandas`;
- `matplotlib`;
- `plotly`.

The scripts were developed on macOS. The sample-size simulation uses Unix-style multicore parallelism where available.

## Repository Layout

Important paths:

- `R/`: shared copies of the proposed binary probit pretraining and refinement algorithms.
- `scripts/ifeval/`: IFEval fitting, tuning, comparison, and visualization scripts.
- `scripts/sample_size/`: simulation, Gibbs comparator, and checkpoint plotting scripts.
- `data/ifeval/`: cleaned IFEval binary matrix and item/model metadata.
- `configs/sample_size_intercepts_centered.env`: the main simulation grid/settings.
- `results/selected_tables/`: selected CSV outputs committed for review.
- `results/selected_plots/`: selected figures committed for review.
- `writeup/`: IFEval writeups and rendered PDFs.
- `CODE_AUDIT.md`: current static audit notes, parse checks, and caveats.

Full regenerated outputs should be written under `results/full/`, which is ignored by git.

## IFEval Analysis

The cleaned input files are:

- `data/ifeval/openeval_ifeval_only_binary_matrix.csv`;
- `data/ifeval/openeval_item_metadata.csv`;
- `data/ifeval/openeval_model_metadata.csv`.

The current IFEval workflow fits only the proposed independent-mixture probit
factor model.  Rank `H`, component counts `G_h`, and the entrywise sparse
loading penalty are tuned by held-out predictive log likelihood on randomly
removed model-by-item cells.  The current component-wise grid allows
`G_h in {1,2,3}` with at most one Gaussian coordinate (`G_h = 1`).

### Run The Current IFEval Pipeline

The easiest way to reproduce the current analysis is:

```sh
bash scripts/ifeval/run_full_analysis.sh
```

This runs held-out-likelihood CV, selected mixture refit, Lambda heatmaps,
marginal mixture plots, loading summaries, and factor-score visualizations.

### Tune H, Component-Wise G, And Sparse Loading Penalty

```sh
OUT_DIR=results/full/ifeval_columnwise_G_cv_atmost1_gaussian \
H_GRID=1:5 \
G_MODE=column_grid \
G_GRID=1,2,3 \
G_COMPONENT_VALUES=1,2,3 \
MAX_GAUSSIAN_COORDS=1 \
LAMBDA_L1_GRID=0,1,2,4,8,12 \
K_FOLDS=3 \
WORKERS=6 \
PRETRAIN_AUG_ITER=20 \
REFINE_ITER=20 \
MIXTURE_MAX_ITER=200 \
FIT_SELECTED_AFTER_CV=TRUE \
SAVE_FITS=FALSE \
RESUME_EXISTING=TRUE \
REFRESH_PLOTS=TRUE \
Rscript scripts/ifeval/cv_ifeval_rank_lambda_models.R
```

Primary outputs:

- `ifeval_rank_lambda_cv_fold_scores.csv`;
- `ifeval_rank_lambda_cv_histories.csv`;
- `ifeval_rank_lambda_cv_summary.csv`;
- `ifeval_rank_lambda_selected_by_heldout_ll.csv`;
- `rank_lambda_top_candidates_by_heldout_ll.csv`;
- `rank_lambda_top_candidates_by_heldout_ll.png`;
- selected full-data mixture refit folder when `FIT_SELECTED_AFTER_CV=TRUE`.

### Refit The Interpretable Component-Wise Model

The current component-wise interpretation uses:

- `H = 4`;
- `G = (3,3,1,3)`, so the third coordinate is Gaussian;
- sparse-loading MAP refinement with `lambda_l1_penalty = 4`;
- `20` maximum pretraining iterations and `20` maximum refinement iterations.

Run:

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

Primary outputs:

- `openeval_item_intercepts_loadings_metadata.csv`;
- `openeval_model_factor_scores_profiles.csv`;
- `openeval_factor_interpretation_summary.csv`;
- `openeval_top_loading_item_examples.csv`;
- `openeval_factor_mixture_groups.csv`;
- loading, marginal-mixture, factor-score, and profile plots.

### Summarize Loadings And Cross-Loadings

```sh
LOADINGS_PATH=results/full/ifeval/reproduced_componentwise_H4_G3313_lambda4/openeval_item_intercepts_loadings_metadata.csv \
OUT_DIR=results/full/ifeval/loadings_crossloadings \
Rscript scripts/ifeval/summarize_ifeval_loadings_crossloadings.R
```

Primary outputs:

- `ifeval_sparse_learned_lambda_all_items_ordered.csv`;
- `ifeval_sparse_lambda_matrix_threshold_0p5.csv`;
- `ifeval_G3_factor_only_exact_items.csv`;
- `ifeval_G3_cross_loading_items.csv`;
- heatmaps for sparse and cross-loading items.

### Regenerate Selected-Fit Visualizations

```sh
FIT_DIR=results/full/ifeval/reproduced_componentwise_H4_G3313_lambda4 \
OUT_DIR=results/full/ifeval/selected_visualizations \
Rscript scripts/ifeval/plot_openeval_mixture_lambda_heatmaps.R

FIT_DIR=results/full/ifeval/reproduced_componentwise_H4_G3313_lambda4 \
OUT_DIR=results/full/ifeval/selected_visualizations \
Rscript scripts/ifeval/plot_factor_marginal_mixtures.R

FIT_DIR=results/full/ifeval/reproduced_componentwise_H4_G3313_lambda4 \
OUT_DIR=results/full/ifeval/ifeval_3d_factor_visualizations \
Rscript scripts/ifeval/plot_ifeval_learned_factors_baseR.R
```

Optional Python/Plotly plots:

```sh
FIT_DIR=results/full/ifeval/reproduced_componentwise_H4_G3313_lambda4 \
OUT_DIR=results/full/ifeval/ifeval_3d_factor_visualizations \
python3 scripts/ifeval/plot_ifeval_learned_factors.py
```

The rendered component-wise writeup is:

```text
writeup/ifeval_componentwise_G3313.pdf
```

It summarizes the selected fit, compact cross-loading examples, and LLM
ability profiles by factor.

## Sample-Size Simulation

The current final simulation uses deterministic EM-SVD likelihood pretraining
for the proposed product-mixture method and compares it against two
Viroli-style Gibbs samplers.  The design and algorithms are summarized in:

```text
writeup/final_simulation_design/final_simulation_design_algorithms.pdf
```

The final launcher is:

```sh
Rscript scripts/sample_size/run_final_product_viroli_simulation.R
```

By default this runs:

- `n in {100, 200, 300, 400}`;
- `p in {500, 1000, 2000, 4000}` for product MAP;
- `p in {500, 1000}` for the Viroli baselines, to control runtime;
- `H in {5, 10, 15, 20}`;
- `G in {2, 3}`, expanded to equal component counts across factors;
- separation values `{1, 2}`;
- 25 Monte Carlo repetitions;
- IFEval-like unbalanced Cross loadings;
- item intercepts using the IFEval-like intercept design;
- loading-based sign/permutation alignment for all recovery metrics;
- Product MAP with EM-SVD likelihood pretraining, rotation, and MAP
  refinement;
- Viroli Gibbs with a Laplace loading prior and penalty 10;
- Viroli Gibbs with a diffuse Gaussian loading prior.

The main output folder is:

```text
results/full/final_cross_ifeval_product_viroli
```

The run is resumable.  Each scenario-method pair is keyed in
`comparison_results_checkpoint.csv`; rerunning the launcher with
`RESUME_EXISTING=TRUE` skips completed fits.  The checkpoint records the full
DGP setting, penalties, mixture priors, convergence diagnostics, elapsed
seconds, ESS summaries for Gibbs draws, and DGP diagnostics such as item
prevalence and observed joint profiles.

### Final Simulation Smoke Test

For a fast code-path check:

```sh
OUT_DIR=results/diagnostics/final_launcher_smoke \
N_VALUES=20 \
P_VALUES_PRODUCT=30 \
P_VALUES_VIROLI=30 \
BLOCK_SIZE_MODES=balanced \
H_VALUES=2 \
G_VALUES=2 \
SEPARATIONS=1 \
REP_VALUES=1 \
EM_SVD_ITER=3 \
ROTATION_ITER=2 \
REFINE_ITER=3 \
MIXTURE_MAX_ITER=10 \
VIROLI_ITER=20 \
VIROLI_BURN=10 \
VIROLI_COMPUTE_PARAMETER_ESS=FALSE \
PARALLEL_OURS=FALSE \
PARALLEL_GIBBS=FALSE \
WRITE_PARAMETER_TABLES=TRUE \
Rscript scripts/sample_size/run_final_product_viroli_simulation.R
```

This smoke test only verifies the product MAP, Viroli-Laplace, and
Viroli-Gaussian code paths.  It is not a scientific simulation.

### Legacy Sampled-Z Simulation

The earlier sampled-Z simulation compares:

- `independent_marginal_mixture`: the proposed method with sampled augmented-probit SVD pretraining, independent marginal mixture rotation, item intercepts, and MAP refinement;
- `joint_mixture_factor_gibbs`: a correctly specified joint-mixture factor model with `G^H` latent profiles and diagonal within-profile covariance, fitted by full Gibbs sampling with probit augmentation.

The main grid is:

- `H in {3, 4}`;
- `G in {2, 3}`;
- `n in {100, 500, 1000, 2000}`;
- `p in {250, 500, 1000, 2000}`;
- `25` Monte Carlo repetitions;
- loading designs `balanced_moderate_few_positive_cross` and `balanced_moderate_dense_signed_cross`;
- block-size modes `balanced` and `moderate_ifeval_like`;
- item intercept mode `ifeval_like`;
- unequal mixture variances;
- proposed-method pretraining max iterations `10`;
- proposed-method refinement max iterations `10`;
- Gibbs iterations `2000`, burn-in `1000`, thin `1`, for `p in {250, 500}` only.

The canonical launcher runs the balanced block grid first, then the moderately
IFEval-like unbalanced block grid:

```sh
PARALLEL_OURS=TRUE \
PARALLEL_GIBBS=FALSE \
PARALLEL_WORKERS=18 \
Rscript scripts/sample_size/run_sampledZ_full_pgrid_smallp_gibbs_MAP_intercepts.R
```

Outputs are written to:

```text
results/full/sampledZ_pgrid_balanced_crossloading_smallp_gibbs_MAP_intercepts
results/full/sampledZ_pgrid_unbalanced_crossloading_smallp_gibbs_MAP_intercepts
```

The launcher uses `RESUME_EXISTING=TRUE`, so rerunning the command resumes from
existing checkpoint files rather than starting over.

Primary outputs:

- `comparison_results_checkpoint.csv`;
- `comparison_results.csv`;
- `comparison_summary.csv`;
- `comparison_convergence.csv`;
- per-repetition parameter-recovery CSV files;
- Gibbs history CSV/PNG files for `p in {250, 500}`;
- checkpoint/final line plots.

### Run A Smoke Test

For a fast check of the code path:

```sh
cd "/Users/joefeldman/Documents/Deep Factor Models/factorial-factor-mixtures"

OUT_DIR=results/full/smoke_sample_size \
H_VALUES=3 \
G_VALUES=2 \
NP_GRID=n100p250:100:250 \
LOADING_DESIGNS=balanced_moderate_few_positive_cross \
REP_VALUES=1 \
PRETRAIN_AUG_ITER=3 \
REFINE_ITER=3 \
MFA_ITER=50 \
MFA_BURN=25 \
RESUME_EXISTING=FALSE \
Rscript scripts/sample_size/compare_original_simulation_joint_mfa_gibbs.R
```

This is only a code-path check, not a scientific simulation.

### Parallelization

The current long-run recommendation is to parallelize the proposed method with
18 workers and keep Gibbs serial:

```sh
PARALLEL_OURS=TRUE \
PARALLEL_GIBBS=FALSE \
PARALLEL_WORKERS=18 \
Rscript scripts/sample_size/run_sampledZ_full_pgrid_smallp_gibbs_MAP_intercepts.R
```

For the proposed product-mixture method, this parallelizes independent
within-iteration tasks: marginal mixture fits, rotation subproblems where
available, itemwise loading regressions, and subject-wise MAP factor-score
updates. The outer pretraining/refinement iterations and SVD are still serial.

For the joint-mixture Gibbs comparator, the implementation can parallelize
subject factor draws and item regression draws, but focused checks found this
slower in the current R implementation because serial MCMC steps and
worker/data-copy overhead dominate. The full grid therefore leaves Gibbs
serial and skips Gibbs entirely for `p > 500`.

### Regenerate Paper-Facing Sample-Size Figures

The selected paper-facing simulation plots are generated from the checkpoint
CSV and the exact DGP loading generator. Set `OUT_DIR` to the full simulation
output directory you want to refresh:

```sh
OUT_DIR=results/full/sampledZ_pgrid_balanced_crossloading_smallp_gibbs_MAP_intercepts
```

For the moderately IFEval-like unbalanced block simulation, use:

```sh
OUT_DIR=results/full/sampledZ_pgrid_unbalanced_crossloading_smallp_gibbs_MAP_intercepts
```

To regenerate representative DGP loading heatmaps for `Loadings = "Sparse"`
and `Loadings = "Cross"` at `H = 3, 4`:

```sh
OUT_DIR=$OUT_DIR \
Rscript scripts/sample_size/plot_dgp_loading_heatmaps.R
```

To regenerate the six-panel recovery figures for each loading/H/G setting:

```sh
OUT_DIR=$OUT_DIR \
Rscript scripts/sample_size/plot_sample_size_rmse_panels.R
```

These figures report RMSE for item intercepts, loadings, mixture means,
mixture variances, mixture weights, and factor scores. The mixture-weight RMSE
is computed by vectorizing the aligned joint-profile weights and taking the
RMSE against the true profile weights.

For method-comparison RMSE panels, colors distinguish methods and line types
distinguish the two dimensions where both methods are fit: dotted lines are
`p = 250` and solid lines are `p = 500`. Larger `p` values are product-mixture
only in the full simulation and are therefore excluded from those direct
Gibbs-vs-product panels.

To regenerate the product-mixture factor-score convergence plots across `n`,
with separate colored lines for each `p`:

```sh
Rscript scripts/sample_size/plot_factor_score_rmse_by_n_product_by_p.R
```

These figures use the committed aggregate result CSVs in
`results/selected_tables/sample_size`, color lines by `p`, and plot
mean factor-score RMSE with `+/- 2` standard deviations across repetitions.

To regenerate runtime figures:

```sh
OUT_DIR=$OUT_DIR \
Rscript scripts/sample_size/plot_sample_size_timing_lines.R
```

These figures plot mean `log(seconds)` with `+/- 2` standard deviations across
Monte Carlo repetitions. The `seconds` field is total elapsed wall-clock time
for each method fit; for the Gibbs comparator it includes the configured 2000
Gibbs iterations and 1000 burn-in iterations.

### Interpret Simulation Metrics

The main recovery metrics are computed after aligning estimated factors and
loadings to the data-generating factors and loadings by best signed
permutation.

Important fields:

- `factor_score_rmse`: RMSE between aligned estimated and true factor scores.
- `lambda_rmse`: RMSE between aligned estimated and true loading entries.
- `alpha_rmse`: RMSE between estimated and true item intercepts.
- `mixture_mu_rmse`: RMSE between estimated and true marginal mixture means.
- `mixture_var_rmse`: RMSE between estimated and true marginal mixture variances.
- `mixture_weight_rmse`: RMSE between vectorized aligned joint-profile weights.
- `seconds`: wall-clock runtime for the method in that repetition.

For `G = 3`, variance recovery can be sensitive to saturated binary item
blocks, so inspect variance RMSE together with the generated DGP Lambda
heatmaps and factor-score RMSE.

## Committed Selected Results

Selected outputs from prior runs are committed for immediate inspection:

- sample-size plots: `results/selected_plots/sample_size`;
- sample-size tables: `results/selected_tables/sample_size`;
- aggregate raw sample-size results:
  `results/selected_tables/sample_size/balanced_sampledZ_pgrid_product_allp_gibbs_smallp_comparison_results.csv`
  and
  `results/selected_tables/sample_size/unbalanced_sampledZ_pgrid_product_allp_gibbs_smallp_comparison_results.csv`;
- IFEval plots: `results/selected_plots/ifeval`;
- IFEval tables: `results/selected_tables/ifeval`;
- IFEval writeups: `writeup/ifeval_analysis_writeup.pdf` and
  `writeup/ifeval_componentwise_G3313.pdf`.

The latest static audit notes are in `CODE_AUDIT.md`. They record which R files
were parsed, which input files were checked, which PDF was rendered, and which
caveats remain before claiming a full fresh end-to-end rerun.

These selected outputs are snapshots. Full reruns should write to `results/full/`.

## Git Notes

Generated `.rds` files, logs, and full output directories are ignored by `.gitignore`. Commit source code, replication instructions, selected summary tables, selected plots, and final writeups.
