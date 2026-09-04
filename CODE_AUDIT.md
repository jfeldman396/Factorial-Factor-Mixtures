# Code Audit Notes

This note records the current reproducibility audit status for the repository.

## 2026-09-04 Final Simulation Cleanup

Scope checked:

- Cleaned the sample-size simulation surface around the final signal-support
  grid and the IFEval scripts/writeups.
- Removed old selected sampled-`Z` simulation plots/tables and scratch
  diagnostic scripts that are not part of the current IFEval or final
  simulation replication path.
- Verified the final simulation launcher defaults to 25 Monte Carlo
  replications per setting.
- Verified the final loading design records weak/strong loading ranges,
  cross-loading probability, block-size mode, and loading-support diagnostics
  in each output row.
- Regenerated all paper-facing DGP loading heatmaps and corresponding Lambda
  CSVs for the final grid.
- Parsed all shared R files and one-level scripts under `scripts/*/*.R`.
- Checked the IFEval shell launchers with `zsh -n`.
- Smoke-tested Product MAP, Viroli-Laplace Gibbs, and Viroli-Gaussian Gibbs
  through the final launcher on a tiny toy grid.
- Re-rendered the final simulation design PDF:
  - `writeup/final_simulation_design/final_simulation_design_algorithms.pdf`.

Commands used for the current parse/smoke audit:

```sh
Rscript -e 'files <- c(Sys.glob("R/*.R"), Sys.glob("scripts/*/*.R")); bad <- character(); for (f in files) { ok <- tryCatch({parse(f); TRUE}, error=function(e) {message("PARSE ERROR ", f, ": ", conditionMessage(e)); FALSE}); if (!ok) bad <- c(bad, f) }; cat("Parsed", length(files), "R files; failures", length(bad), "\n"); if (length(bad)) quit(status=1)'
zsh -n scripts/ifeval/run_full_analysis.sh
zsh -n scripts/ifeval/run_ifeval_threshold_analyses.sh
N_VALUES=12 P_VALUES_PRODUCT=20 P_VALUES_GIBBS=20 H_VALUES=2 G_VALUES=2 LOADING_STRENGTHS=weak CROSS_LOADING_PROBS=0.075 BLOCK_SIZE_MODES=balanced REP_VALUES=1 TASK_WORKERS_PRODUCT=1 PRODUCT_INTERNAL_WORKERS=2 TASK_WORKERS_GIBBS=2 VIROLI_ITER=6 VIROLI_BURN=3 VIROLI_COMPUTE_PARAMETER_ESS=FALSE RUN_LABEL=signal_support_grid_smoke_gibbs2 Rscript scripts/sample_size/run_final_product_viroli_simulation.R
```

Result:

```text
Parsed 28 R files; failures 0
Fresh tiny launcher smoke completed with 3 result rows.
```

## 2026-08-15 Static Audit

Scope checked:

- Parsed all shared R files in `R/*.R`.
- Parsed all R scripts one level under `scripts/*/*.R`, including IFEval and
  sample-size simulation entry points.
- Verified the cleaned IFEval inputs exist:
  - `data/ifeval/openeval_ifeval_only_binary_matrix.csv`;
  - `data/ifeval/openeval_item_metadata.csv`;
  - `data/ifeval/openeval_model_metadata.csv`.
- Verified key replication entry points exist:
  - `scripts/ifeval/cv_ifeval_rank_lambda_models.R`;
  - `scripts/ifeval/fit_interpret_ifeval_mixture.R`;
  - `scripts/sample_size/run_final_product_viroli_simulation.R`;
  - `scripts/sample_size/plot_signal_support_simulation_progress.R`;
  - `scripts/sample_size/plot_dgp_loading_heatmaps.R`.
- Smoke-tested the cleaned IFEval CV script with a tiny `H=1, G=1` two-fold run.
- Smoke-tested the cleaned selected-fit script with a tiny `H=1, G=1` run.
- Checked `scripts/ifeval/run_full_analysis.sh` with `zsh -n`.
- Re-rendered the component-wise IFEval writeup PDF:
  - `writeup/ifeval_componentwise_G3313.pdf`.
- Visually inspected the rendered PDF pages containing:
  - component-wise CV plot;
  - marginal mixture plot;
  - LLM factor-score and MAP-profile heatmaps;
  - Lambda heatmap;
  - primary-loading and cross-loading examples.

Commands used for the R parse audit:

```sh
Rscript -e 'files <- c(Sys.glob("R/*.R"), Sys.glob("scripts/*/*.R")); bad <- character(); for (f in files) { ok <- tryCatch({parse(f); TRUE}, error=function(e) {message("PARSE ERROR ", f, ": ", conditionMessage(e)); FALSE}); if (!ok) bad <- c(bad,f) }; cat("Parsed", length(files), "R files; failures", length(bad), "\n"); if (length(bad)) quit(status=1)'
```

Result:

```text
Parsed 17 R files; failures 0
```

## Remaining Audit Caveats

- This was a static code-path and documentation audit, not a fresh end-to-end rerun of the full IFEval cross-validation grid or full sample-size simulation.
- The full column-wise IFEval CV grid is computationally heavier than the selected-model refits. Use `RESUME_EXISTING=TRUE` and checkpointed output directories for full reruns.
- The current sample-size simulation is resumable by task chunks under
  `results/full/signal_support_grid/chunks`; full raw outputs are ignored by
  git and should be regenerated from the launcher.
- The selected-fit script `fit_interpret_ifeval_mixture.R` accepts arbitrary
  `H_FIXED` and scalar or vector `G_FIXED`.
- Full generated outputs under `results/full/` or temporary analysis folders are intentionally not committed. Paper-facing snapshots and writeups are committed.
