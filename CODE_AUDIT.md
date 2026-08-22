# Code Audit Notes

This note records the current reproducibility audit status for the repository.

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
  - `scripts/sample_size/run_sampledZ_full_pgrid_smallp_gibbs_MAP_intercepts.R`;
  - `scripts/sample_size/plot_sample_size_rmse_panels.R`.
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
- The selected-fit script `fit_interpret_ifeval_mixture.R` accepts arbitrary
  `H_FIXED` and scalar or vector `G_FIXED`.
- Full generated outputs under `results/full/` or temporary analysis folders are intentionally not committed. Paper-facing snapshots and writeups are committed.
