# Reproduction Guide

Run all commands from the repository root. The code was developed with R 4.x on macOS and uses multicore processing through `parallel::mclapply`.

## Dependencies

Install the R packages used by the retained workflows:

```r
install.packages(c(
  "clue", "expm", "fastICA", "ggplot2", "glmnet", "gridExtra",
  "jsonlite"
))
```

The optional interactive IFEval plots also use Python packages `numpy`, `pandas`, `matplotlib`, and `plotly`. Rebuilding the IFEval matrices from OpenEval additionally requires `pyarrow` and `huggingface_hub`.

## 1. Recovery Study

The design is documented in `docs/fixed_ifeval_lambda_simulation_design.md`. The authoritative launcher is:

```sh
zsh scripts/sample_size/run_three_arm_recovery_study.sh
```

It runs 25 replications of the following grid:

- `n = {100, 200, 400}`;
- Product MAP `p = {500, 1000, 1500, 2000}`;
- Laplace Gibbs `p = {500, 1000}`;
- `H = {5, 10}`;
- all factors have either two or three mixture components;
- separated, asymmetric-weight, and moderate-overlap mixture settings.

The DGP uses fixed IFEval-like unbalanced loading blocks within each design cell, loading magnitudes from `Uniform(1, 2)`, cross-loading probability `0.05`, random cross-loading signs, and 25 independently simulated datasets. The shared loading penalty is 3 at `n = 100` and 5 otherwise. Product MAP uses 18 internal workers; Gibbs uses four launcher workers and four internal workers by default. Gibbs retains 1,000 draws after 1,000 burn-in iterations, canonicalizes and aligns retained draws before averaging, and records parameter ESS.

Resumable chunks and logs are written to `results/full/` and are ignored by git. Each completed arm exports its combined CSV to:

```text
results/saved/simulation/separated/recovery_results.csv
results/saved/simulation/asymmetric_pi/recovery_results.csv
results/saved/simulation/moderate_overlap/recovery_results.csv
```

Observation-level subtype metrics require retained allocation summaries and are reproduced with:

```sh
zsh scripts/sample_size/run_component_profile_recovery_addendum.sh
```

This writes `subtype_results.csv` beside each arm's recovery CSV. Validate the retained recovery grid with:

```sh
Rscript scripts/sample_size/build_simulation_release_results.R
Rscript scripts/sample_size/validate_three_arm_recovery_study.R
```

Generate every manuscript and supplement figure/table with:

```sh
Rscript scripts/sample_size/make_simulation_section_artifacts.R
```

The outputs are under:

```text
results/selected_plots/sample_size/paper_simulation_section/
results/selected_tables/sample_size/paper_simulation_section/
```

## 2. Rotation Ablation

Run the full three-setting ablation:

```sh
zsh scripts/sample_size/run_rotation_ablation_three_arms.sh
```

Per-replication results are written directly to:

```text
results/saved/rotation_ablation/{separated,asymmetric_pi,moderate_overlap}/
```

Validate coverage and regenerate the paper summaries:

```sh
Rscript scripts/sample_size/summarize_rotation_ablation_three_arms.R
Rscript scripts/sample_size/plot_rotation_ablation_estimated_z_boxplots.R
```

The figures and summary tables are under:

```text
results/selected_plots/sample_size/rotation_ablation_three_arms/
results/selected_tables/sample_size/rotation_ablation_three_arms/
```

## 3. IFEval

The committed model-by-item matrices and metadata are in:

```text
data/ifeval/
data/ifeval_threshold_0p5/
data/ifeval_threshold_0p67/
data/ifeval_threshold_1/
```

The `0p67` label uses the exact cutoff `2/3`. To rebuild these matrices from a local or downloaded OpenEval snapshot:

```sh
zsh scripts/ifeval/build_ifeval_threshold_matrices.sh
```

Set `LOCAL_SNAPSHOT_DIR` if the OpenEval snapshot is not at the default Hugging Face cache location.

Run rank/component/penalty selection, selected-model refitting, and interpretation for all three thresholds with:

```sh
OUT_BASE=results/full/ifeval_threshold_sensitivity_current \
  WORKERS=18 \
  zsh scripts/ifeval/run_ifeval_threshold_analyses.sh
```

The default search uses `H = 2, ..., 8`, factor-specific component counts in
`{1, 2, 3}` with at most one Gaussian coordinate, and loading penalties
`{0, 1, 2, 4, 8, 12}` under three-fold response-level cross-validation.

The workflow is resumable. Intermediate fits are written to `results/full/`; final analysis-ready outputs belong in `results/saved/ifeval/`, and paper artifacts belong in `results/selected_plots/ifeval/` and `results/selected_tables/ifeval/`.

To export a completed run manually, use:

```sh
IFEVAL_FULL_ROOT=results/full/ifeval_threshold_sensitivity_current \
  Rscript scripts/ifeval/export_ifeval_release.R
```

## Tests

Run the retained Gibbs identification and subtype checks with:

```sh
Rscript tests/test_viroli_draw_alignment.R
Rscript tests/test_viroli_mixture_updates_and_subtype_mode.R
```

## Saved Results Contract

`results/saved/` contains analysis-ready results only. Resumable chunks, logs, RDS objects, smoke tests, and intermediate diagnostics are deliberately excluded. `results/selected_plots/` and `results/selected_tables/` contain only artifacts used by the paper or supplement.
