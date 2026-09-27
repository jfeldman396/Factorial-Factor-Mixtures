# Saved Analysis Results

This directory contains the compact, analysis-ready results retained for the paper.

- `simulation/<arm>/recovery_results.csv`: per-replication recovery and runtime metrics for Product MAP and matched Gibbs fits.
- `simulation/<arm>/subtype_results.csv`: matched per-replication results with observation-level component-profile metrics.
- `simulation/<arm>/analysis_results.csv`: authoritative paper-facing combination of the matched rerun and large-`p` Product MAP rows.
- `simulation/dgp/`: numeric loading and mixture parameters used in the DGP figures.
- `rotation_ablation/<arm>/`: per-replication rotation results, estimated-signal diagnostics, and rank-selection results.
- `ifeval/`: final threshold-specific model-selection and fitted-model summaries from the current IFEval analysis.

The three mixture arms are `separated`, `asymmetric_pi`, and `moderate_overlap`. Resumable chunks, logs, and fitted R objects are not release artifacts and are excluded.
