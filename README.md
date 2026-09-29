# Factorial Factor Mixtures

This repository is the reproducibility release for the P-IFA paper. It contains only the code, saved results, and manuscript artifacts for:

1. the three-setting Product MAP versus Gibbs recovery study;
2. the FastICA/mixture-rotation ablation; and
3. the row-wise few-shot IFEval model-selection and Gaussian comparison.

## Repository Layout

```text
R/                         Shared estimation and evaluation functions
scripts/sample_size/       Recovery-study and rotation-ablation entry points
scripts/ifeval/            IFEval data, tuning, fitting, and plotting scripts
scripts/data/              OpenEval-to-IFEval formatter
data/                      Committed IFEval analysis matrices and metadata
results/saved/             Analysis-ready per-replication results
results/selected_plots/    Figures used in the paper or supplement
results/selected_tables/   Tables used in the paper or supplement
docs/                      Simulation design and manuscript-ready LaTeX
tests/                     Gibbs alignment and subtype-recovery checks
```

Historical smoke tests, sensitivity detours, checkpoints, and presentation-only
exports are excluded from the reproducibility workflow and curated result
directories.

## Current Studies

The recovery study uses `n = {100, 200, 400}`, Product MAP `p = {500, 1000, 1500, 2000}`, Gibbs `p = {500, 1000}`, `H = {5, 10}`, homogeneous `G_h = 2` or `G_h = 3`, and 25 replications. It evaluates separated symmetric mixtures, asymmetric weights, and moderate overlap. The loading penalty is 3 at `n = 100` and 5 at `n = {200, 400}` for both Product MAP and Laplace Gibbs.

The rotation ablation compares FastICA, FastICA followed by mixture rotation, and identity initialization followed by mixture rotation under both estimated and oracle latent Gaussian signals.

The primary IFEval study uses the strict threshold-1 binary matrix and selects
rank `H = 2, ..., 8`, factor-specific component counts, and loading penalty by
five-fold row-wise few-shot predictive likelihood. Entire LLM rows are held out;
20% of each held-out row estimates its factor score and the remaining 80% is
scored. The selected rank is `H=4` with `lambda=4`. The paper interprets the
near-tied parsimonious `G=(2,2,3,1)` orientation and compares it with a matched
rank-4 Gaussian probit Gibbs model. See `docs/ifeval_rowwise_fewshot_cv.md`.

See [REPRODUCE.md](REPRODUCE.md) for exact commands and output locations.
