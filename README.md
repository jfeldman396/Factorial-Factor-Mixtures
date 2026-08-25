# Factorial Factor Mixtures

This repository collects the reproducible code, selected results, and writeups for factorial factor mixture analyses.

The current contents are organized around two analysis tracks:

1. `experiments/sample_size_intercepts_centered`: simulations comparing the proposed binary probit independent-mixture factor method against a correctly specified joint-mixture factor analyzer with probit augmentation and full Gibbs sampling.
2. `experiments/ifeval`: IFEval rank, component-count, and sparsity selection by held-out predictive likelihood, followed by factor interpretation/visualization for the selected mixture model.

The main implementation lives in `R/`. Runnable scripts are kept under `scripts/` and source the shared implementation so the simulation and IFEval analyses use the same fitting code.

See `REPLICATION.md` for step-by-step instructions to rerun the IFEval
analysis and the sample-size simulation. See `CODE_AUDIT.md` for the latest
static audit notes and reproducibility caveats.

## Current Simulation

The active sample-size experiment uses:

- `H in {3, 4}`
- `G in {2, 3}`
- `n in {100, 500, 1000, 2000}`
- `p in {250, 500, 1000, 2000}`
- `25` Monte Carlo repetitions per setting
- two loading designs: `few-positive-cross` and `dense-signed-cross`
- balanced and moderately IFEval-like unbalanced item blocks
- item intercepts generated in an IFEval-like pattern
- centered augmented probit SVD initialization
- sampled augmented `Z` updates during pretraining
- MAP refinement for the proposed method
- 2000 Gibbs iterations for the joint-mixture comparator, with 1000 burn-in draws, only for `p in {250, 500}`

Selected design heatmaps, RMSE panels, timing summaries, and aggregate raw
result CSVs are stored in `results/selected_plots/sample_size` and
`results/selected_tables/sample_size`.

## IFEval Analysis

The IFEval analysis includes:

- missing-aware rank/component/penalty selection by held-out predictive likelihood
- tuned sparse loading penalty for the mixture model
- selected-model loading interpretation, including column-specific mixture sizes when selected by CV
- cross-loading summaries
- 3D factor visualizations
- LaTeX writeup and rendered PDF in `writeup/`

The current compact component-wise IFEval writeup is available as a rendered
PDF and as source Markdown:

- [Rendered IFEval component-wise PDF](writeup/ifeval_componentwise_G3313.pdf)
- [IFEval component-wise Markdown source](writeup/ifeval_componentwise_G3313.md)

It summarizes the `H = 4`, `G = (3,3,1,3)` fit, representative solo-loading
and cross-loading items, joint LLM ability profiles, the loading heatmap,
marginal mixture fits, and factor-score visualizations.

The cleaned IFEval data used by the scripts live in `data/ifeval`.

## Remote GitHub Repository

The configured remote is:

```sh
https://github.com/jfeldman396/Factorial-Factor-Mixtures.git
```
