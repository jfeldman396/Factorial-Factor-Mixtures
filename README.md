# Factorial Factor Mixtures

This repository collects the reproducible code, selected results, and writeups for factorial factor mixture analyses.

The current contents are organized around two analysis tracks:

1. `experiments/sample_size_intercepts_centered`: simulations comparing the proposed binary probit independent-mixture factor method against a correctly specified joint-mixture factor analyzer with probit augmentation and full Gibbs sampling.
2. `experiments/ifeval`: IFEval rank and sparsity selection by held-out predictive likelihood, predictive comparison against ordinary binary probit factor models, and factor interpretation/visualization for the selected model.

The main implementation lives in `R/`. Runnable scripts are kept under `scripts/` and source the shared implementation so the simulation and IFEval analyses use the same fitting code.

See `REPLICATION.md` for step-by-step instructions to rerun the IFEval analysis and the sample-size simulation.

## Current Simulation

The active sample-size experiment uses:

- `H in {3, 4}`
- `G in {2, 3}`
- `n in {100, 500, 1000, 2000}`
- `p = 500`
- `25` Monte Carlo repetitions per setting
- two loading designs: `few-positive-cross` and `dense-signed-cross`
- item intercepts generated in an IFEval-like pattern
- centered augmented probit SVD initialization
- MAP refinement for the proposed method
- 2000 Gibbs iterations for the joint-mixture comparator, with 1000 burn-in draws

Selected checkpoint plots and tables are stored in `results/selected_plots/sample_size` and `results/selected_tables/sample_size`.

## IFEval Analysis

The IFEval analysis includes:

- singular-value shelf diagnostics for the augmented probit matrix
- Bai-Ng `ICp2` rank diagnostics
- missing-aware rank selection by held-out probit likelihood, held-out BIC, and training BIC
- ordinary binary probit factor comparison
- tuned sparse loading penalty for the mixture model
- selected-model loading interpretation, including column-specific mixture sizes when selected by CV
- cross-loading summaries
- 3D factor visualizations
- LaTeX writeup and rendered PDF in `writeup/`

The cleaned IFEval data used by the scripts live in `data/ifeval`.

## Remote GitHub Repository

The configured remote is:

```sh
https://github.com/jfeldman396/Factorial-Factor-Mixtures.git
```
