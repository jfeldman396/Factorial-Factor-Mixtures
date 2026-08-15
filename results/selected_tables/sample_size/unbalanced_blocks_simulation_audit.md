# Unbalanced-Block Sample-Size Simulation Audit

This note records the reproducibility checks for the unbalanced-block simulation
stored in `results/selected_tables/sample_size/unbalanced_blocks_comparison_results.csv`.

## Scope

- Data-generating settings: `H in {3, 4}`, `G in {2, 3}`, `p = 500`, `n in {100, 500, 1000, 2000}`, 25 replications per setting.
- Loading designs:
  - `Loadings = Sparse`: few positive cross-loadings.
  - `Loadings = Cross`: denser signed cross-loadings.
- Methods compared:
  - independent marginal-mixture probit factor model with MAP refinement.
  - joint mixture factor probit Gibbs sampler with `K = G^H` joint profiles.

## Completion Checks

- The selected results table has 1600 rows: 800 product-mixture fits and 800 joint-Gibbs fits.
- Every design-by-method cell has the expected 25 replications at each sample size.
- Product-mixture refinement converged in all 800 fits.
- Product-mixture pretraining ran for the requested 10 iterations in all 800 fits. The pretraining convergence flag is usually `FALSE` because the stage hits the iteration cap rather than the strict pretraining tolerance; this should be described as capped pretraining, not a failed fit.
- The Gibbs comparison ran 2000 iterations with 1000 retained posterior draws in all 800 fits.

## No-Truth-Leakage Check

The simulation code generates truth in `simulate_original_binary_probit()` and passes only the observed binary matrix `X`, target dimensions, tuning settings, and random seeds into the two fitting routines:

- `fit_ours()` receives `X`, `H`, `G`, and tuning values.
- `fit_joint_class_probit_mfa_gibbs()` receives `X`, `H`, `K = G^H`, prior/tuning values, and a seed.

The true factors, loadings, intercepts, and mixture parameters are used only after fitting in the evaluation functions.

## Alignment Check

Both methods use the same factor-axis post-processing before computing factor and loading metrics:

- Factors are aligned to the true factors by absolute factor-score correlations.
- Signs are chosen from the matched factor correlations.
- Loadings are then transformed by the same signed permutation.
- Joint mixture profiles are aligned after factor-axis alignment by solving a profile matching problem based on component means and log variances.

This makes the reported factor/loadings comparisons apples-to-apples across the two methods. A stricter paper-ready version should replace the current greedy factor matching with a global assignment on absolute correlations, then verify that the reported conclusions are unchanged.

## Main Diagnostic Caveats

- Gibbs component labels are aligned after posterior averaging, not relabeled draw-by-draw. That is a conservative summary for the Gibbs sampler and may worsen mixture-parameter recovery when label switching occurs.
- The Gibbs sampler typically occupies far fewer than `K = G^H` joint profiles by the final draw:
  - `K = 8`: median 5 occupied classes.
  - `K = 16`: median 6 occupied classes.
  - `K = 27`: median 5 occupied classes.
  - `K = 81`: median 4 occupied classes.
- Therefore, poor Gibbs recovery in the largest `G = 3, H = 4` settings is not caused by early stopping. It reflects the practical difficulty of exploring and estimating a very large joint mixture in this `n, p` regime.

## Summary at n = 2000

At the largest sample size, the product-mixture method retains high factor recovery and small loading/intercept error across all settings. The joint Gibbs sampler recovers factors reasonably in easier settings but degrades as `G^H` grows, especially for `H = 4, G = 3`.

The strongest evidence from this simulation is computational-statistical: the independent marginal-mixture structure can estimate factor axes and item parameters accurately without needing to populate or sample all `G^H` joint profiles.
