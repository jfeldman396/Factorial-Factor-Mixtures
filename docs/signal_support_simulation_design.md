# Signal-Support Simulation Design

This is the current paper-facing simulation study.  It is meant to isolate when
the proposed product-mixture MAP estimator succeeds or fails as the binary
signal becomes easier or harder to estimate.

## Data-Generating Model

Each replication generates binary responses from

```text
X_ij = 1{Z_ij > 0}
Z_ij = alpha_j + Lambda_j' F_i + epsilon_ij
epsilon_ij iid N(0, 1).
```

The factor coordinates are independent product mixtures:

```text
F_ih ~ sum_g pi_hg N(mu_hg, sigma_hg^2).
```

The same simulated data set is used for Product MAP, Viroli Laplace Gibbs, and
Viroli diffuse Gaussian Gibbs within a scenario/replication.

## Grid

The default launcher uses:

```text
n = 100, 200
p = 500, 1000, 1500, 2000 for Product MAP
p = 500, 1000 for the Viroli Gibbs baselines
H = 5, 10, 15, 20
G = 2 or 3 on every factor coordinate
separation = 1
loading strength = weak or strong
cross-loading probability = 0.075 or 0.20
block sizes = balanced or IFEval-like unbalanced
replications = 25 per setting
```

Weak nonzero loadings are drawn from `Uniform(1.25, 1.75)`.  Strong nonzero
loadings are drawn from `Uniform(2.50, 3.00)`.  Primary loadings are positive.
Cross-loadings are randomly signed.

The factor mixture parameters use `MIXTURE_PARAM_MODE=viroli_smoke`:

```text
G = 2: pi = (0.5, 0.5), mu = (-1, 1), sd = (0.55, 0.85)
G = 3: pi = (0.30, 0.40, 0.30), mu = (-1.35, 0, 1.35), sd = (0.45, 0.65, 0.45)
```

Factors are standardized after generation, so the reported mixture parameters
are on the canonical mean-zero, unit-variance factor scale.

## Methods

`independent_marginal_mixture` is the proposed Product MAP procedure:

1. Estimate the rank-`H` centered probit signal by EM-SVD.
2. Rotate the left singular-vector scores toward independent marginal mixtures,
   with an L1 loading penalty.
3. Refine factor scores, intercepts, loadings, and marginal mixture parameters
   by MAP.

The default penalties are:

```text
PRETRAIN_LOADING_PENALTY = 10
ROTATION_LOADING_L1_PENALTY = 10
LAMBDA_L1_PENALTY = 10
```

`viroli_laplace_gibbs` is a probit-augmented independent-mixture Gibbs sampler
with the same Laplace loading penalty, run for 2000 iterations with 1000
burn-in draws.

`viroli_gaussian_gibbs` is the same Gibbs sampler with a diffuse Gaussian
loading prior, also run for 2000 iterations with 1000 burn-in draws.

## Recorded Quantities

Each result row records:

- all DGP settings and penalties;
- elapsed wall-clock seconds;
- convergence flags and completed iterations;
- factor-score RMSE and mean absolute factor correlation;
- loading RMSE;
- item-intercept RMSE;
- marginal mixture mean, variance, and weight RMSE;
- Gibbs ESS summaries when enabled;
- item-prevalence diagnostics;
- loading-support diagnostics, including cross-loading density and minimum
  loading support by factor;
- Product MAP first-stage signal and subspace diagnostics.

The first-stage diagnostics are central for interpreting this study.  They
include the relative Frobenius error of the centered probit signal and the
operator-norm sine of the principal-angle error between the estimated and true
rank-`H` signal subspaces.

## Replication

From the repository root:

```sh
Rscript scripts/sample_size/run_final_product_viroli_simulation.R
```

Generate representative loading heatmaps:

```sh
Rscript scripts/sample_size/plot_dgp_loading_heatmaps.R
```

Generate progress plots from completed chunks:

```sh
Rscript scripts/sample_size/plot_signal_support_simulation_progress.R
```

Raw full outputs are written to `results/full/signal_support_grid`.  Selected
GitHub-facing plots and tables are written to `results/selected_plots` and
`results/selected_tables`.
