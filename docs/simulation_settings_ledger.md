# Simulation Settings Ledger

This document tracks the data-generating parameters and experimental knobs used
in the sample-size simulations.  Its main purpose is to prevent accidental
comparisons across simulations that have the same short label but different
generating assumptions.

## Generative Model

All simulations use the binary probit factor model

```text
X_ij = 1{Z_ij > 0}
Z_ij = alpha_j + Lambda_j' F_i + epsilon_ij
epsilon_ij iid N(0, 1).
```

The latent factor coordinates are independent across dimensions:

```text
F_ih ~ sum_g pi_hg N(mu_hg, sigma_hg^2),
h = 1, ..., H.
```

Unless otherwise stated, each replication redraws factor scores, mixture class
labels, probit noise, and binary responses from the specified DGP using the
scenario seed.  The current fixed IFEval-like main simulation instead holds
the loading matrix, mixture parameters, and intercept design fixed within each
design cell, then redraws the latent sample and binary responses across
replications.

## Main Experimental Knobs

The shared low-level simulation engine is:

```text
scripts/sample_size/compare_original_simulation_joint_mfa_gibbs.R
```

The current paper-facing launcher that sets the main grid is:

```text
scripts/sample_size/run_fixed_ifeval_lambda_simulation.R
```

Canonical DGP utilities live in:

```text
R/sample_size_dgp.R
```

Simulation scripts should call this file for block sizes, loading designs, and
item intercepts rather than carrying local copies of those generators.

## Product-Mixture Estimator

The core product-mixture estimator should use
`OURS_PRETRAINING_METHOD=em_svd`.  This stage does not sample augmented
probit variables.  It maximizes the low-rank probit likelihood over the item
intercepts `alpha` and the rank-`H` signal matrix.  The writeup denotes this
signal by `B`; the R code calls it `L`.  The pretraining step then computes

```text
L_hat = S Lambda_svd'
```

by SVD, with `S = sqrt(n) U` the left singular-vector scores.  The next
pretraining step rotates `S` over orthogonal matrices and fits independent
one-dimensional mixtures to the rotated coordinates.  The sampled-`Z`
pretraining path is retained only as a legacy ablation and should be requested
explicitly with `OURS_PRETRAINING_METHOD=sampled_z`.

The important environment variables are:

| Variable | Role |
| --- | --- |
| `H_VALUES` | Factor dimensions to simulate. |
| `G_VALUES` or `G_CONFIGS` | Number of mixture components per factor coordinate. |
| `NP_GRID` | Sample-size and item-count grid, e.g. `n100p500:100:500`. |
| `LOADING_DESIGNS` | Loading matrix pattern. |
| `BLOCK_SIZE_MODE` | Balanced or IFEval-like block sizes. |
| `MIXTURE_PARAM_MODE` | Mixture centers, variances, and weights. |
| `MIXTURE_VARIANCE_MODE` | Equal or unequal component variances when applicable. |
| `INTERCEPT_MODE` | Item intercept generator. |
| `LOADING_SIGN_MODE` | Block-level signed loadings versus smoke-test positive primaries. |
| `ALIGNMENT_MODE` | Post-processing alignment for RMSE/correlation metrics. |
| `RUN_OURS`, `RUN_JOINT_MFA`, `RUN_VIROLI` | Which methods to fit. |
| `VIROLI_SEED` | Optional fixed random seed for Viroli.  If unset, the driver uses the scenario seed. |

Every result row from the main driver records these settings so summaries do
not average across incompatible DGPs or alignment rules.

## Mixture Parameter Modes

### `MIXTURE_PARAM_MODE=equal`

For `G_h = 2`:

```text
pi = (0.5, 0.5)
mu = (-sep, sep)
sd = (0.25, 0.65) when MIXTURE_VARIANCE_MODE=unequal
```

For `G_h = 3`:

```text
pi = (0.30, 0.40, 0.30)
mu = (-sep, 0, sep)
sd = G3_RAW_SD, default (0.25, 0.50, 0.75)
```

### `MIXTURE_PARAM_MODE=random`

Weights, means, and variances are randomly generated around the same basic
component layout, with separation controlled by `sep`.

### `MIXTURE_PARAM_MODE=viroli_smoke`

This mode preserves the compact DGP used in the earlier Viroli smoke checks.

For `G_h = 2`:

```text
pi = (0.5, 0.5)
mu = sep * (-1, 1)
sd = (0.55, 0.85)
```

For `G_h = 3`:

```text
pi = (0.30, 0.40, 0.30)
mu = sep * (-1.35, 0, 1.35)
sd = (0.45, 0.65, 0.45)
```

## Block Size Modes

### `BLOCK_SIZE_MODE=balanced`

Items are split as evenly as possible across the `H` primary loading blocks.

### `BLOCK_SIZE_MODE=ifeval_like`

For `H=3` and `H=4`, this uses the original IFEval-inspired proportions from
the simulation figures.  For larger `H`, it uses a smooth decreasing profile
rather than silently reverting to balanced blocks.

### `BLOCK_SIZE_MODE=moderate_ifeval_like`

This blends the IFEval-like proportions 50/50 with a balanced allocation.  For
larger `H`, this remains genuinely unbalanced by blending the generic
IFEval-like profile with the balanced profile.

## Loading Designs

The user-facing designs for the main sample-size study are:

```text
Sparse = balanced_moderate_few_positive_cross
Cross  = balanced_moderate_dense_signed_cross
```

The aliases `sparse`, `few_positive_cross`, `cross`, `dense_signed_cross`, and
`block_cross` are resolved by `normalize_sample_size_loading_design()` in
`R/sample_size_dgp.R`.  Historical note: `block_cross` previously meant a
two-neighbor cross-loading design in some diagnostic scripts.  In the cleaned
sample-size code, `block_cross` now resolves to the main-study Cross design.
Use `neighbor_cross` for the older two-neighbor pattern.

### `block_sparse`

Items are assigned to balanced or unbalanced primary blocks.  Each item receives
one primary loading and rare weak cross-loadings.

Current main-driver primary loading range:

```text
primary loading magnitude ~ Uniform(0.75, 1.25)
cross-loading probability = 0.035
cross-loading magnitude ~ Uniform(0.12, 0.28)
```

### `block_sparse_multisigned`

This is the older two-neighbor cross-loading design used in exploratory runs.
It is now accepted as `neighbor_cross` on the command line.

```text
primary loading magnitude ~ Uniform(0.75, 1.25)
each item cross-loads on up to two neighboring factors
cross-loading magnitude ~ Uniform(0.55, 0.95)
```

### `balanced_moderate_few_positive_cross`

Positive primary loadings with rare positive weak cross-loadings:

```text
primary loading magnitude ~ Uniform(0.75, 1.25)
cross-loading probability = 0.035
cross-loading magnitude ~ Uniform(0.12, 0.28)
```

### `balanced_moderate_dense_signed_cross`

Positive primary loadings with more frequent signed cross-loadings:

```text
primary loading magnitude ~ Uniform(0.75, 1.25)
cross-loading probability = 0.25
cross-loading magnitude ~ Uniform(0.20, 0.60)
```

## Loading Sign Modes

### `LOADING_SIGN_MODE=block`

The old block simulations use block-level signs.  Primary and structured
cross-loading signs are generated from a block sign matrix.

### `LOADING_SIGN_MODE=smoke`

Primary loadings are positive, matching the earlier compact smoke-check DGP.
Cross-loadings remain randomly signed.

## Intercept Modes

### `INTERCEPT_MODE=none`

All item intercepts are zero.

### `INTERCEPT_MODE=viroli_smoke`

This matches the earlier compact Viroli smoke-check intercept DGP:

```text
block_shift = seq(-0.65, 0.65, length.out = H)
alpha_j = block_shift[block_j] + N(0, 0.20^2)
alpha_j clipped to [-1.50, 1.50]
```

Other modes in the driver include `block`, `random`, and `ifeval_like`.

## Alignment Modes

### `ALIGNMENT_MODE=factors`

Columns are matched by factor-score correlation.  This is useful for assessing
latent-score recovery, but it can make a method look good on factor scores even
when the estimated loading matrix is poorly oriented.

### `ALIGNMENT_MODE=loadings`

Columns are matched by loading-vector distance, with sign correction.  This is
the stricter setting used in the current final simulation.  It is
the better diagnostic when the question is whether the whole parameterization is
coherently recovered.

## Fixed IFEval-Like Main Simulation Grid

The current paper-facing launcher is:

```text
scripts/sample_size/run_fixed_ifeval_lambda_simulation.R
```

The canonical run label is:

```text
fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10
```

Default design:

```text
n = 100, 200, 400
p = 500, 1000, 1500, 2000 for Product MAP
p = 500, 1000 for the Viroli Gibbs baselines
H = 5, 10
G_h = 2 for all factor coordinates, or G_h = 3 for all coordinates
separation = 2
replications = 25
loading design = Cross/IFEval-like
block size mode = ifeval_min30
smallest primary block size >= 30
primary loading magnitude = Uniform(2, 3)
cross-loading magnitude = Uniform(2, 3)
cross-loading probability = 0.05
cross-loading signs = random
intercept mode = ifeval_like
mixture parameter mode = viroli_smoke
alignment mode = loadings
```

The loading matrix and mixture parameters are fixed within each design cell.
Each replication draws a new set of factor scores, mixture classes, probit
noise, and binary responses from that fixed population model.  Smaller `p`
settings use nested block-wise subsets of a master `p = 2000` loading matrix,
so increasing `p` increases item support without changing to an unrelated
population loading design.

Methods:

```text
independent_marginal_mixture:
  EM-SVD probit signal pretraining, Riemannian rotation, MAP refinement
  PRODUCT_INTERNAL_WORKERS = 18
  PRETRAIN_LOADING_PENALTY, ROTATION_LOADING_L1_PENALTY,
  and LAMBDA_L1_PENALTY follow LAMBDA_L1_PENALTY_BY_N

viroli_laplace_gibbs:
  probit-augmented independent-mixture Gibbs
  Bayesian lasso loading prior
  VIROLI_LAMBDA_L1_PENALTY follows LAMBDA_L1_PENALTY_BY_N
  VIROLI_ITER = 2000, VIROLI_BURN = 1000

viroli_gaussian_gibbs:
  same Gibbs sampler with diffuse Gaussian loading prior
  VIROLI_LAMBDA_L1_PENALTY = 0
  VIROLI_ITER = 2000, VIROLI_BURN = 1000
```

Default n-dependent penalty schedule:

```text
LAMBDA_L1_PENALTY_BY_N = 100=5,200=5,400=8
```

All methods use the same generated data for a scenario and replication.  Dense
joint mixture parameter tables are skipped when `G^H` exceeds
`MAX_JOINT_PARAMETER_K`; marginal mixture RMSEs are always recorded.

The completed main run contains 2400 rows:

```text
1200 Product MAP rows
 600 Viroli Laplace Gibbs rows
 600 Viroli Gaussian Gibbs rows
```

Raw output, ignored by git:

```text
results/full/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10/comparison_results.csv
```

Committed selected summaries:

```text
results/selected_tables/sample_size/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10_*.csv
results/selected_plots/sample_size/fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10/
```
