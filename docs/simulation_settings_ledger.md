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

Unless otherwise stated, each simulated data set redraws factor scores, mixture
class labels, item loadings, intercepts, and probit noise from the specified
DGP using the scenario seed.

## Main Experimental Knobs

The main simulation driver is:

```text
scripts/sample_size/compare_original_simulation_joint_mfa_gibbs.R
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

This mode matches the standalone Viroli smoke-test script.

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

Primary loadings are positive, matching
`scripts/sample_size/test_viroli_probit_independent_gibbs.R`.  Cross-loadings
remain randomly signed.

## Intercept Modes

### `INTERCEPT_MODE=none`

All item intercepts are zero.

### `INTERCEPT_MODE=viroli_smoke`

This matches the standalone Viroli smoke-test script:

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
the stricter setting used to reproduce the standalone Viroli smoke test.  It is
the better diagnostic when the question is whether the whole parameterization is
coherently recovered.

## Final Product/Viroli Simulation Grid

The current paper-facing launcher is:

```text
scripts/sample_size/run_final_product_viroli_simulation.R
```

Default design:

```text
n = 100, 200, 300, 400
p = 500, 1000, 2000, 4000 for product MAP
p = 500, 1000 for the Viroli Gibbs baselines
H = 5, 10, 15, 20
G = 2 or 3 on every factor coordinate
separation = 1, 2
replications = 25
loading design = Cross
block size mode = moderate_ifeval_like
intercept mode = ifeval_like
mixture parameter mode = viroli_smoke
alignment mode = loadings
```

Methods:

```text
independent_marginal_mixture:
  EM-SVD probit signal pretraining, sparse rotation, MAP refinement
  PRETRAIN_LOADING_PENALTY = 10
  ROTATION_LOADING_L1_PENALTY = 10
  LAMBDA_L1_PENALTY = 10

viroli_laplace_gibbs:
  probit-augmented independent-mixture Gibbs
  Bayesian lasso loading prior with VIROLI_LAMBDA_L1_PENALTY = 10

viroli_gaussian_gibbs:
  same Gibbs sampler with diffuse Gaussian loading prior
  VIROLI_LAMBDA_L1_PENALTY = 0
```

All methods use the same generated data for a scenario and replication.  Dense
joint mixture parameter tables are skipped when `G^H` exceeds
`MAX_JOINT_PARAMETER_K`; marginal mixture RMSEs are always recorded.

## Smoke-Test Reproduction Command

To reproduce the standalone Viroli smoke-test DGP inside the main driver, use:

```sh
OUT_DIR=results/diagnostics/main_driver_viroli_smoke_reproduction \
H_VALUES=10 \
G_VALUES=2,3 \
NP_GRID=n100p500:100:500 \
LOADING_DESIGNS=block_sparse \
BLOCK_SIZE_MODE=balanced \
REP_VALUES=1,2 \
MIXTURE_PARAM_MODE=viroli_smoke \
INTERCEPT_MODE=viroli_smoke \
LOADING_SIGN_MODE=smoke \
ALIGNMENT_MODE=loadings \
RUN_OURS=TRUE \
RUN_JOINT_MFA=FALSE \
RUN_VIROLI=TRUE \
OURS_PRETRAINING_METHOD=em_svd \
EM_SVD_ITER=30 \
ROTATION_ITER=30 \
ROTATION_N_MIX_STARTS=1 \
ROTATION_GRID_SIZE=11 \
REFINE_ITER=30 \
MIXTURE_UPDATE=map \
LAMBDA_L1_PENALTY=10 \
REFINE_MU_PRIOR_KAPPA=0.25 \
REFINE_WEIGHT_PRIOR_ALPHA=5 \
PARALLEL_OURS=TRUE \
VIROLI_ITER=2000 \
VIROLI_BURN=1000 \
VIROLI_THIN=1 \
VIROLI_LAMBDA_L1_PENALTY=10 \
VIROLI_SEED=1 \
VIROLI_NORMALIZE_EACH_DRAW=TRUE \
VIROLI_COMPUTE_PARAMETER_ESS=FALSE \
Rscript scripts/sample_size/compare_original_simulation_joint_mfa_gibbs.R
```

The original standalone smoke script is:

```text
scripts/sample_size/test_viroli_probit_independent_gibbs.R
```

It is useful as a compact diagnostic, but the main driver should be preferred
for paper-facing simulation results because it records the settings in every
output row.
