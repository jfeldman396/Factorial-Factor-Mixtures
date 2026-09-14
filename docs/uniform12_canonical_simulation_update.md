# Uniform(1,2) Canonical Simulation Update

This run repeats the fixed-DGP IFEval-like simulation after weakening the
nonzero loadings. The motivation is diagnostic: with primary loadings sampled
from `Uniform(2,3)`, many primary-block responses become nearly deterministic
for low/high latent classes, which can make factor scores look collapsed inside
classes. The updated design uses milder nonzero loading magnitudes while keeping
the same IFEval-like block structure.

## Design Changes

The scientific grid remains:

- `n in {100, 200, 400}`
- Product MAP: `p in {500, 1000, 1500, 2000}`
- Gibbs baselines: `p in {500, 1000}`
- `H in {5, 10}`
- `G_h = 2` for all factors or `G_h = 3` for all factors
- `25` replications per setting

The DGP changes are:

- Primary nonzero loadings are sampled from `Uniform(1,2)`.
- Cross-loading magnitudes are sampled from `Uniform(1,2)`.
- Cross-loading probability remains `0.05`.
- Mixture separation remains `sep = 2`.
- The IFEval-like unbalanced block design still enforces at least 30 primary
  items in the smallest block.

The shared shrinkage schedule is reduced by two units relative to the previous
`Uniform(2,3)` run:

- `n = 100`: loading/Laplace penalty `3`
- `n = 200`: loading/Laplace penalty `3`
- `n = 400`: loading/Laplace penalty `6`

For Viroli diffuse Gaussian Gibbs, the loading penalty is zero by construction.

## Canonicalization

All paper-facing parameter RMSEs are computed after mapping estimates to a
canonical factor scale. For each fitted factor coordinate, the marginal mixture
is transformed to have mean zero and variance one. The corresponding item
intercepts and loadings are transformed so that the probit linear predictor

```text
alpha_j + f_i' lambda_j
```

is preserved. This is essential because the latent scale and location are not
identified without a convention. Product MAP applies this final canonicalization
before scoring. Viroli Gibbs normalizes each retained posterior draw using the
same helper before posterior averaging, so the comparison is on the same scale.

## Computation Plan

The run is intentionally ordered in two phases.

First, Product MAP is run for all settings. Launcher-level chunks are serial
(`TASK_WORKERS_PRODUCT=1`), and each fit uses `18` internal workers. Those
workers are used in the conditionally independent pieces of the Product MAP
algorithm: item loading/intercept updates, factor-score updates, marginal
mixture fits, rotation starts/mixture evaluations, and EM-SVD pretraining steps
where available.

Second, Gibbs baselines are run. Independent Gibbs chunks are run in parallel
over replications/H/G cells (`TASK_WORKERS_GIBBS=4`). Inside each Gibbs fit,
`4` workers are used for conditionally independent update blocks; the Markov
chain iteration order itself remains serial.

## Reproduction

From the repository root:

```bash
zsh scripts/sample_size/run_uniform12_canonical_product_first_simulation.sh
```

The output directory is:

```text
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_6_h5_h10_canonical_productfirst_gibbsparallel/
```

The combined raw result file is:

```text
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_6_h5_h10_canonical_productfirst_gibbsparallel/comparison_results.csv
```

The run is resumable. Completed chunks are skipped when the wrapper is relaunched.
