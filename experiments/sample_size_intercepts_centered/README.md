# Sample-Size Simulation: Probit IMFM vs Joint-Mixture Gibbs

This experiment compares:

- Proposed method: binary probit independent-mixture factor model, initialized by centered augmented-`Z` SVD and refined by MAP updates.
- Comparator: full Gibbs sampler for a joint-mixture factor model with `G^H` latent profiles and diagonal within-profile factor covariance.

Each repetition simulates a new dataset, including new latent factor draws, class labels, probit augmentation noise, binary responses, and item intercepts/loadings generated from the design.

Post-processing uses the same alignment idea for both methods: choose the best factor permutation and sign pattern to match the data-generating loadings, then compare flattened parameter vectors and parameter blocks.

The current long run is launched by:

```sh
Rscript scripts/sample_size/run_sampledZ_full_pgrid_smallp_gibbs_MAP_intercepts.R
```

Paper-facing selected plots can be regenerated from the completed checkpoint
with:

```sh
OUT_DIR=results/full/sampledZ_pgrid_balanced_crossloading_smallp_gibbs_MAP_intercepts \
Rscript scripts/sample_size/plot_dgp_loading_heatmaps.R

OUT_DIR=results/full/sampledZ_pgrid_balanced_crossloading_smallp_gibbs_MAP_intercepts \
Rscript scripts/sample_size/plot_sample_size_rmse_panels.R

OUT_DIR=results/full/sampledZ_pgrid_balanced_crossloading_smallp_gibbs_MAP_intercepts \
Rscript scripts/sample_size/plot_sample_size_timing_lines.R
```

Repeat those plotting commands with
`OUT_DIR=results/full/sampledZ_pgrid_unbalanced_crossloading_smallp_gibbs_MAP_intercepts`
for the moderately IFEval-like unbalanced block setting.

These produce the DGP Lambda heatmaps, six-panel RMSE recovery plots, and
`log(seconds)` timing plots for every loading/H/G setting.

## Parallelization

The main comparison driver is:

```sh
scripts/sample_size/compare_original_simulation_joint_mfa_gibbs.R
```

Parallel behavior is controlled by:

```sh
PARALLEL_OURS=TRUE/FALSE
PARALLEL_GIBBS=TRUE/FALSE
PARALLEL_WORKERS=<integer>
```

If `PARALLEL_WORKERS` is omitted, the shared helper uses a conservative default
based on `parallel::detectCores(logical = TRUE) - 1`.

### Proposed Product-Mixture Method

When `PARALLEL_OURS=TRUE`, the proposed method passes `parallel=TRUE` and
`workers=PARALLEL_WORKERS` into the shared pretraining and MAP-refinement
functions. The parallelized work is within a single model fit:

- marginal mixture fitting across factor coordinates during rotation and refinement;
- candidate rotation work where the rotation routine can evaluate independent starts or pairwise updates concurrently;
- itemwise probit/lasso loading regressions across item columns;
- subject-wise MAP factor-score updates across rows during refinement.

The following steps remain serial or only partly parallel:

- the outer pretraining augmentation loop;
- the outer MAP-refinement loop;
- SVD of the current augmented latent matrix;
- global objective evaluation and convergence checks.

Thus the proposed method is parallelized over conditionally independent
subproblems inside each iteration, but the iteration sequence itself is still
serial.

### Joint-Mixture Gibbs Comparator

When `PARALLEL_GIBBS=TRUE`, the Gibbs sampler parallelizes the two largest
conditionally independent Gaussian updates inside each iteration:

- subject factor draws, `eta_i | Z_i, c_i, alpha, Lambda, mixture parameters`,
  across rows `i = 1, ..., n`;
- item regression draws, `(alpha_j, lambda_j) | Z_j, eta`, across item columns
  `j = 1, ..., p`.

The following Gibbs steps remain serial:

- sampling the augmented probit matrix `Z`;
- computing joint-profile probabilities for all `K = prod_h G_h` profiles;
- sampling profile labels and mixture weights;
- sampling profile-specific mixture means and variances;
- normalization of factor location/scale;
- posterior averaging and convergence diagnostics.

These serial pieces are substantial, so parallel Gibbs can be slower than
serial Gibbs when worker start-up and data-copying overhead dominate. In the
initial focused timing check with `H=3`, `G=(3,3,1)`, `p=500`, `n in {100,500}`,
and three repetitions, four workers sped up the proposed method by about `2.2x`
but slowed the Gibbs sampler.

The current recommended configuration for the full sample-size study is to
parallelize only the proposed method:

```sh
PARALLEL_OURS=TRUE \
PARALLEL_GIBBS=FALSE \
PARALLEL_WORKERS=18 \
Rscript scripts/sample_size/run_sampledZ_full_pgrid_smallp_gibbs_MAP_intercepts.R
```

This reflects the focused benchmark result: 18 workers substantially sped up
the proposed method, but slowed the Gibbs comparator because the serial MCMC
steps and worker communication overhead dominated the conditionally independent
Gaussian draws.
