# Fixed IFEval-Like Loading Simulation

This simulation is designed to test whether the product MAP estimator gives
factor and parameter recovery competitive with Viroli-style probit Gibbs
samplers, while being substantially faster in larger item and factor settings.

## Scientific Question

The target question is:

> Holding the population model fixed, does recovery improve as the number of
> items loading on each factor increases, and how does product MAP compare with
> Viroli Gibbs under the same data-generating model?

The simulation therefore fixes the population parameters within each design
cell and generates a new binary response dataset in each replication.

## Data-Generating Model

For observation `i` and item `j`,

```text
X_ij = 1{Z_ij > 0}
Z_ij = alpha_j + f_i' lambda_j + epsilon_ij
epsilon_ij ~ N(0, 1)
```

The latent factors are independent marginal mixtures:

```text
f_ih ~ sum_g pi_hg N(mu_hg, sigma_hg^2)
```

The mixture parameters are generated with `MIXTURE_PARAM_MODE=viroli_smoke`.
For `sep = 1`, the raw component parameters are:

- `G_h = 2`: weights `(0.5, 0.5)`, means `(-1, 1)`, sds `(0.55, 0.85)`.
- `G_h = 3`: weights `(0.3, 0.4, 0.3)`, means `(-1.35, 0, 1.35)`,
  sds `(0.45, 0.65, 0.45)`.

Each factor coordinate is then standardized before generating responses, so the
latent factor scale is comparable across settings.

## Fixed Loading Design

Only one loading design is used in the main simulation:

- IFEval-like unbalanced item blocks.
- The smallest primary-loading block has at least 30 items.
- Nonzero loading magnitudes are sampled from `Uniform(2, 3)`.
- Cross-loading probability is `0.05`.
- Cross-loading signs are random.
- Item intercepts use the IFEval-like intercept design.

For each `H`, a master loading matrix is generated at `p_max = 2000`. Smaller
`p` settings use nested block-wise subsets of this master matrix. Thus increasing
`p` increases item support without redrawing the population loading matrix.

## Simulation Grid

The main grid is:

- `n in {100, 200, 400}`
- `p in {500, 1000, 1500, 2000}` for product MAP
- `p in {500, 1000}` for the Gibbs baselines
- `H in {5, 10}`
- `G_h = 2` for every factor, or `G_h = 3` for every factor
- `25` replications per design cell

The alternating or mixed `G_h` setting is intentionally excluded to keep the
interpretation focused.

## Methods

Three methods are compared.

1. Product MAP:
   - EM-SVD probit signal pretraining.
   - Riemannian rotation toward independent marginal mixtures.
   - MAP refinement with lasso-penalized loading updates.
   - Uses 18 internal workers.

2. Viroli Laplace Gibbs:
   - Probit-augmented independent-mixture Gibbs sampler.
   - Laplace loading prior with penalty 10.
   - Uses 4 internal workers.
   - Run for 2000 iterations with 1000 burn-in draws.

3. Viroli Gaussian Gibbs:
   - Same independent-mixture Gibbs sampler.
   - Diffuse Gaussian loading prior.
   - Uses 4 internal workers.
   - Run for 2000 iterations with 1000 burn-in draws.

Replications are run serially at the launcher level. Parallelization occurs
within each fitted method.

## Metrics

Each replication records:

- factor score RMSE after alignment
- raw factor score RMSE
- loading RMSE
- intercept RMSE
- marginal mixture mean RMSE
- marginal mixture variance RMSE
- marginal mixture weight RMSE
- probability RMSE
- stage-one signal and subspace diagnostics for product MAP
- end-to-end runtime in seconds
- Gibbs effective sample size summaries when available

The output also records the fixed DGP seeds:

- `loading_parameter_seed`
- `mixture_parameter_seed`
- `data_seed`

This makes it possible to verify that population parameters are fixed across
replications while the sampled datasets change.

## Reproduction

From the repository root:

```bash
Rscript scripts/sample_size/plot_fixed_ifeval_lambda_heatmaps.R
Rscript scripts/sample_size/run_fixed_ifeval_lambda_simulation.R
```

For a short smoke run:

```bash
REP_VALUES=1,2 \
N_VALUES=100 \
P_VALUES_PRODUCT=500 \
P_VALUES_GIBBS=500 \
H_VALUES=5 \
G_CONFIG_TYPES=all2 \
  RUN_LABEL=fixed_ifeval_lambda_min30_u2_3_cp0_05_h5_h10_smoke \
Rscript scripts/sample_size/run_fixed_ifeval_lambda_simulation.R
```

Main outputs:

- Raw results:
  `results/full/fixed_ifeval_lambda_min30_u2_3_cp0_05_h5_h10/comparison_results.csv`
- Chunk-level logs:
  `results/full/fixed_ifeval_lambda_min30_u2_3_cp0_05_h5_h10/chunks/`
- True loading heatmaps:
  `results/selected_plots/sample_size/fixed_ifeval_lambda_min30_u2_3_cp0_05_h5_h10/true_lambda_heatmaps/`
- True loading matrices:
  `results/selected_tables/sample_size/fixed_ifeval_lambda_min30_u2_3_cp0_05_h5_h10/true_lambda/`
