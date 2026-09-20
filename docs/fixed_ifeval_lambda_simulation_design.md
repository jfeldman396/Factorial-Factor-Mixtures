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
For `sep = 2`, the raw component parameters are:

- `G_h = 2`: weights `(0.5, 0.5)`, means `(-2.7, 2.7)`, sds `(0.45, 0.45)`.
- `G_h = 3`: weights `(0.3, 0.4, 0.3)`, means `(-1.35, 0, 1.35)`,
  multiplied by `2`, and sds `(0.45, 0.65, 0.45)`.

Each factor coordinate is then standardized before generating responses, so the
latent factor scale is comparable across settings.

## Fixed Loading Design

Only one loading design is used in the main simulation:

- IFEval-like unbalanced item blocks.
- The smallest primary-loading block has at least 30 items.
- Nonzero loading magnitudes are sampled from `Uniform(1, 2)`.
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
   - The pretraining stage terminates when the left singular subspace stabilizes
     at tolerance `2e-3`, or when the usual log-likelihood/loading-change
     criteria or iteration cap are met.
   - Riemannian rotation toward independent marginal mixtures.
   - MAP refinement with lasso-penalized loading updates.
   - Refinement uses bounded factor-score updates and does not force scale
     normalization inside each MAP sweep.
   - Final estimates are put into the canonical factor parameterization before
     RMSE is computed: for each factor coordinate, the fitted marginal mixture
     has mean zero and variance one, and the corresponding `alpha` and
     `Lambda` are transformed so that `alpha_j + f_i' lambda_j` is unchanged.
   - Uses 18 internal workers.

2. Viroli Laplace Gibbs:
   - Probit-augmented independent-mixture Gibbs sampler.
   - Laplace loading prior using the same loading penalty as Product MAP:
     `lambda = 5` for all `n`, `H`, and `G` settings.
   - Each posterior draw is transformed to the same canonical factor
     parameterization before posterior averaging using the shared
     `canonical_normalize_factor_parameters()` helper also used by Product MAP.
     Both methods pass the same `CANONICAL_MIN_SCALE` constant to that helper.
   - Uses 4 internal workers.
   - Run for 2000 iterations with 1000 burn-in draws.

3. Viroli Gaussian Gibbs:
   - Same independent-mixture Gibbs sampler.
   - Diffuse Gaussian loading prior.
   - Each posterior draw is transformed to the same canonical factor
     parameterization before posterior averaging using the shared
     `canonical_normalize_factor_parameters()` helper also used by Product MAP.
     Both methods pass the same `CANONICAL_MIN_SCALE` constant to that helper.
   - Uses 4 internal workers.
   - Run for 2000 iterations with 1000 burn-in draws.

Replications are run serially at the launcher level. Parallelization occurs
within each fitted method.

## Metrics

Each replication records:

- factor score RMSE after alignment and empirical standardization
- raw factor score RMSE as a diagnostic only, not a paper-facing estimand
- loading RMSE
- intercept RMSE
- marginal mixture mean RMSE
- marginal mixture variance RMSE
- marginal mixture weight RMSE
- probability RMSE
- stage-one signal and subspace diagnostics for product MAP
- end-to-end runtime in seconds
- Gibbs effective sample size summaries when available

Parameter RMSEs are intended to be read on the canonical scale. Raw latent
scale/location parameters are not identified and should not be used for
method comparisons.

The output also records the fixed DGP seeds:

- `loading_parameter_seed`
- `mixture_parameter_seed`
- `data_seed`

This makes it possible to verify that population parameters are fixed across
replications while the sampled datasets change.

For matched method comparisons, `STABLE_SCENARIO_SEEDS=TRUE` should be used.
Then `data_seed` is a deterministic hash of the scientific design cell
(`rep`, `n`, `p`, `H`, `G`, loading design, cross-loading probability,
separation, and related DGP settings), not the row position in a launcher
chunk. This makes Gibbs-only reruns directly comparable to Product MAP chunks
for the same scenario.

## Reproduction

From the repository root:

```bash
zsh scripts/sample_size/run_full_grid_lambda5_subspace2e3_simulation.sh

RUN_LABEL=fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_subspace2e3_seedfix_full_grid \
Rscript scripts/sample_size/plot_fixed_ifeval_grouped_boxplot_panels.R
```

The default run is intentionally resumable.  Existing task chunks are skipped,
so rerunning the launcher continues from the last completed chunk rather than
restarting the grid.

For a short smoke run:

```bash
REP_VALUES=1,2 \
N_VALUES=100 \
P_VALUES_PRODUCT=500 \
P_VALUES_GIBBS=500 \
H_VALUES=5 \
G_CONFIG_TYPES=all2 \
  RUN_LABEL=fixed_ifeval_lambda_min30_u2_3_cp0_05_sep2_npenalty5_8_h5_h10_smoke \
Rscript scripts/sample_size/run_fixed_ifeval_lambda_simulation.R
```

### n=100 lambda-3 small-sample rerun

The final small-sample diagnostic reruns the `n = 100` portion of the same
fixed-DGP design with a weaker shared Laplace/loading penalty, `lambda = 3`.
This was added after the lambda sweep showed that the weaker penalty improves
small-sample factor recovery while retaining the same canonical parameter
scoring used in the main run.

Run all three `n = 100` arms with:

```bash
zsh scripts/sample_size/run_n100_lambda3_study.sh
```

The grid is:

- `n = 100`
- `p in {500, 1000, 1500, 2000}` for Product MAP
- `p in {500, 1000}` for Viroli Laplace Gibbs
- `H in {5, 10}`
- all-2 and all-3 marginal component configurations
- `25` replications per cell
- nonzero loading magnitudes sampled from `Uniform(1, 2)`
- cross-loading probability `0.05`
- EM-SVD subspace stopping tolerance `2e-3`
- `lambda = 3` for Product MAP pretraining, rotation, refinement, and
  Viroli-Laplace Gibbs

The three mixture arms are:

1. Separated/default:
   - `G_h = 2`: weights `(0.50, 0.50)`, means `(-2.7, 2.7)`,
     sds `(0.45, 0.45)`.
   - `G_h = 3`: weights `(0.30, 0.40, 0.30)`, means `(-2.7, 0, 2.7)`,
     sds `(0.45, 0.65, 0.45)`.
2. Asymmetric mixture probabilities:
   - `G_h = 2`: weights `(0.65, 0.35)`, means `(-2.7, 2.7)`,
     sds `(0.45, 0.45)`.
   - `G_h = 3`: weights `(0.20, 0.50, 0.30)`, means `(-2.7, 0, 2.7)`,
     sds `(0.45, 0.65, 0.45)`.
3. Moderate mixture overlap:
   - `G_h = 2`: weights `(0.50, 0.50)`, means `(-2.2, 2.2)`,
     sds `(0.60, 0.60)`.
   - `G_h = 3`: weights `(0.30, 0.40, 0.30)`, means `(-2.2, 0, 2.2)`,
     sds `(0.60, 0.75, 0.60)`.

The output roots are:

```text
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda3_n100_separated/
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda3_n100_asymmetric_pi/
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda3_n100_moderate_overlap/
```

The corresponding presentation boxplots are written under:

```text
results/selected_plots/sample_size/<run-label>/grouped_boxplots/
```

All three arms completed successfully with `600` rows apiece (`400` Product
MAP rows and `200` Viroli-Laplace Gibbs rows), `25` replications in every
expected cell, and no duplicate method/cell/replication keys. Compact
cell-level and headline summaries are committed as:

```text
results/selected_tables/sample_size/n100_lambda3_three_arm_method_cell_summary.csv
results/selected_tables/sample_size/n100_lambda3_three_arm_headline_summary.csv
results/selected_tables/sample_size/n100_lambda3_three_arm_integrity_counts.csv
```

The committed full-result snapshots are the top-level
`comparison_results.csv` files in the three output roots above. Chunk-level
checkpoints and logs remain local working artifacts and are not required to
reproduce the summary figures.

## Rotation Ablation and Rank Diagnostic

After the recovery comparison finishes, run the paired rotation ablation with:

```bash
zsh scripts/sample_size/run_rotation_ablation_three_arms.sh
```

The ablation uses the same three mixture scenarios, fixed IFEval-like loading
matrices, nested item subsets, and `25` replications. Its grid is:

- `n in {100, 200, 400}`
- `p in {500, 1000, 1500, 2000}`
- `H in {5, 10}`
- all-2 and all-3 marginal component configurations
- loading penalty `3` at `n = 100` and `5` at `n in {200, 400}`

For every simulated dataset, the same stage-one signal is passed to three
rotation arms:

1. FastICA only.
2. FastICA initialization followed by mixture-criterion rotation.
3. Identity/SVD initialization followed by mixture-criterion rotation.

Each arm is evaluated with both the EM-SVD estimated latent-response signal and
the oracle latent Gaussian response matrix. Recovery is recorded before and
after the identical MAP refinement, along with end-to-end runtime and the
estimated-versus-oracle signal discrepancy.

The rank diagnostic does not compute an eigengap from the usual rank-`H`
stage-one estimate, because that would build the known rank into the answer.
Instead, it separately fits an overcomplete rank-15 EM-SVD signal and examines
candidate ranks `1` through `14`. For each dataset it records:

- the normalized estimated-signal eigenvalue profile;
- raw and relative adjacent eigengaps and singular-value ratios;
- the rank selected by the largest value of each criterion;
- each criterion at the true `H`; and
- indicators for whether the selected rank equals the true `H`.

Dataset-level and setting-level eigengap outputs are written beside the
rotation results as `rotation_fastica_eigengap_dataset_results.csv` and
`rotation_fastica_eigengap_setting_summary.csv`.

Main output root:

```text
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_subspace2e3_seedfix_full_grid/
```

Within that directory:

- `comparison_results.csv`: combined raw results after phase checkpoints.
- `chunks/`: chunk-level logs and checkpoint result files.
- `task_status_*.csv`: completion status for each launcher phase.

The completed main run contains `2400` result rows:

- `1200` Product MAP rows;
- `600` Viroli Laplace Gibbs rows;
- `600` Viroli Gaussian Gibbs rows.

Interim and final plots are written by:

```text
scripts/sample_size/plot_fixed_ifeval_lambda_progress.R
```

Single-cell loading-recovery examples are written by:

```text
scripts/sample_size/plot_example_lambda_recovery.R
```

Deterministic DGP artifacts use the same run label:

```text
fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_subspace2e3_seedfix_full_grid
```

- Selected heatmaps and recovery plots live under:

```text
results/selected_plots/sample_size/
  fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_subspace2e3_seedfix_full_grid/
```

- Selected summary CSVs live under the following directory, with the same
  run-label prefix:

```text
results/selected_tables/sample_size/
```

- Full raw outputs under `results/full/` are intentionally ignored by git
  because they contain large chunk-level logs and checkpoint files.

## Sensitivity Arms

After the main run, two small robustness arms are run with the same loading
design, penalties, canonicalization, `n` grid, and `H/G` grid. To keep these as
"check-box" sensitivity analyses rather than a second full production study,
they use `p in {500, 1000}` and `5` replications. They retain Product MAP and
Viroli Laplace Gibbs, but skip the diffuse-Gaussian Gibbs baseline because the
main simulation already shows that it is dominated.

1. Asymmetric mixture probabilities:
   - `G_h = 2`: weights `(0.65, 0.35)`, means `(-2.7, 2.7)`,
     sds `(0.45, 0.45)`.
   - `G_h = 3`: weights `(0.20, 0.50, 0.30)`, means `(-2.7, 0, 2.7)`,
     sds `(0.45, 0.65, 0.45)`.

2. Moderate mixture overlap:
   - `G_h = 2`: weights `(0.50, 0.50)`, means `(-2.2, 2.2)`,
     sds `(0.60, 0.60)`.
   - `G_h = 3`: weights `(0.30, 0.40, 0.30)`, means `(-2.2, 0, 2.2)`,
     sds `(0.60, 0.75, 0.60)`.

Run both arms with:

```bash
zsh scripts/sample_size/run_lambda5_sensitivity_asym_overlap.sh
```

The sensitivity output roots are:

```text
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_asymmetric_pi_check/
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_moderate_overlap_check/
```
