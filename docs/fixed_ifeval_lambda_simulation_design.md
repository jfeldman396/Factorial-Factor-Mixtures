# Fixed IFEval-Like Recovery Simulation

## Purpose

This study compares Product MAP with a correctly specified Viroli-style
probit Gibbs sampler. It asks whether Product MAP gives competitive recovery
of latent scores and structural parameters while reducing end-to-end runtime,
and whether that comparison changes under asymmetric or overlapping mixtures.

The population parameters are fixed within a design cell. Each Monte Carlo
replication draws a new factor sample, latent Gaussian response, and binary
response matrix from those fixed parameters.

## Data-Generating Model

For person `i` and item `j`,

```text
X_ij = 1{Z_ij > 0}
Z_ij = alpha_j + f_i' lambda_j + epsilon_ij
epsilon_ij ~ N(0, 1)
```

Factor coordinates are independent marginal mixtures,

```text
f_ih ~ sum_g pi_hg N(mu_hg, sigma_hg^2).
```

The raw mixture moments are transformed to marginal mean zero and variance
one before data generation. Fitted results are evaluated in the same canonical
parameterization.

## Loading And Intercept Design

- IFEval-like unbalanced primary-loading blocks.
- At least 30 primary items in the smallest block.
- Primary-loading magnitudes sampled once from `Uniform(1, 2)`.
- Cross-loading probability `0.05`.
- Cross-loading magnitudes sampled once from `Uniform(1, 2)`.
- Random cross-loading signs and block-level primary signs.
- IFEval-like item-intercept design.
- One master loading matrix at `p_max = 2000` for each fixed population cell.
- Smaller `p` values use nested block-wise subsets of that master matrix.

The nested construction is important: changes across `p` reflect additional
items under the same loading design rather than independently redrawn loading
matrices.

## Mixture Arms

All arms use `sep = 2`. The mean multiplier acts on base locations `(-2, 2)`
for `G=2` and `(-2, 0, 2)` for `G=3`, before canonicalization.

| Arm | G | Weights | Mean multiplier | Standard deviations |
|---|---:|---|---:|---|
| separated | 2 | `(0.50, 0.50)` | `1.35` | `(0.45, 0.45)` |
| separated | 3 | `(0.30, 0.40, 0.30)` | `1.35` | `(0.45, 0.65, 0.45)` |
| asymmetric | 2 | `(0.65, 0.35)` | `1.35` | `(0.45, 0.45)` |
| asymmetric | 3 | `(0.20, 0.50, 0.30)` | `1.35` | `(0.45, 0.65, 0.45)` |
| overlap | 2 | `(0.50, 0.50)` | `1.10` | `(0.60, 0.60)` |
| overlap | 3 | `(0.30, 0.40, 0.30)` | `1.10` | `(0.60, 0.75, 0.60)` |

## Simulation Grid

- `n in {100, 200, 400}`.
- Product MAP: `p in {500, 1000, 1500, 2000}`.
- Gibbs: `p in {500, 1000}`.
- `H in {5, 10}`.
- Either `G_h = 2` for every coordinate or `G_h = 3` for every coordinate.
- 25 replications per method and design cell.
- Shared loading penalty `lambda = 3` at `n = 100`.
- Shared loading penalty `lambda = 5` at `n in {200, 400}`.

The Product MAP grid contains 1200 rows per mixture arm. The Laplace-Gibbs
grid contains 600 rows per arm. The default two-method study therefore expects
1800 unique result rows per arm. Enabling Gaussian Gibbs adds 600 rows.

## Estimators

### Product MAP

1. Estimate the rank-`H` probit signal by EM-SVD.
2. Stop pretraining when the left singular subspace change is below `2e-3`,
   the loading change is below `1e-4`, the log-likelihood change is below
   `1e-5`, or the iteration cap is reached.
3. Rotate the estimated factor subspace using the Riemannian independent-
   mixture criterion.
4. Refine factor scores, loadings, intercepts, and marginal mixture parameters
   by blockwise MAP updates.
5. Use the same per-`n` loading penalty in pretraining, rotation, and
   refinement.
6. Retain the best posterior-objective iterate while enforcing monotone
   accepted refinement steps.
7. Canonicalize the final factor parameterization before scoring.

The production worker policy is one outer Product MAP task and 18 internal
workers.

### Gibbs With Laplace Loadings

The comparator is a probit-augmented independent-mixture Gibbs sampler with a
Gaussian scale-mixture representation of the Laplace loading prior. It uses
the same per-`n` loading penalty as Product MAP, 2000 draws, 1000 burn-in draws,
and thinning one. Every retained draw is canonicalized before posterior
averaging. Before accumulation, each retained copy is also aligned to a
running posterior-mean loading reference by a globally optimal signed
permutation, restricted to factors with the same component count. The same
transformation is applied to scores, loadings, allocations, and mixture
parameters, after which component means are re-sorted. This post-processing
does not alter the Gibbs chain, and ESS is calculated from the aligned traces.

The recorded production policy permits four outer Gibbs tasks and four
workers inside each fit. Runtime tables must report this policy with wall time.

### Optional Gaussian Gibbs

Set `RUN_VIROLI_GAUSSIAN=TRUE` to add the same Gibbs sampler with a diffuse
Gaussian loading prior. It is not required for the default three-arm figures.

## Canonical Parameterization

For every factor coordinate, let the fitted marginal mixture have mean `m_h`
and standard deviation `s_h`. The canonical transformation is

```text
f*_ih = (f_ih - m_h) / s_h
lambda*_jh = s_h lambda_jh
alpha*_j = alpha_j + sum_h m_h lambda_jh
mu*_hg = (mu_hg - m_h) / s_h
sigma*^2_hg = sigma^2_hg / s_h^2.
```

This preserves `alpha_j + f_i' lambda_j` and therefore preserves every fitted
probit probability. Product MAP is transformed after fitting. Gibbs applies
the same transformation to each retained draw before averaging. Raw
location/scale parameters are diagnostic only.

## Recorded Outcomes

- factor-score RMSE after signed-permutation alignment;
- probability RMSE;
- loading RMSE;
- intercept RMSE;
- marginal mixture-mean RMSE;
- marginal mixture-variance RMSE;
- marginal mixture-weight RMSE;
- observation-level component-profile Hamming accuracy after factor and
  component-label alignment;
- exact component-profile recovery, mean marginal adjusted Rand index, and
  soft-assignment probability, Brier-score, log-loss, and entropy diagnostics;
- stage-one signal and subspace errors for Product MAP;
- end-to-end wall time;
- Gibbs parameter ESS summaries;
- convergence flags, iteration counts, and DGP support diagnostics.

For Product MAP and Viroli Gibbs, component-profile Hamming accuracy is the
mean, over observations, of the fraction of factor-wise mixture labels that
are recovered correctly. Factor axes are first aligned by the loading-based
evaluation permutation, then the mixture labels within each factor are
permuted to maximize agreement with the simulated labels. The same aligned
assignments define exact-profile accuracy and marginal adjusted Rand indices.
Product MAP and independent-mixture Gibbs additionally report the posterior
probability assigned to the true component, Brier score, log loss, and
normalized assignment entropy. Hard joint-profile Gibbs assignments support
the hard-label metrics but not these soft calibration metrics.

## Reproducibility

`STABLE_SCENARIO_SEEDS=TRUE` makes each data seed a deterministic function of
the scientific design cell, not its position in a launcher chunk. Matched
methods therefore receive the same generated dataset. The output records the
loading, mixture, and data seeds.

Run the complete, resumable study from the repository root:

```bash
zsh scripts/sample_size/run_three_arm_recovery_study.sh
```

Historical runs saved aggregate recovery metrics but not observation-level
responsibilities. To reproduce the matched `p = 500, 1000` cells and add the
component-profile metrics for Product MAP and Viroli-Laplace Gibbs, run:

```bash
zsh scripts/sample_size/run_component_profile_recovery_addendum.sh
```

This addendum retains the same three DGP arms, scientific seeds, convergence
settings, loading penalties, and 25 replications as the main study. It writes
to separate `*_component_profile_<arm>` result roots so the original outputs
remain immutable.

The three output roots are:

```text
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full_separated/
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full_asymmetric_pi/
results/full/fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full_moderate_overlap/
```

Each root contains chunk checkpoints, a combined `comparison_results.csv`, and
a `replication_manifest.txt` with the realized settings. Rerunning the command
skips completed chunks.

Recreate the eight grouped boxplots without refitting:

```bash
RUN_FITS=FALSE RUN_PLOTS=TRUE \
zsh scripts/sample_size/run_three_arm_recovery_study.sh
```

Every metric is written as both a 300-dpi PNG for slides and a vector PDF for
papers. All three arms use the same dimensions, facet order, method colors, and
labels.

For the completed historical robustness jobs, run:

```bash
zsh scripts/sample_size/finalize_three_arm_recovery_study.sh
```

This preserves the source run on every row, validates the intended 25-rep
Product MAP coverage, records the five-rep `n=200,400` Gibbs robustness scope,
and verifies all 48 paper-ready PNG/PDF exports.

Figures are written below:

```text
results/selected_plots/sample_size/three_mixture_settings/
```

Cell means are written below:

```text
results/selected_tables/sample_size/three_mixture_settings/
```

Use `ARM_FILTER=separated`, `ARM_FILTER=asymmetric_pi`, or
`ARM_FILTER=moderate_overlap` to fit or plot one arm.

## Interrupted-Run Utilities

These scripts exist to recover historical partial runs and are not the
preferred entry point for a clean replication:

- `run_n100_lambda3_study.sh`: the complete `n = 100` block for all arms.
- `run_product_map_robustness_completion.sh`: Product MAP at `n = 200, 400`
  for the asymmetric and overlap arms, with `p_max = 2000`.
- `run_lambda5_sensitivity_asym_overlap.sh`: superseded five-replication pilot.
- `assemble_three_arm_recovery_results.R` and
  `finalize_three_arm_recovery_study.sh`: provenance-preserving assembly and
  final figure export for historical partial runs.

## Rotation Ablation And Rank Diagnostic

After the recovery study, run:

```bash
zsh scripts/sample_size/run_rotation_ablation_three_arms.sh
```

This uses the same three DGP arms and full `n,p,H,G` grid. On each generated
dataset it compares FastICA only, FastICA initialized mixture rotation, and
identity/SVD initialized mixture rotation under both estimated and oracle
latent Gaussian responses. It records factor RMSE before and after identical
refinement, parameter recovery, runtime, estimated-versus-oracle signal error,
and eigengap rank diagnostics from a separate overcomplete rank-15 fit.
The loading penalty is `3` at `n = 100` and `5` at `n in {200, 400}`, and the
ablation uses 18 internal workers.

After all three arms finish, validate coverage and create tracked compact
tables and paper-ready PNG/PDF figures with:

```bash
Rscript scripts/sample_size/summarize_rotation_ablation_three_arms.R
```

The completed findings and output map are documented in
`docs/rotation_ablation_results.md`.
