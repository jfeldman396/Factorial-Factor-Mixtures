# Rotation Ablation Results

## Design

The rotation ablation uses the same fixed IFEval-like loading design and three
mixture scenarios as the recovery study:

- sample sizes `n = 100, 200, 400`;
- item counts `p = 500, 1000, 1500, 2000`;
- ranks `H = 5, 10`;
- either two or three components on every factor;
- 25 independently generated data sets per design cell; and
- separated, asymmetric-weight, and moderate-overlap mixture scenarios.

For every data set, the same low-rank signal is passed to three rotation
strategies: FastICA only, FastICA followed by mixture rotation, and identity/SVD
initialization followed by mixture rotation. Each strategy is evaluated using
both the estimated latent Gaussian signal and its oracle low-rank counterpart.
Factor recovery is recorded before and after the same MAP refinement.

The loading penalty is 3 when `n = 100` and 5 when `n = 200, 400`. Stage-one
EM-SVD terminates using the log-likelihood, low-rank signal, and spectral
subspace convergence checks documented in the main simulation design.

## Coverage

All three arms completed successfully. Each arm contains 1,200 generated data
sets and 14,400 unique scientific result rows:

```text
3 n values x 4 p values x 2 H values x 2 G settings x 25 replications
  = 1,200 data sets per arm

1,200 x 2 signal sources x 3 rotation methods x 2 stages
  = 14,400 result rows per arm
```

Across the three arms, validation found 43,200 result rows, 3,600 unique data
sets, no duplicated scientific keys, and no missing core recovery metrics.

## Main Findings

FastICA followed by mixture rotation gave the best average post-refinement
factor recovery in every estimated-signal scenario. Its mean factor-score RMSE
was 0.113 in the separated arm, 0.118 under asymmetric weights, and 0.124 under
moderate overlap. The corresponding FastICA-only values were 0.132, 0.197, and
0.151. Identity-initialized mixture rotation was close, but slightly worse than
the FastICA-initialized version, at 0.118, 0.121, and 0.128.

The same ordering generally holds for loading and probability recovery. The
advantage of mixture rotation is largest in the harder high-rank and
three-component cells, especially under asymmetric weights. Oracle-signal fits
show a similar pattern, indicating that mixture-informed rotation contributes
beyond improvements in the estimated binary signal alone.

The refinement stage sharply reduces factor error for every rotation strategy.
However, the quality of the rotation remains consequential after refinement:
FastICA-only starts exhibit substantially poorer recovery in the difficult
cells, whereas both mixture rotations are more stable.

For rank selection, the raw eigenvalue-difference criterion is unreliable. Its
correct-selection rate ranges from approximately 1% to 13% after averaging
over `n`, `p`, and component counts. The relative eigengap and singular-value
ratio criteria select the true rank essentially perfectly: 100% at `H = 5`
and at least 99.7% at `H = 10` in every mixture scenario.

## Reproduction

Run the full ablation with:

```bash
zsh scripts/sample_size/run_rotation_ablation_three_arms.sh
```

After all arms finish, validate and regenerate the compact results with:

```bash
Rscript scripts/sample_size/summarize_rotation_ablation_three_arms.R
```

The full per-replication outputs are written below `results/diagnostics/` and
remain gitignored. Compact tables are tracked below:

```text
results/selected_tables/sample_size/rotation_ablation_three_arms/
```

Paper-ready PNG and PDF figures are tracked below:

```text
results/selected_plots/sample_size/rotation_ablation_three_arms/
```

The coverage table is the authoritative validation record. The cell summary
contains means, medians, standard deviations, standard errors, and replication
counts for every scientific cell. The paired contrast table reports each
mixture rotation minus FastICA only, so negative RMSE differences favor the
mixture rotation.
