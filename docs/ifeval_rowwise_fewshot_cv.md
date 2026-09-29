# Few-Shot Row-Wise IFEval Cross-Validation

This is the primary IFEval model-selection workflow for the paper. It selects
the P-IFA rank `H`, factor-specific component counts
`G`, and Laplace loading penalty `lambda` for prediction on previously unseen
LLMs. It complements the older response-level cross-validation by changing the
generalization target from missing responses of observed LLMs to responses of
entirely held-out LLMs.

## Design

- Five mutually exclusive row folds are used by default.
- A training-fold P-IFA model is fitted by the existing spectral pretraining,
  mixture rotation, and blockwise MAP refinement pipeline.
- No response from a test LLM contributes to the training-fold intercepts,
  loadings, factor distribution, or mixture parameters.
- For each test LLM, 20% of its observed responses are revealed and used only
  to estimate its factor vector by MAP with population parameters fixed.
- The remaining 80% are scored by plug-in Bernoulli predictive log likelihood.
- Three independently generated reveal sets are used per test LLM by default.
- Row folds and reveal masks are generated before fitting and shared by every
  hyperparameter candidate.

The tuning grid is `H = 2:8`, `G_h in {1,2,3}` with at most one `G_h = 1`
and permutation-equivalent configurations removed, and
`lambda in {0,1,2,4,8,12}`. This produces 462 candidates and 2,310 training
fits under five-fold CV.

## Run

```sh
WORKERS=18 \
MATRIX_PATH=data/ifeval_threshold_1/openeval_ifeval_only_binary_matrix.csv \
OUT_DIR=results/full/ifeval_rowwise_fewshot_cv \
zsh scripts/ifeval/run_ifeval_rowwise_fewshot_cv.sh
```

The run is checkpointed after each parallel task batch and resumes by default.
Set `RESUME_EXISTING=FALSE` to start a clean run with different settings.

The reusable R interface is:

```r
source("R/pifa_rowwise_cv.R")

out <- cv_pifa_rowwise(
  X = X,
  H_grid = 2:8,
  G_values = 1:3,
  lambda_grid = c(0, 1, 2, 4, 8, 12),
  K = 5,
  reveal_fraction = 0.20,
  n_reveal_repeats = 3,
  seed = 1,
  parallel = TRUE,
  workers = 18,
  out_dir = "results/full/ifeval_rowwise_fewshot_cv"
)

out$best_model_spec
out$cv_results
out$final_full_data_fit
```

## Outputs

- `rowwise_cv_results.csv`: candidate-level rankings and uncertainty.
- `rowwise_cv_selected_model.csv`: selected `H`, `G`, and `lambda`.
- `rowwise_cv_fold_scores.csv`: fold-level log score, accuracy, and Brier score.
- `rowwise_cv_factor_diagnostics.csv`: every held-out factor optimization.
- `rowwise_cv_fold_assignments.rds`: fixed row-fold assignment.
- `rowwise_cv_reveal_masks.rds`: fixed reveal and prediction indices.
- `rowwise_cv_selected_full_data_fit.rds`: final full-data MAP fit.
- `rowwise_cv_task_checkpoint.rds`: resumable intermediate task results.

Model selection uses total held-out predictive log likelihood divided by the
number of evaluated prediction responses. Accuracy and Brier score are retained
as descriptive metrics only.

## Paper Analysis

The completed threshold-1 search evaluated 462 candidates and 2,310 training
fits. The predictive winner was `H = 4`, component-count multiset
`G = (1,3,3,3)`, and `lambda = 4`, with mean held-out log score `-0.30584`.
The second-ranked `H = 4`, `G = (1,2,2,3)`, `lambda = 4` model had log score
`-0.30618`; the difference of `0.00033` per response is negligible relative to
the fold uncertainty. The latter component-count multiset was attached to the
displayed factor orientation as `G = (2,2,3,1)` for the paper interpretation:

```sh
BASE_FIT=results/full/ifeval_rowwise_fewshot_cv/rowwise_cv_selected_full_data_fit.rds \
MATRIX_PATH=data/ifeval_threshold_1/openeval_ifeval_only_binary_matrix.csv \
OUT_DIR=results/full/ifeval_rowwise_fewshot_cv/alternative_H4_G2-2-3-1_fixed_orientation_lambda4 \
G_FIXED=2,2,3,1 \
LAMBDA_L1_PENALTY=4 \
WORKERS=18 \
Rscript scripts/ifeval/refit_ifeval_fixed_orientation_components.R
```

The fixed-orientation fit preserves F1--F4 while canonicalizing location,
scale, and sign. It does not rerun the unconstrained rotation, which would be
free to move the three-component slot to another coordinate.

## Gaussian Comparison

The rank-4 Gaussian probit factor model is evaluated using the exact saved row
folds and reveal masks, an identical Laplace loading penalty of 4, 2,000 Gibbs
iterations, and 1,000 burn-in draws:

```sh
PIFA_CV_DIR=results/full/ifeval_rowwise_fewshot_cv \
MATRIX_PATH=data/ifeval_threshold_1/openeval_ifeval_only_binary_matrix.csv \
OUT_DIR=results/full/ifeval_rowwise_fewshot_cv/gaussian_H4_lambda4_rowwise_fewshot_cv \
FOLD_WORKERS=5 \
Rscript scripts/ifeval/run_ifeval_gaussian_rowwise_fewshot_cv.R
```

The Gaussian model obtains log score `-0.30879`, accuracy `0.87008`, and Brier
score `0.09389`. The interpretable P-IFA configuration obtains `-0.30618`,
`0.87050`, and `0.09389`, respectively. Thus P-IFA has a small predictive
log-score advantage while accuracy and Brier score are effectively tied.
