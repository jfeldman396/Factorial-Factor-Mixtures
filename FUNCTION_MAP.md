# Function Map

This file is a starting point for code review. It records the current conceptual boundaries before a deeper refactor.

## Binary Probit Pretraining

Source files:

- `R/binary_probit_pretraining.R`
- `scripts/sample_size/binary_probit_pretraining_algorithm_commented.R`
- `scripts/ifeval/binary_probit_pretraining_algorithm_commented.R`

Main responsibilities:

- initialize item intercepts from empirical pass rates
- initialize augmented probit latent values `Z`
- update expected or sampled truncated-normal latent values
- center augmented `Z` when intercepts are estimated
- compute spectral factor/loadings initialization
- select rank by spectral criteria or external cross-validation wrappers
- estimate preliminary factors/loadings/intercepts

## Binary Probit Refinement

Source files:

- `R/binary_probit_refinement.R`
- `scripts/sample_size/binary_probit_refinement_algorithm_commented.R`
- `scripts/ifeval/binary_probit_refinement_algorithm_commented.R`

Main responsibilities:

- update binary probit item intercepts and loading rows by itemwise GLM/surrogate regression
- rotate factors toward independent one-dimensional mixture coordinates
- estimate marginal Gaussian mixture parameters for each factor
- update factor scores by MAP refinement
- normalize factor locations and adjust intercepts by the invariant transformation
- compute observed binary probit likelihood and posterior objectives

## Joint-Mixture Gibbs Comparator

Source file:

- `scripts/sample_size/compare_original_simulation_joint_mfa_gibbs.R`

Main responsibilities:

- generate binary probit data from independent mixture factors
- fit the proposed method with MAP refinement
- fit the comparison model whose full latent factor vector follows a joint mixture with `G^H` diagonal-covariance components
- sample augmented probit `Z`, factors, class labels, loadings, intercepts, mixture means, mixture variances, and mixture weights
- align fitted factors/loadings to truth by best permutation and sign
- compute parameter correlations, RMSEs, predictive RMSE, and runtime

## IFEval Scripts

Directory:

- `scripts/ifeval`

Main responsibilities:

- reproduce rank-selection and held-out predictive comparisons
- tune the sparse loading penalty
- fit and interpret the selected `H = 3, G = 3` model
- compare ordinary binary probit and mixture binary probit fits
- summarize loadings, cross-loadings, item examples, and factor profiles
- produce 2D/3D factor visualizations
