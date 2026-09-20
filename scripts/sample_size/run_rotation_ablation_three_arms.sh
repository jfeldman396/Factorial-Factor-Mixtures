#!/bin/zsh

set -u

export PATH="/usr/local/bin:/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin"

REPO_ROOT="/Users/joefeldman/Documents/Deep Factor Models/factorial-factor-mixtures"
RSCRIPT="${RSCRIPT:-$(command -v Rscript)}"
cd "${REPO_ROOT}" || exit 1

COMMON_LABEL_PREFIX="${COMMON_LABEL_PREFIX:-rotation_ablation_ifeval_u1_2_cp0_05_sep2_penalty3_5}"
COMMON_N_VALUES="${COMMON_N_VALUES:-100,200,400}"
COMMON_P_VALUES="${COMMON_P_VALUES:-500,1000,1500,2000}"
COMMON_H_VALUES="${COMMON_H_VALUES:-5,10}"
COMMON_G_TYPES="${COMMON_G_TYPES:-all2,all3}"
COMMON_REP_VALUES="${COMMON_REP_VALUES:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}"

run_arm() {
  local arm_label="$1"
  local g2_pi="$2"
  local g2_mu_multiplier="$3"
  local g2_sd="$4"
  local g3_pi="$5"
  local g3_mu_multiplier="$6"
  local g3_sd="$7"

  export RUN_LABEL="${COMMON_LABEL_PREFIX}_${arm_label}"
  export OUT_DIR="${REPO_ROOT}/results/diagnostics/${RUN_LABEL}"
  /bin/mkdir -p "${OUT_DIR}"

  export RESUME_EXISTING="${RESUME_EXISTING:-TRUE}"
  export MIXTURE_SCENARIO="${arm_label}"

  export N_VALUES="${COMMON_N_VALUES}"
  export P_VALUES="${COMMON_P_VALUES}"
  export H_VALUES="${COMMON_H_VALUES}"
  export G_TYPES="${COMMON_G_TYPES}"
  export REP_VALUES="${COMMON_REP_VALUES}"

  export BLOCK_SIZE_MODES=ifeval_min30
  export LOADING_DESIGNS=balanced_moderate_dense_signed_cross
  export PRIMARY_LOADING_RANGE=1,2
  export CROSS_LOADING_RANGE=1,2
  export CROSS_LOADING_PROB=0.05
  export SEPARATION=2

  export VIROLI_SMOKE_G2_PI="${g2_pi}"
  export VIROLI_SMOKE_G2_MU_MULTIPLIER="${g2_mu_multiplier}"
  export VIROLI_SMOKE_G2_SD="${g2_sd}"
  export VIROLI_SMOKE_G3_PI="${g3_pi}"
  export VIROLI_SMOKE_G3_MU_MULTIPLIER="${g3_mu_multiplier}"
  export VIROLI_SMOKE_G3_SD="${g3_sd}"

  # Match the final recovery study: lambda=3 at n=100 and lambda=5 otherwise.
  export PENALTY_BY_N=100:3,200:5,400:5
  export PRETRAIN_LOADING_PENALTY=5
  export ROTATION_LOADING_L1_PENALTY=5
  export LAMBDA_L1_PENALTY=5

  export PARALLEL=TRUE
  export WORKERS="${WORKERS:-18}"

  export EM_SVD_ITER="${EM_SVD_ITER:-50}"
  export EM_SVD_TOL_LOGLIK=1e-5
  export EM_SVD_TOL_L=1e-4
  export EM_SVD_TOL_SUBSPACE=2e-3
  export EM_SVD_INIT=both
  export EM_SVD_INIT_Z=expectation

  export ROTATION_METHODS=fastica_only,mixture_fastica_start,mixture_identity_start
  export ROTATION_ITER="${ROTATION_ITER:-20}"
  export ROTATION_MIN_ITER=3
  export RIEMANNIAN_ROTATION_STEPS="${RIEMANNIAN_ROTATION_STEPS:-8}"
  export ROTATION_OBJECTIVE_TOLERANCE=1e-4
  export FASTICA_STARTS=1
  export FASTICA_SELECTION=first
  export N_MIX_STARTS="${N_MIX_STARTS:-3}"
  export MIXTURE_MAX_ITER="${MIXTURE_MAX_ITER:-100}"

  export REFINE_ITER="${REFINE_ITER:-50}"
  export REFINE_MIN_ITER=3
  export REFINE_OBJECTIVE_TOLERANCE=1e-3
  export LASSO_BACKEND=glmnet

  # Rank selection is computed from a deliberately overfitted estimated-Z
  # signal. This avoids baking the known true H into the eigengap diagnostic.
  export RUN_EIGENGAP_DIAGNOSTIC=TRUE
  export EIGENGAP_FIT_RANK="${EIGENGAP_FIT_RANK:-15}"
  export EIGENGAP_MAX_CANDIDATE_RANK="${EIGENGAP_MAX_CANDIDATE_RANK:-14}"

  {
    echo "Rotation ablation arm: ${arm_label}"
    echo "Started: $(/bin/date)"
    echo "Grid: n=${N_VALUES}; p=${P_VALUES}; H=${H_VALUES}; G=${G_TYPES}; reps=${REP_VALUES}"
    echo "Penalties by n: ${PENALTY_BY_N}"
    echo "G2 pi=${g2_pi}; multiplier=${g2_mu_multiplier}; sd=${g2_sd}"
    echo "G3 pi=${g3_pi}; multiplier=${g3_mu_multiplier}; sd=${g3_sd}"
    echo "Eigengap fit rank=${EIGENGAP_FIT_RANK}; candidate max=${EIGENGAP_MAX_CANDIDATE_RANK}"
  } > "${OUT_DIR}/launcher.log"

  "${RSCRIPT}" scripts/sample_size/compare_rotation_fastica_diagnostic.R >> "${OUT_DIR}/launcher.log" 2>&1
  local rc=$?
  {
    echo "Finished: $(/bin/date)"
    echo "Exit status: ${rc}"
  } >> "${OUT_DIR}/launcher.log"
  return "${rc}"
}

run_arm \
  separated \
  "0.50,0.50" 1.35 "0.45,0.45" \
  "0.30,0.40,0.30" 1.35 "0.45,0.65,0.45" || exit $?

run_arm \
  asymmetric_pi \
  "0.65,0.35" 1.35 "0.45,0.45" \
  "0.20,0.50,0.30" 1.35 "0.45,0.65,0.45" || exit $?

run_arm \
  moderate_overlap \
  "0.50,0.50" 1.10 "0.60,0.60" \
  "0.30,0.40,0.30" 1.10 "0.60,0.75,0.60" || exit $?

exit 0
