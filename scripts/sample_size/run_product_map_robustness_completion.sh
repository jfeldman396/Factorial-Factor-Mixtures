#!/bin/zsh

set -u

export PATH="/usr/local/bin:/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin"

REPO_ROOT="/Users/joefeldman/Documents/Deep Factor Models/factorial-factor-mixtures"
RSCRIPT="${RSCRIPT:-$(command -v Rscript)}"
COMMON_LABEL_PREFIX="fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5_pmax2000_product_completion"

cd "${REPO_ROOT}" || exit 1

run_arm() {
  local arm_label="$1"
  local g2_pi="$2"
  local g2_mu_multiplier="$3"
  local g2_sd="$4"
  local g3_pi="$5"
  local g3_mu_multiplier="$6"
  local g3_sd="$7"

  export RUN_LABEL="${COMMON_LABEL_PREFIX}_${arm_label}"
  export OUT_DIR="${REPO_ROOT}/results/full/${RUN_LABEL}"
  /bin/mkdir -p "${OUT_DIR}"

  {
    echo "Product MAP robustness completion: ${arm_label}"
    echo "Started: $(/bin/date)"
    echo "Repository: ${REPO_ROOT}"
    echo "Output: ${OUT_DIR}"
    echo "Rscript: ${RSCRIPT}"
    echo "G2 pi=${g2_pi}; G2 multiplier=${g2_mu_multiplier}; G2 sd=${g2_sd}"
    echo "G3 pi=${g3_pi}; G3 multiplier=${g3_mu_multiplier}; G3 sd=${g3_sd}"
  } > "${OUT_DIR}/launcher.log"

  export RESUME_EXISTING=TRUE
  export STABLE_SCENARIO_SEEDS=TRUE

  export RUN_PRODUCT_MAP=TRUE
  export RUN_VIROLI_LAPLACE=FALSE
  export RUN_VIROLI_GAUSSIAN=FALSE

  # Include p=500/1000 so every n=200/400 result is generated from the same
  # p_max=2000 master loading matrix as the p=1500/2000 cells.
  export N_VALUES=200,400
  export P_VALUES_PRODUCT=500,1000,1500,2000
  export P_VALUES_GIBBS=500,1000
  export H_VALUES=5,10
  export G_CONFIG_TYPES=all2,all3
  export REP_VALUES=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25

  export PRIMARY_LOADING_RANGE=1,2
  export CROSS_LOADING_RANGE=1,2
  export CROSS_LOADING_PROB=0.05
  export CROSS_SIGN_MODE=random
  export SEPARATIONS=2

  export MIXTURE_PARAM_MODE=viroli_smoke
  export MIXTURE_VARIANCE_MODE=unequal
  export VIROLI_SMOKE_G2_PI="${g2_pi}"
  export VIROLI_SMOKE_G2_MU_MULTIPLIER="${g2_mu_multiplier}"
  export VIROLI_SMOKE_G2_SD="${g2_sd}"
  export VIROLI_SMOKE_G3_PI="${g3_pi}"
  export VIROLI_SMOKE_G3_MU_MULTIPLIER="${g3_mu_multiplier}"
  export VIROLI_SMOKE_G3_SD="${g3_sd}"

  export LAMBDA_L1_PENALTY_BY_N=200=5,400=5
  export PRETRAIN_LOADING_PENALTY=5
  export ROTATION_LOADING_L1_PENALTY=5
  export LAMBDA_L1_PENALTY=5

  export EM_SVD_TOL_SUBSPACE=2e-3
  export EM_SVD_TOL_L=1e-4
  export EM_SVD_TOL_LOGLIK=1e-5

  export TASK_WORKERS_PRODUCT=1
  export PRODUCT_INTERNAL_WORKERS=18

  export CANONICAL_NORMALIZE_OURS=TRUE
  export CANONICAL_MIN_SCALE=1e-4
  export REFINE_NORMALIZE_FACTOR_SCALE=FALSE
  export REFINE_NORMALIZE_FACTOR_LOCATION=TRUE
  export REFINE_FACTOR_SCORE_BOUND=3

  "${RSCRIPT}" scripts/sample_size/run_fixed_ifeval_lambda_simulation.R >> "${OUT_DIR}/launcher.log" 2>&1
  local rc=$?
  {
    echo "Finished: $(/bin/date)"
    echo "Exit status: ${rc}"
  } >> "${OUT_DIR}/launcher.log"
  return "${rc}"
}

run_arm \
  "asymmetric_pi" \
  "0.65,0.35" \
  "1.35" \
  "0.45,0.45" \
  "0.20,0.50,0.30" \
  "1.35" \
  "0.45,0.65,0.45"

rc=$?
if [ "${rc}" -ne 0 ]; then
  exit "${rc}"
fi

run_arm \
  "moderate_overlap" \
  "0.50,0.50" \
  "1.10" \
  "0.60,0.60" \
  "0.30,0.40,0.30" \
  "1.10" \
  "0.60,0.75,0.60"

exit "$?"
