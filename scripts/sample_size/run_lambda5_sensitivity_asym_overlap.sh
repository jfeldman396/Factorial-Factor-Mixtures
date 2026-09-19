#!/bin/zsh

set -u

export PATH="/usr/local/bin:/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin"

REPO_ROOT="/Users/joefeldman/Documents/Deep Factor Models/factorial-factor-mixtures"
RSCRIPT="${RSCRIPT:-$(command -v Rscript)}"

cd "${REPO_ROOT}" || exit 1

COMMON_LABEL_PREFIX="fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_lambda5"

# These arms are intentionally lighter than the main production grid: they are
# robustness checks for mixture-shape perturbations, not a second full study.
SENS_N_VALUES="${SENS_N_VALUES:-100,200,400}"
SENS_P_VALUES_PRODUCT="${SENS_P_VALUES_PRODUCT:-500,1000}"
SENS_P_VALUES_GIBBS="${SENS_P_VALUES_GIBBS:-500,1000}"
SENS_H_VALUES="${SENS_H_VALUES:-5,10}"
SENS_G_CONFIG_TYPES="${SENS_G_CONFIG_TYPES:-all2,all3}"
SENS_REP_VALUES="${SENS_REP_VALUES:-1,2,3,4,5}"

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
    echo "Sensitivity arm: ${arm_label}"
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
  export RUN_VIROLI_LAPLACE=TRUE
  export RUN_VIROLI_GAUSSIAN=FALSE

  export N_VALUES="${SENS_N_VALUES}"
  export P_VALUES_PRODUCT="${SENS_P_VALUES_PRODUCT}"
  export P_VALUES_GIBBS="${SENS_P_VALUES_GIBBS}"
  export H_VALUES="${SENS_H_VALUES}"
  export G_CONFIG_TYPES="${SENS_G_CONFIG_TYPES}"
  export REP_VALUES="${SENS_REP_VALUES}"

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

  export LAMBDA_L1_PENALTY_BY_N=100=5,200=5,400=5
  export PRETRAIN_LOADING_PENALTY=5
  export ROTATION_LOADING_L1_PENALTY=5
  export LAMBDA_L1_PENALTY=5
  export VIROLI_LAPLACE_L1_PENALTY=5

  export EM_SVD_TOL_SUBSPACE=2e-3
  export EM_SVD_TOL_L=1e-4
  export EM_SVD_TOL_LOGLIK=1e-5

  export TASK_WORKERS_PRODUCT=1
  export PRODUCT_INTERNAL_WORKERS=18

  export TASK_WORKERS_GIBBS=4
  export GIBBS_PARALLEL_P_MIN=0
  export GIBBS_INTERNAL_WORKERS_SERIAL=4
  export GIBBS_INTERNAL_WORKERS_PARALLEL=4

  export CANONICAL_NORMALIZE_OURS=TRUE
  export CANONICAL_MIN_SCALE=1e-4
  export VIROLI_NORMALIZE_EACH_DRAW=TRUE

  export REFINE_NORMALIZE_FACTOR_SCALE=FALSE
  export REFINE_NORMALIZE_FACTOR_LOCATION=TRUE
  export REFINE_FACTOR_SCORE_BOUND=3

  export VIROLI_ITER=2000
  export VIROLI_BURN=1000
  export VIROLI_THIN=1
  export VIROLI_COMPUTE_PARAMETER_ESS=TRUE
  export VIROLI_VERBOSE=FALSE

  "${RSCRIPT}" scripts/sample_size/run_fixed_ifeval_lambda_simulation.R >> "${OUT_DIR}/launcher.log" 2>&1
  local rc=$?
  {
    echo "Finished: $(/bin/date)"
    echo "Exit status: ${rc}"
  } >> "${OUT_DIR}/launcher.log"
  return "${rc}"
}

run_arm \
  "asymmetric_pi_check" \
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
  "moderate_overlap_check" \
  "0.50,0.50" \
  "1.10" \
  "0.60,0.60" \
  "0.30,0.40,0.30" \
  "1.10" \
  "0.60,0.75,0.60"

exit "$?"
