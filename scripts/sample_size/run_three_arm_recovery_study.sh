#!/bin/zsh

# Authoritative replication launcher for the three-arm recovery study.
# The run is chunk-resumable: completed method/H/G/replication chunks are
# reused when the same RUN_LABEL_PREFIX and scientific settings are supplied.

set -eu

export PATH="/usr/local/bin:/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin"

SCRIPT_DIR="${0:A:h}"
REPO_ROOT="${SCRIPT_DIR:h:h}"
RSCRIPT="${RSCRIPT:-$(command -v Rscript)}"

cd "${REPO_ROOT}"

RUN_LABEL_PREFIX="${RUN_LABEL_PREFIX:-fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_full}"
ARM_FILTER="${ARM_FILTER:-separated,asymmetric_pi,moderate_overlap}"
RUN_FITS="${RUN_FITS:-TRUE}"
RUN_PLOTS="${RUN_PLOTS:-TRUE}"
VALIDATE_RESULTS="${VALIDATE_RESULTS:-TRUE}"
export RUN_LABEL_PREFIX ARM_FILTER

N_VALUES="${N_VALUES:-100,200,400}"
P_VALUES_PRODUCT="${P_VALUES_PRODUCT:-500,1000,1500,2000}"
P_VALUES_GIBBS="${P_VALUES_GIBBS:-500,1000}"
H_VALUES="${H_VALUES:-5,10}"
G_CONFIG_TYPES="${G_CONFIG_TYPES:-all2,all3}"
REP_VALUES="${REP_VALUES:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}"

RUN_PRODUCT_MAP="${RUN_PRODUCT_MAP:-TRUE}"
RUN_VIROLI_LAPLACE="${RUN_VIROLI_LAPLACE:-TRUE}"
RUN_VIROLI_GAUSSIAN="${RUN_VIROLI_GAUSSIAN:-FALSE}"

TASK_WORKERS_PRODUCT="${TASK_WORKERS_PRODUCT:-1}"
PRODUCT_INTERNAL_WORKERS="${PRODUCT_INTERNAL_WORKERS:-18}"
TASK_WORKERS_GIBBS="${TASK_WORKERS_GIBBS:-4}"
GIBBS_INTERNAL_WORKERS_SERIAL="${GIBBS_INTERNAL_WORKERS_SERIAL:-4}"
GIBBS_INTERNAL_WORKERS_PARALLEL="${GIBBS_INTERNAL_WORKERS_PARALLEL:-4}"
GIBBS_PARALLEL_P_MIN="${GIBBS_PARALLEL_P_MIN:-0}"

arm_selected() {
  [[ ",${ARM_FILTER}," == *",$1,"* ]]
}

run_arm() {
  local arm_label="$1"
  local plot_label="$2"
  local g2_pi="$3"
  local g2_mu_multiplier="$4"
  local g2_sd="$5"
  local g3_pi="$6"
  local g3_mu_multiplier="$7"
  local g3_sd="$8"

  if ! arm_selected "${arm_label}"; then
    return 0
  fi

  export RUN_LABEL="${RUN_LABEL_PREFIX}_${arm_label}"
  export OUT_DIR="${REPO_ROOT}/results/full/${RUN_LABEL}"
  /bin/mkdir -p "${OUT_DIR}"

  export RESUME_EXISTING=TRUE
  export STABLE_SCENARIO_SEEDS=TRUE
  export RUN_PRODUCT_MAP RUN_VIROLI_LAPLACE RUN_VIROLI_GAUSSIAN

  export N_VALUES P_VALUES_PRODUCT P_VALUES_GIBBS H_VALUES G_CONFIG_TYPES REP_VALUES
  export LAMBDA_L1_PENALTY_BY_N="${LAMBDA_L1_PENALTY_BY_N:-100=3,200=5,400=5}"

  export PRIMARY_LOADING_RANGE="${PRIMARY_LOADING_RANGE:-1,2}"
  export CROSS_LOADING_RANGE="${CROSS_LOADING_RANGE:-1,2}"
  export CROSS_LOADING_PROB="${CROSS_LOADING_PROB:-0.05}"
  export CROSS_SIGN_MODE="${CROSS_SIGN_MODE:-random}"
  export SEPARATIONS="${SEPARATIONS:-2}"

  export MIXTURE_PARAM_MODE=viroli_smoke
  export MIXTURE_VARIANCE_MODE=unequal
  export VIROLI_SMOKE_G2_PI="${g2_pi}"
  export VIROLI_SMOKE_G2_MU_MULTIPLIER="${g2_mu_multiplier}"
  export VIROLI_SMOKE_G2_SD="${g2_sd}"
  export VIROLI_SMOKE_G3_PI="${g3_pi}"
  export VIROLI_SMOKE_G3_MU_MULTIPLIER="${g3_mu_multiplier}"
  export VIROLI_SMOKE_G3_SD="${g3_sd}"

  # The per-n schedule supersedes these fallback values inside each task.
  export PRETRAIN_LOADING_PENALTY=5
  export ROTATION_LOADING_L1_PENALTY=5
  export LAMBDA_L1_PENALTY=5
  export VIROLI_LAPLACE_L1_PENALTY=5

  export EM_SVD_TOL_SUBSPACE="${EM_SVD_TOL_SUBSPACE:-2e-3}"
  export EM_SVD_TOL_L="${EM_SVD_TOL_L:-1e-4}"
  export EM_SVD_TOL_LOGLIK="${EM_SVD_TOL_LOGLIK:-1e-5}"

  export TASK_WORKERS_PRODUCT PRODUCT_INTERNAL_WORKERS
  export TASK_WORKERS_GIBBS GIBBS_PARALLEL_P_MIN
  export GIBBS_INTERNAL_WORKERS_SERIAL GIBBS_INTERNAL_WORKERS_PARALLEL

  export CANONICAL_NORMALIZE_OURS=TRUE
  export CANONICAL_MIN_SCALE="${CANONICAL_MIN_SCALE:-1e-4}"
  export VIROLI_NORMALIZE_EACH_DRAW=TRUE
  export REFINE_NORMALIZE_FACTOR_SCALE=FALSE
  export REFINE_NORMALIZE_FACTOR_LOCATION=TRUE
  export REFINE_FACTOR_SCORE_BOUND="${REFINE_FACTOR_SCORE_BOUND:-3}"

  export VIROLI_ITER="${VIROLI_ITER:-2000}"
  export VIROLI_BURN="${VIROLI_BURN:-1000}"
  export VIROLI_THIN="${VIROLI_THIN:-1}"
  export VIROLI_COMPUTE_PARAMETER_ESS="${VIROLI_COMPUTE_PARAMETER_ESS:-TRUE}"
  export VIROLI_VERBOSE="${VIROLI_VERBOSE:-FALSE}"

  {
    echo "Three-arm recovery study: ${plot_label}"
    echo "Run label: ${RUN_LABEL}"
    echo "Started/resumed: $(/bin/date)"
    echo "n: ${N_VALUES}"
    echo "Product MAP p: ${P_VALUES_PRODUCT}"
    echo "Gibbs p: ${P_VALUES_GIBBS}"
    echo "H: ${H_VALUES}; G types: ${G_CONFIG_TYPES}; reps: ${REP_VALUES}"
    echo "Penalty by n: ${LAMBDA_L1_PENALTY_BY_N}"
    echo "Loadings: primary=${PRIMARY_LOADING_RANGE}; cross=${CROSS_LOADING_RANGE}; cross probability=${CROSS_LOADING_PROB}"
    echo "G=2: pi=${g2_pi}; mean multiplier=${g2_mu_multiplier}; sd=${g2_sd}"
    echo "G=3: pi=${g3_pi}; mean multiplier=${g3_mu_multiplier}; sd=${g3_sd}"
    echo "Methods: Product MAP=${RUN_PRODUCT_MAP}; Gibbs Laplace=${RUN_VIROLI_LAPLACE}; Gibbs Gaussian=${RUN_VIROLI_GAUSSIAN}"
    echo "Workers: Product outer/internal=${TASK_WORKERS_PRODUCT}/${PRODUCT_INTERNAL_WORKERS}; Gibbs outer/internal=${TASK_WORKERS_GIBBS}/${GIBBS_INTERNAL_WORKERS_PARALLEL}"
  } >> "${OUT_DIR}/replication_manifest.txt"

  if [[ "${RUN_FITS}" == "TRUE" ]]; then
    echo "Launching ${plot_label}: ${RUN_LABEL}"
    "${RSCRIPT}" scripts/sample_size/run_fixed_ifeval_lambda_simulation.R \
      >> "${OUT_DIR}/launcher.log" 2>&1
  fi

  if [[ "${RUN_PLOTS}" == "TRUE" ]]; then
    local results_file="${OUT_DIR}/comparison_results.csv"
    if [[ ! -f "${results_file}" ]]; then
      echo "Skipping plots for ${plot_label}: ${results_file} does not exist."
      return 0
    fi
    export RESULTS_FILE="${results_file}"
    export PLOT_DIR="${REPO_ROOT}/results/selected_plots/sample_size/three_mixture_settings/${arm_label}"
    export TABLE_DIR="${REPO_ROOT}/results/selected_tables/sample_size/three_mixture_settings"
    export METHOD_FILTER="independent_marginal_mixture,viroli_laplace_gibbs"
    export G_COMPONENT_FILTER="2,3"
    export P_FILTER="500,1000,1500,2000"
    export OUTPUT_TAG="${arm_label}_product_map_vs_gibbs"
    export PLOT_SUBTITLE="${plot_label}; Gibbs was run for p=500/1000"
    "${RSCRIPT}" scripts/sample_size/plot_fixed_ifeval_grouped_boxplot_panels.R
  fi
}

run_arm \
  separated \
  "Separated symmetric mixtures" \
  "0.50,0.50" "1.35" "0.45,0.45" \
  "0.30,0.40,0.30" "1.35" "0.45,0.65,0.45"

run_arm \
  asymmetric_pi \
  "Asymmetric mixture probabilities" \
  "0.65,0.35" "1.35" "0.45,0.45" \
  "0.20,0.50,0.30" "1.35" "0.45,0.65,0.45"

run_arm \
  moderate_overlap \
  "Moderately overlapping mixtures" \
  "0.50,0.50" "1.10" "0.60,0.60" \
  "0.30,0.40,0.30" "1.10" "0.60,0.75,0.60"

if [[ "${VALIDATE_RESULTS}" == "TRUE" ]]; then
  "${RSCRIPT}" scripts/sample_size/validate_three_arm_recovery_study.R
fi

echo "Three-arm recovery workflow complete."
