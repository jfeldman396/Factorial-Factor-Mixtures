#!/bin/zsh

set -u

export PATH="/usr/local/bin:/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin"

REPO_ROOT="/Users/joefeldman/Documents/Deep Factor Models/factorial-factor-mixtures"
RUN_LABEL="fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_6_h5_h10_canonical_productfirst_gibbsparallel"
OUT_DIR="${REPO_ROOT}/results/full/${RUN_LABEL}"
RSCRIPT="/usr/local/bin/Rscript"

/bin/mkdir -p "${OUT_DIR}"
cd "${REPO_ROOT}"

exec > "${OUT_DIR}/launcher.log" 2>&1

echo "Uniform(1,2) fixed-DGP IFEval-like canonical simulation"
echo "Started: $(/bin/date)"
echo "Repository: ${REPO_ROOT}"
echo "Output: ${OUT_DIR}"
echo "PID: $$"
echo "Rscript: ${RSCRIPT}"

export RUN_LABEL="${RUN_LABEL}"
export OUT_DIR="${OUT_DIR}"
export RESUME_EXISTING=TRUE
export STABLE_SCENARIO_SEEDS=TRUE

export RUN_PRODUCT_MAP=TRUE
export RUN_VIROLI_LAPLACE=TRUE
export RUN_VIROLI_GAUSSIAN=TRUE

export N_VALUES=100,200,400
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

# Shared loading/Laplace penalty schedule. The launcher passes the same value to
# Product MAP pretraining, rotation, refinement, and Viroli-Laplace Gibbs.
export LAMBDA_L1_PENALTY_BY_N=100=3,200=3,400=6

# Product MAP runs one independent chunk at a time, with all available cores
# assigned to the within-fit updates.
export TASK_WORKERS_PRODUCT=1
export PRODUCT_INTERNAL_WORKERS=18

# Gibbs chunks run in parallel over independent rep/H/G cells. Each Gibbs fit
# also parallelizes its conditionally independent update blocks using 4 workers.
export TASK_WORKERS_GIBBS=4
export GIBBS_PARALLEL_P_MIN=0
export GIBBS_INTERNAL_WORKERS_SERIAL=4
export GIBBS_INTERNAL_WORKERS_PARALLEL=4

export CANONICAL_NORMALIZE_OURS=TRUE
export CANONICAL_MIN_SCALE=1e-4
export VIROLI_NORMALIZE_EACH_DRAW=TRUE

export VIROLI_ITER=2000
export VIROLI_BURN=1000
export VIROLI_THIN=1
export VIROLI_COMPUTE_PARAMETER_ESS=TRUE
export VIROLI_VERBOSE=FALSE

"${RSCRIPT}" scripts/sample_size/run_fixed_ifeval_lambda_simulation.R
status=$?

echo "Finished: $(/bin/date)"
echo "Exit status: ${status}"
exit "${status}"
