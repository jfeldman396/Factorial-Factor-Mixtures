#!/bin/zsh

# Refit the matched Product MAP and Viroli-Laplace Gibbs cells needed for
# observation-level component-profile recovery. Historical aggregate outputs
# do not contain factor responsibilities, so these metrics require refitting.

set -eu

SCRIPT_DIR="${0:A:h}"
REPO_ROOT="${SCRIPT_DIR:h:h}"

cd "${REPO_ROOT}"

export RUN_LABEL_PREFIX="${RUN_LABEL_PREFIX:-fixed_ifeval_lambda_min30_u1_2_cp0_05_sep2_npenalty3_5_component_profile}"
export ARM_FILTER="${ARM_FILTER:-separated,asymmetric_pi,moderate_overlap}"
export N_VALUES="${N_VALUES:-100,200,400}"
export P_VALUES_PRODUCT="${P_VALUES_PRODUCT:-500,1000}"
export P_VALUES_GIBBS="${P_VALUES_GIBBS:-500,1000}"
export H_VALUES="${H_VALUES:-5,10}"
export G_CONFIG_TYPES="${G_CONFIG_TYPES:-all2,all3}"
export REP_VALUES="${REP_VALUES:-1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25}"

export RUN_PRODUCT_MAP=TRUE
export RUN_VIROLI_LAPLACE=TRUE
export RUN_VIROLI_GAUSSIAN=FALSE
export RUN_FITS=TRUE
export RUN_PLOTS=TRUE
export VALIDATE_RESULTS=FALSE
export VIROLI_NORMALIZE_EACH_DRAW=TRUE
export VIROLI_ALIGN_RETAINED_DRAWS=TRUE

zsh scripts/sample_size/run_three_arm_recovery_study.sh
