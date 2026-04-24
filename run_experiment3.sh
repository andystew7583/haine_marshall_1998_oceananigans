#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/Users/astewart/Documents/haine_marshall_1998_oceananigans"

# Override any of these before invoking the script, for example:
# HM98_EXP3_SIMULATION_DAYS=2 HM98_EXP3_SNAPSHOT_INTERVAL_DAYS=0.05 ./run_experiment3.sh
: "${HM98_EXP3_SIMULATION_DAYS:=10}"
: "${HM98_EXP3_SNAPSHOT_INTERVAL_DAYS:=0.05}"
: "${HM98_EXP3_INITIAL_DT:=60}"
: "${HM98_EXP3_MAX_DT:=600}"

export HM98_EXP3_SIMULATION_DAYS
export HM98_EXP3_SNAPSHOT_INTERVAL_DAYS
export HM98_EXP3_INITIAL_DT
export HM98_EXP3_MAX_DT

julia --project="${PROJECT_DIR}" "${PROJECT_DIR}/haine_marshall_1998_experiment3_run.jl"
julia --project="${PROJECT_DIR}" "${PROJECT_DIR}/haine_marshall_1998_experiment3_plot.jl"
