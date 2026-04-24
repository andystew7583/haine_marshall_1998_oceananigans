#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/Users/astewart/Documents/haine_marshall_1998_oceananigans"

# Override any of these before invoking the script, for example:
# HM98_EXP4_SIMULATION_DAYS=0.1 HM98_EXP4_SNAPSHOT_INTERVAL_DAYS=0.05 ./run_experiment4.sh
: "${HM98_EXP4_SIMULATION_DAYS:=10}"
: "${HM98_EXP4_SNAPSHOT_INTERVAL_DAYS:=0.05}"
: "${HM98_EXP4_INITIAL_DT:=60}"
: "${HM98_EXP4_MAX_DT:=600}"

export HM98_EXP4_SIMULATION_DAYS
export HM98_EXP4_SNAPSHOT_INTERVAL_DAYS
export HM98_EXP4_INITIAL_DT
export HM98_EXP4_MAX_DT

julia --project="${PROJECT_DIR}" "${PROJECT_DIR}/haine_marshall_1998_experiment4_run.jl"
julia --project="${PROJECT_DIR}" "${PROJECT_DIR}/haine_marshall_1998_experiment4_plot.jl"
