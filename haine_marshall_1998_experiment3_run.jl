include("hm98_common.jl")
using .HM98Common

experiment = experiment_spec(3)
runtime = load_runtime(experiment)

run_case(@__DIR__, experiment, runtime)
