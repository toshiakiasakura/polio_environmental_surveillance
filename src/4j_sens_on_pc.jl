# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.7
#   kernelspec:
#     display_name: Julia 1.8.3
#     language: julia
#     name: julia-1.8
# ---

include("utils.jl")
include("geo_ana.jl")
include("model_meta_pop.jl")

# # Run all patterns

base_kwds = (n_sim=10_000, R0=14.0, α=0.05)

# # Sensitivity analysis for pc

println("Sensitivity analysis for pc ================")
ES_pattern = "ES_population_size"
summary_info = []
pattern = ""
if ARGS[1] == "1"
    pattern = "population_size"
elseif ARGS[1] == "2"
    pattern = "airport"
elseif ARGS[1] == "3"
    pattern = "mozambique"
else
    pattern = "pass"
end

for pc_sens in [0.01, 0.05, 0.50, 1.0]
    if pattern == "pass"
        continue
    end
    rec = run_trans_detection_process(;
        pc=pc_sens,
        pattern=pattern, ES_pattern=ES_pattern,
        base_kwds...
    )
    push!(summary_info, rec)
end
println(summary_info)

println("Sensitivity analysis for pc ================")
if ARGS[1] == "4"
    pattern = "population_size"
elseif ARGS[1] == "5"
    pattern = "airport"
elseif ARGS[1] == "6"
    pattern = "mozambique"
else
    pattern = "pass"
end
ES_pattern = "ES_mozambique_imp_risk"
summary_info = []
for pc_sens in [0.01, 0.05, 0.50, 1.0]
    if pattern == "pass"
        continue
    end
    rec = run_trans_detection_process(;
        pc=pc_sens,
        pattern=pattern, ES_pattern=ES_pattern,
        base_kwds...
    )
    push!(summary_info, rec)
end
println(summary_info)




