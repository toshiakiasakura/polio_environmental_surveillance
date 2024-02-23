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

base_kwds = (n_sim=10, R0=14.0, α=0.05, pc=0.25)

ES_pattern = "ES_population_size"
summary_info = []
for pattern in ["population_size", "airport", "mozambique"]
    rec = run_trans_detection_process(;
        pattern=pattern, ES_pattern=ES_pattern,
        base_kwds...
    )
    push!(summary_info, rec)
end
println(summary_info)

# + jupyter={"outputs_hidden": true}
ES_pattern = "ES_mozambique_imp_risk"
summary_info = []
for pattern in ["population_size", "airport", "mozambique"]
    rec = run_trans_detection_process(;
        pattern=pattern, ES_pattern=ES_pattern,
        base_kwds...
    )
    push!(summary_info, rec)
end
println(summary_info)
# -









