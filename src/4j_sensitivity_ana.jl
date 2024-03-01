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

# # Prepare baseline values

R0 = 14
α = 0.05
pc = 0.25
pattern = "population_size"
ES_pattern = "ES_population_size"
n_sim = 10_000
base_kwds = (
    n_sim = n_sim,
    pattern = pattern,
    ES_pattern = ES_pattern,
    pc=pc,
)

# ## Sensitivity analysis for R0
if ARGS[1] == "1"
    println("Sensitivity analysis for R0 ================")
    summary_info = []
    for R0_sens in [10, 12, 16, 18]
        rec = run_trans_detection_process(;
            R0=R0_sens, α=α,
            base_kwds...
        )
        push!(summary_info, rec)
    end
    println(summary_info)

    println("Sensitivity analysis for α ================")
    summary_info = []
    for α_sens in [0.005, 0.010, 0.100, 0.500]
        rec = run_trans_detection_process(;
            R0=R0, α=α_sens,
            base_kwds...
        )
        push!(summary_info, rec)
    end
    println(summary_info)
end

# # Sensitivity analysis for ES related parameters

if ARGS[1] == "2"
    println("Sensitivity analysis for detection sensitivity==========================")
    summary_info = []
    path_trans = run_transmission_model(;
        R0=R0, α=α, imp_ws=[1.0],
        base_kwds...
    )
    ES_σ = 1.450
    for ES_μ in [3.121, 1.917, -0.281, -1.485]
        par_AFP, par_ES = set_par_AFP_ES(
            pc=pc, pattern=pattern, ES_pattern=ES_pattern,
            ES_μ=ES_μ, ES_σ=ES_σ
        )
        path_save = save_sensitivity_ES_catchment_area(
            par_AFP, par_ES, path_trans;
            inc_prop=0.01, pattern=pattern, ES_pattern=ES_pattern,
        )
        rec = (R0=R0, α=α, pc=pc, n_sim=n_sim,
            pattern=pattern, ES_pattern=ES_pattern,
            path=path_save, ES_μ=ES_μ, ES_σ=ES_σ)
        push!(summary_info, rec)
    end
    remove_all_transmission_results(path_trans)
    print(summary_info)

    println("Sensitivity analysis for sampling frequency==========================")
    summary_info = []
    path_trans = run_transmission_model(;
        R0=R0, α=α, imp_ws=[1.0],
        base_kwds...
    )
    for n_freq in [1, 7, 14, 60]
        par_AFP, par_ES = set_par_AFP_ES(
            pc=pc, pattern=pattern, ES_pattern=ES_pattern,
            n_freq=n_freq
        )
        path_save = save_sensitivity_ES_catchment_area(
            par_AFP, par_ES, path_trans;
            inc_prop=0.01, pattern=pattern, ES_pattern=ES_pattern,
        )
        rec = (R0=R0, α=α, pc=pc, n_sim=n_sim,
            pattern=pattern, ES_pattern=ES_pattern,
            path=path_save, n_freq=n_freq)
        push!(summary_info, rec)
    end
    remove_all_transmission_results(path_trans)
    print(summary_info)
end



