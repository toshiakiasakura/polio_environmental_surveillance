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

# +
using ArchGDAL
using ColorSchemes
using FreqTables
using LinearAlgebra
using DataFrames
using Dates
using Pipe
using Plots
using PyFormattedStrings
using Rasters
using Serialization

include("util.jl")
include("model_meta_pop.jl")
# -

# ## Baseline for 5000

# +
path_res1 = "../dt_tmp_res/20230822_043508.ser" # 5000 simulations, R0=14.0, α=0.05, 
path_res2 = "../dt_tmp_res/20230822_002803.ser" # 5000 simulations, R0=14.0, α=0.05, International
path_res_unvac = "../dt_tmp_res/20230822_010159.ser" # 5000 simulations, R0=14.0, α=0.05, unvaccinated order

path_R10 = "../dt_tmp_res/20230820_024708.ser" # 5000 simulations, R0=10.0, α=0.05, 
path_R12 =  "../dt_tmp_res/20230817_090219.ser" # 5000 simulations, R0=12.0, α=0.05, 
path_R16 = "../dt_tmp_res/20230821_174121.ser" # 5000 simulations, R0=16.0, α=0.05, 

path_α0005 = "../dt_tmp_res/20230821_164605.ser" # 5000 simulations, R0=14.0, α=0.005, 
path_α001 = "../dt_tmp_res/20230821_165819.ser" # 5000 simulations, R0=14.0, α=0.01,
path_α01= "../dt_tmp_res/20230821_173031.ser" # 5000 simulations, R0=14.0, α=0.1, 
path_α05 = "../dt_tmp_res/20230821_175937.ser" # 5000 simulations, R0=14.0, α=0.5, 

path_lis_R0 = [
    path_R10, path_R12, path_res1, path_R16
]
path_lis_α = [
    path_α0005, path_α001, path_res1, path_α01, path_α05
]
# -

path = path_α05
path = path_res1
path_params, res_all1 = deserialize(path)
paths = fetch_sim_paths(path_params)
res = deserialize(paths[1])
n_sim = length(paths)
println("$n_sim simulations, R0=$(res.pars.R0), α=$(res.pars.α), ")
res.pars |> dump

# +
path_unvac = "../dt_tmp/spatial_params_agg230_unvac.ser"
sp_pars_unvac = deserialize(path_unvac)

path = "../dt_tmp/spatial_params_agg230.ser"
sp_pars = deserialize(path)
nothing
# -

path_params, res_all1 = deserialize(path_res1)
paths = fetch_sim_paths(path_params)
res = deserialize(paths[1])
n_sim = length(paths)
println("$n_sim simulations, R0=$(res.pars.R0), α=$(res.pars.α), ")
res.pars |> dump

# ## Visualise main figures

path_params, res_all1 = deserialize(path_res1)
path_params, res_all2 = deserialize(path_res2)
println(path_params)
nothing

function vis_ES_population_coverage(res_all1, res_all2)
    df_res1, df_res2, df_res3 = res_all1
    xlabel = "ES population coverage (%)"

    # ES population coverage
    per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
    sens_index = obtain_ES_sensitivity_index(sp_pars.pop, 0.01)
    # ES cachment area
    df_res = df_res3
    tab = detection_pattern_sensitivity(df_res, :ind_site)
    tab[:, :per_pop] = per_pop[sens_index]
    pl1 = vis_detection_pattern(tab, n_sim, :per_pop; 
        xlabel=xlabel, title="",
    )
    
    df_diff = leadtime_diff_sensitivity(df_res, :ind_site)
    df_diff[:, :per_pop] = per_pop[sens_index]
    pl2 = vis_leadtime_diff_sensitivity(df_diff, :per_pop;  
        xlabel=xlabel, title="",
        ylim=[-300, 300], legend=(0.9,0.3),
        )
        df_det = early_detect_prob_sensitivity(df_res, :ind_site) 

    df_det = early_detect_prob_sensitivity(df_res, :ind_site) 
    df_det[:, :per_pop] = per_pop[sens_index]
    pl3 = vis_early_detect_prob(df_det, :per_pop,
        xlabel=xlabel, title="", legend=(0.7, 0.2),
    )
    #plot!(pl3, [0, 100.0], [0,100.0], color=:black, label=:none, linestyle=:dot)
    
    df_res1, df_res2, df_res3 = res_all2
    airport_cov = per_pop[[11, 7, 62]]
    # ES population coverage
    per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
    sens_index = obtain_ES_sensitivity_index(sp_pars.pop, 0.01)
    # ES cachment area
    df_res = df_res3
    tab = detection_pattern_sensitivity(df_res, :ind_site)
    tab[:, :per_pop] = per_pop[sens_index]
    pl4 = vis_detection_pattern(tab, n_sim, :per_pop; 
        xlabel=xlabel, title="",
        legend=(0.6, 0.5),
    )
    vline!(pl4, airport_cov, color=:black, linestyle=:dot, alpha=0.7, label=:none)
    
    df_diff = leadtime_diff_sensitivity(df_res, :ind_site)
    df_diff[:, :per_pop] = per_pop[sens_index]
    pl5 = vis_leadtime_diff_sensitivity(df_diff, :per_pop;  
        xlabel=xlabel, title="",
        ylim=[-350, 350], legend=(0.9,0.3),
        )
    vline!(pl5, airport_cov, color=:black, linestyle=:dot, alpha=0.7, label=:none)

    df_det = early_detect_prob_sensitivity(df_res, :ind_site) 
    df_det[:, :per_pop] = per_pop[sens_index]
    pl6 = vis_early_detect_prob(df_det, :per_pop,
        xlabel=xlabel, title="", legend=(0.7, 0.2),
    )
    #plot!(pl6, [0, 100.0], [0,100.0], color=:black, label=:none, linestyle=:dot)
    vline!(pl6, airport_cov, color=:black, linestyle=:dot, alpha=0.7, label=:none)
    
    annotate!(pl1, (0.08, 0.95), "(A)")
    annotate!(pl2, (0.08, 0.95), "(B)")
    annotate!(pl3, (0.08, 0.95), "(C)")
    annotate!(pl4, (0.08, 0.95), "(D)")
    annotate!(pl5, (0.08, 0.95), "(E)")
    annotate!(pl6, (0.08, 0.95), "(F)")
    pls = [pl1, pl2, pl3, pl4, pl5, pl6]
    pl = plot(pls...,
        size=(1200, 300*2),  # 400 * 2
        layout=(2,3),
        dpi=300,
        left_margin=5Plots.mm, bottom_margin=5Plots.mm,
    )
    savefig("../res/fig_es_population_cov.png")
    display(pl)
end

include("model_meta_pop.jl")
vis_ES_population_coverage(res_all1, res_all2)


