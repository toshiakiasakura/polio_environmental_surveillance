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
include("model_meta_pop.jl")
include("visualise_fig.jl")

using Test

# ## Check real values for the simulation results

path_res1 = "../dt_tmp_hpc/sens_ES_catchment_20240302_050027059.jld2"
path_res2 = "../dt_tmp_hpc/sens_ES_catchment_20240302_091756993.jld2"
path_res3 = "../dt_tmp_hpc/sens_ES_catchment_20240302_190345473.jld2"
path_res1_moz = "../dt_tmp_hpc/sens_ES_catchment_20240303_020804257.jld2"
path_res2_moz = "../dt_tmp_hpc/sens_ES_catchment_20240303_064157726.jld2"
path_res3_moz = "../dt_tmp_hpc/sens_ES_catchment_20240303_160629569.jld2"

function check_early_detection_prob(path)
@unpack ES_pattern, inc_prop, sim_res, path_trans, pars = load(path)
    df_fil, bin_labels = create_lead_time_category(sim_res)
    y, grp_50 = proportion_each_cate_by_group(df_fil, :ind_site)
    display(grp_50[1:10, :])
    #display(grp_50[25:32, :])
end

check_early_detection_prob(path_res1)
check_early_detection_prob(path_res2)
check_early_detection_prob(path_res3)

check_early_detection_prob(path_res1_moz)
check_early_detection_prob(path_res2_moz)
check_early_detection_prob(path_res3_moz)

# ## Check the empirical 8.96% value corresponds to which index.

# +
path_res1 = "../dt_tmp_res/sens_ES_catchment_20240219_145929.jld2"
sp_pars = read_spatial_params_file("ES_population_size")
per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100

pc = 0.25
ES_nat_cov = 11.3 # Old value, 8.96
ind = argmin(abs.(per_pop .- (ES_nat_cov/pc)))
@test ind == 31
# -

# ## Check 90 ES covered sites correspoinds to what ES population coverage. 

p = per_pop[90]

per_pop[90] |> println
@test per_pop[90] == 52.457756f0

sp_pars = read_spatial_params_file("ES_mozambique_imp_risk")
per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
per_pop[90] |> println
@test per_pop[90] == 34.28197242709779

# # Main figure 

["R0=14, α=0.05, pc=0.25, n_sim=500, pattern=population_size, ES_pattern=ES_population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_114746.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=airport, ES_pattern=ES_population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_115551.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=mozambique, ES_pattern=ES_population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_121622.jld2"] |> println
["R0=14, α=0.05, pc=0.25, n_sim=500, pattern=population_size, ES_pattern=ES_mozambique_imp_risk, path: ../dt_tmp_res/sens_ES_catchment_20240217_123152.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=airport, ES_pattern=ES_mozambique_imp_risk, path: ../dt_tmp_res/sens_ES_catchment_20240217_124140.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=mozambique, ES_pattern=ES_mozambique_imp_risk, path: ../dt_tmp_res/sens_ES_catchment_20240217_130352.jld2"] |> println

path_res1 = "../dt_tmp_res/sens_ES_catchment_20240219_145929.jld2"
path_res2 = "../dt_tmp_res/sens_ES_catchment_20240217_115551.jld2"
path_res3 = "../dt_tmp_res/sens_ES_catchment_20240217_111406.jld2"
path_res1_moz = "../dt_tmp_res/sens_ES_catchment_20240217_123152.jld2"
path_res2_moz = "../dt_tmp_res/sens_ES_catchment_20240217_124140.jld2"
path_res3_moz = "../dt_tmp_res/sens_ES_catchment_20240217_130352.jld2"


kwds = (x_var="site", xlim=[0,90], legend=false)
fontsize = 8 
xlabel = "Number of ES covered sites"
ylabel = "Probability of early detection (%)"
vis_kwds = (
    left_margin=0Plots.mm, right_margin=0Plots.mm,
    xlabelfontsize=fontsize, ylabelfontsize=fontsize, tickfontsize=fontsize -2, 
)
nothing

# +
add_annotation!(pl, text_) =  annotate!(pl, 5, 90, text(text_, :white, :left, 11))

pl1 = single_stacked_heatmap(path_res1; vis_kwds=vis_kwds, 
    ylabel=ylabel, kwds...)
plot!(pl1, left_margin=5Plots.mm)
#add_annotation!(pl1, "(A) Population size")
add_annotation!(pl1, "Scenario 1")

pl2 = single_stacked_heatmap(path_res2; vis_kwds=vis_kwds, kwds...)
#add_annotation!(pl2, "(B) Airport")
add_annotation!(pl2, "Scenario 2")

pl3 = single_stacked_heatmap(path_res3; vis_kwds=vis_kwds, kwds...)
#add_annotation!(pl3, "(C) Mozambique")
add_annotation!(pl3, "Scenario 3")

pl4 = single_stacked_heatmap(path_res1_moz; 
    vis_kwds=vis_kwds, ylabel=ylabel, xlabel=xlabel, kwds...)
plot!(pl4, bottom_margin=3Plots.mm)
#add_annotation!(pl4, "(D) Population size")
add_annotation!(pl4, "Scenario 4")

pl5 = single_stacked_heatmap(path_res2_moz; 
    vis_kwds=vis_kwds, xlabel=xlabel, kwds...)
#add_annotation!(pl5, "(E) Airport")
add_annotation!(pl5, "Scenario 5")

pl6 = single_stacked_heatmap(path_res3_moz; 
    vis_kwds=vis_kwds, xlabel=xlabel, kwds...)
#add_annotation!(pl6, "(F) Mozambique")
add_annotation!(pl6, "Scenario 6")

pls = [pl1, pl2, pl3, pl4, pl5, pl6]
plot(pls..., layout = @layout[a b c; d e f],  
    fmt=:png, dpi=300, size=(1000,500),)
# -
# # 50% stacked probs

include("visualise_fig.jl")

fontsize = 8
line5_set = (
    lw = [1, 2, 2, 1, 1],
    ls = [:dash, :solid, :solid, :solid, :dashdot],
    color = [:black, :grey, :black, :grey, :grey],
    xlabelfontsize=fontsize, ylabelfontsize=fontsize, tickfontsize=fontsize - 2,
)

#["R0=10, α=0.05, pc=0.25, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_114256.jld2", "R0=12, α=0.05, pc=0.25, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_115210.jld2", "R0=16, α=0.05, pc=0.25, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_120721.jld2"] |> println
path_R0_10 = "../dt_tmp_res/sens_ES_catchment_20240217_114256.jld2"
path_R0_12 = "../dt_tmp_res/sens_ES_catchment_20240217_115210.jld2"
path_R0_16 = "../dt_tmp_res/sens_ES_catchment_20240217_120721.jld2"
path_R0_18 = "../dt_tmp_res/sens_ES_catchment_20240219_145322.jld2"
path_Rs = [path_R0_10, path_R0_12, path_res1, path_R0_16, path_R0_18]

#["R0=14, α=0.005, pc=0.25, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_121844.jld2", "R0=14, α=0.01, pc=0.25, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_123143.jld2", "R0=14, α=0.1, pc=0.25, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_124506.jld2", "R0=14, α=0.5, pc=0.25, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_130002.jld2"] |> println
path_α_0005 = "../dt_tmp_res/sens_ES_catchment_20240217_121844.jld2"
path_α_001 = "../dt_tmp_res/sens_ES_catchment_20240217_123143.jld2"
path_α_01 = "../dt_tmp_res/sens_ES_catchment_20240217_124506.jld2"
path_α_05 = "../dt_tmp_res/sens_ES_catchment_20240217_130002.jld2"
path_αs = [path_α_0005, path_α_001, path_res1, path_α_01, path_α_05]

#["R0=14, α=0.05, pc=0.25, n_sim=500, pattern=population_size, n_freq=1, path: ../dt_tmp_res/sens_ES_catchment_20240217_201017.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=population_size, n_freq=7, path: ../dt_tmp_res/sens_ES_catchment_20240217_201642.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=population_size, n_freq=14, path: ../dt_tmp_res/sens_ES_catchment_20240217_202012.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=population_size, n_freq=60, path: ../dt_tmp_res/sens_ES_catchment_20240217_202145.jld2"] |> println
path_freq_1 = "../dt_tmp_res/sens_ES_catchment_20240217_201017.jld2"
path_freq_7 = "../dt_tmp_res/sens_ES_catchment_20240217_201642.jld2"
path_freq_14 = "../dt_tmp_res/sens_ES_catchment_20240217_202012.jld2"
path_freq_60 = "../dt_tmp_res/sens_ES_catchment_20240217_202145.jld2"
path_freqs = [path_freq_1, path_freq_7, path_freq_14, path_res1,  
    path_freq_60]

Any["R0=14, α=0.05, pc=0.25, n_sim=500, pattern=population_size, ES_μ=-0.994, path: ../dt_tmp_res/sens_ES_catchment_20240217_203237.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=population_size, ES_μ=-0.301, path: ../dt_tmp_res/sens_ES_catchment_20240217_203443.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=population_size, ES_μ=2.918, path: ../dt_tmp_res/sens_ES_catchment_20240217_203657.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=population_size, ES_μ=3.661, path: ../dt_tmp_res/sens_ES_catchment_20240217_203911.jld2"] |> println
path_det_10low = "../dt_tmp_res/sens_ES_catchment_20240217_203237.jld2"
path_det_5low = "../dt_tmp_res/sens_ES_catchment_20240217_203443.jld2"
path_det_5high = "../dt_tmp_res/sens_ES_catchment_20240217_203657.jld2"
path_det_10high = "../dt_tmp_res/sens_ES_catchment_20240217_203911.jld2"
path_dets = [path_det_10low, path_det_5low, path_res1, 
    path_det_5high, path_det_10high]

pl_R0 = plot_sens_adaptor(
    path_Rs, single_vis_R0!; 
    legendtitle="R0",
    line5_set...
)
display(pl_R0)

pl_α = plot_sens_adaptor(
    path_αs, single_vis_α!; 
    legendtitle="α",
    line5_set...
)
display(pl_α)

pl_freq = plot_sens_adaptor(
    path_freqs, single_vis_sampling_freq!; 
    legendtitle="Sampling freq.",
    line5_set...
)
display(pl_freq)

pl_freq = plot_sens_adaptor(
    path_dets, single_vis_ES_det!; 
    legendtitle="ES sensitivity",
    labels = [
        "10 times lower", "5 times lower", "5 times lower",
        "5 times higher", "10 times higher"
    ],
    line5_set...
)
display(pl_freq)








