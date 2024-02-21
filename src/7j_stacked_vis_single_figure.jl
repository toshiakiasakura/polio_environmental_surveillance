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

# # Data preparation

["R0=14, α=0.05, pc=0.25, n_sim=500, pattern=population_size, ES_pattern=ES_population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_114746.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=airport, ES_pattern=ES_population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_115551.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=mozambique, ES_pattern=ES_population_size, path: ../dt_tmp_res/sens_ES_catchment_20240217_121622.jld2"] |> println
["R0=14, α=0.05, pc=0.25, n_sim=500, pattern=population_size, ES_pattern=ES_mozambique_imp_risk, path: ../dt_tmp_res/sens_ES_catchment_20240217_123152.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=airport, ES_pattern=ES_mozambique_imp_risk, path: ../dt_tmp_res/sens_ES_catchment_20240217_124140.jld2", "R0=14, α=0.05, pc=0.25, n_sim=500, pattern=mozambique, ES_pattern=ES_mozambique_imp_risk, path: ../dt_tmp_res/sens_ES_catchment_20240217_130352.jld2"] |> println

path_res1 = "../dt_tmp_res/sens_ES_catchment_20240219_145929.jld2"
path_res2 = "../dt_tmp_res/sens_ES_catchment_20240219_150715.jld2"
path_res3 = "../dt_tmp_res/sens_ES_catchment_20240219_152555.jld2"
path_res1_moz = "../dt_tmp_res/sens_ES_catchment_20240217_123152.jld2"
path_res2_moz = "../dt_tmp_res/sens_ES_catchment_20240217_124140.jld2"
path_res3_moz = "../dt_tmp_res/sens_ES_catchment_20240217_130352.jld2"


fontsize = 8 
xlabel = "Number of ES covered site"
ylabel = "Proportion of detection pattern (%)"
vis_kwds = (
    left_margin=0Plots.mm, right_margin=0Plots.mm,
    xlabelfontsize=fontsize, ylabelfontsize=fontsize, tickfontsize=fontsize -2, 
)
kwds = (x_var="site", xlim=[0,90], legend=false)
nothing

include("visualise_fig.jl")

# ## For Figure 1

fontsize1 = 20
vis_kwds1 = (
    xlabelfontsize=fontsize1, ylabelfontsize=fontsize1, 
    tickfontsize=fontsize1-4,
)
pl1 = single_stacked_heatmap(path_res1; 
    vis_kwds=vis_kwds1,
    xlabel=xlabel, ylabel=ylabel, kwds...)
plot!(pl1, dpi=200, fmt=:png, 
    right_margin=5Plots.mm, left_margin=5Plots.mm,
    top_margin=5Plots.mm, size=(700, 600),
)

# ## Main Figure

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

pl3 = single_stacked_heatmap(path_res3; vis_kwds=vis_kwds, 
    x_var="site", xlim=[0,90], legend=true)
#add_annotation!(pl3, "(C) Mozambique")
add_annotation!(pl3, "Scenario 3")
plot!(pl3, right_margin=25Plots.mm)

pl4 = single_stacked_heatmap(path_res1_moz; 
    vis_kwds=vis_kwds, ylabel=ylabel, xlabel=xlabel, kwds...)
plot!(pl4, bottom_margin=3Plots.mm, )
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
# ## Main figure with a different scale

xlabel = "ES population coverage (%)"
kwds = (x_var="coverage", xlim=[0,100], legend=false)
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

pl3 = single_stacked_heatmap(path_res3; vis_kwds=vis_kwds, 
    x_var="coverage", xlim=[0,100], legend=true)
#add_annotation!(pl3, "(C) Mozambique")
add_annotation!(pl3, "Scenario 3")
plot!(pl3, right_margin=25Plots.mm)

pl4 = single_stacked_heatmap(path_res1_moz; 
    vis_kwds=vis_kwds, ylabel=ylabel, xlabel=xlabel, kwds...)
plot!(pl4, bottom_margin=3Plots.mm, )
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



# # Sensitivity results

fontsize = 8
line5_set = (
    lw = [1, 1, 2.3, 2, 1],
    ls = [:dot, :solid, :solid, :solid, :dot],
    #ls = [:solid, :solid, :solid, :solid, :solid],
    color = [:black, :black, :black, :grey, :grey],
    #color = cm[1:5],
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
path_det_10low = "../dt_tmp_res/sens_ES_catchment_20240217_203911.jld2"
path_det_5low = "../dt_tmp_res/sens_ES_catchment_20240217_203657.jld2"
path_det_5high = "../dt_tmp_res/sens_ES_catchment_20240217_203443.jld2"
path_det_10high = "../dt_tmp_res/sens_ES_catchment_20240217_203237.jld2"
path_dets = [path_det_10low, path_det_5low, path_res1, 
    path_det_5high, path_det_10high]

# +
### Define summarised plots.
pl_R0 = plot_sens_adaptor(
    path_Rs, single_vis_R0!; 
    xlabel="",
    legendtitle="R0",
    line5_set...
)
plot!(pl_R0, left_margin=5Plots.mm)

pl_α = plot_sens_adaptor(
    path_αs, single_vis_α!; 
    legendtitle="α",
    line5_set...
)

pl_freq = plot_sens_adaptor(
    path_freqs, single_vis_sampling_freq!; 
    xlabel="",
    legendtitle="Sampling freq.",
    line5_set...
)

pl_det = plot_sens_adaptor(
    path_dets, single_vis_ES_det!; 
    legendtitle="ES sensitivity",
    labels = [
        "10x lower", "5x lower", "5x lower",
        "5x higher", "10x higher"
    ],
    line5_set...
)
nothing

# +
#####  Individual stacked plot.
xlabel = "Number of ES covered sites" 
ylabel = "Proportion of detection pattern (%)"
### R0 
pl_R0_10= single_stacked_heatmap(path_R0_10; 
    ylabel=ylabel,
    vis_kwds=vis_kwds, kwds...)
add_annotation!(pl_R0_10, "R0 = 10")
plot!(pl_R0_10, left_margin=3Plots.mm)

pl_R0_18= single_stacked_heatmap(path_R0_18; vis_kwds=vis_kwds, 
    x_var="site", xlim=[0, 90], legend=true)
add_annotation!(pl_R0_18, "R0 = 18")
plot!(pl_R0_18, right_margin=25Plots.mm)

### α
pl_α_0005 = single_stacked_heatmap(path_α_0005; 
    ylabel=ylabel,
    vis_kwds=vis_kwds, xlabel=xlabel, kwds...)
add_annotation!(pl_α_0005, "α = 0.005")

pl_α_05 = single_stacked_heatmap(path_α_05; 
    vis_kwds=vis_kwds, xlabel=xlabel, kwds...)
add_annotation!(pl_α_05, "α = 0.500")
plot!(pl_α, left_margin=5Plots.mm)

### Sampling freq
pl_freq_1 = single_stacked_heatmap(path_freq_1; 
    ylabel=ylabel,
    vis_kwds=vis_kwds, kwds...)
#add_annotation!(pl_freq_1, "Sampling Freq. = 1 day")
annotate!(pl_freq_1, 10, 8, text("Sampling Freq. = 1 day", :white, :left, 11))

pl_freq_60= single_stacked_heatmap(path_α_05; vis_kwds=vis_kwds,
    x_var="site", xlim=[0, 90], legend=true)
add_annotation!(pl_freq_60, "Sampling freq. = 60 day")

### ES detection sensitivity
pl_det_10low = single_stacked_heatmap(path_det_10low; 
    ylabel=ylabel, 
    vis_kwds=vis_kwds, xlabel=xlabel, kwds...)
add_annotation!(pl_det_10low, "10x lower ES sensitivity")

pl_det_10high = single_stacked_heatmap(path_det_10high; 
    vis_kwds=vis_kwds, xlabel=xlabel, kwds...)
add_annotation!(pl_det_10high, "10x higehr ES sensitivity")

nothing
# -

### Synthesize above all.
pls = [
    pl_R0, pl_R0_10, pl_R0_18,
    pl_α, pl_α_0005, pl_α_05,
    pl_freq, pl_freq_1, pl_freq_60,
    pl_det, pl_det_10low, pl_det_10high,
]
layout = @layout [a b c; d e f; g h i; j k l]
plot(pls..., layout=layout, size=(1100,250*4),
    fmt=:png, dpi=300,
)

# # Results for pc. 

["R0=14, α=0.05, pc=0.01, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240219_180248.jld2", "R0=14, α=0.05, pc=0.05, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240219_181547.jld2", "R0=14, α=0.05, pc=0.5, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240219_182854.jld2", "R0=14, α=0.05, pc=1.0, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240219_184204.jld2"] |> println
path_pc001 = "../dt_tmp_res/sens_ES_catchment_20240219_180248.jld2"
path_pc005 = "../dt_tmp_res/sens_ES_catchment_20240219_181547.jld2"
path_pc050 = "../dt_tmp_res/sens_ES_catchment_20240219_182854.jld2"
path_pc100 = "../dt_tmp_res/sens_ES_catchment_20240219_184204.jld2"
path_pcs = [path_pc001, path_pc005, path_res1, path_pc050, path_pc100]

# +
pl_pc_site = plot_sens_adaptor(
    path_pcs, single_vis_pc_site!; 
    xlabel="Number of ES covered site",
    legendtitle="pc",
    line5_set...
)
plot!(pl_pc_site, left_margin=5Plots.mm, bottom_margin=3Plots.mm)

pl_pc_cov = pl_pc_coverage = plot_sens_adaptor(
    path_pcs, single_vis_pc_coverage!;
    xlabel="ES population coverage (%)", 
    legendtitle="pc",
    legend=:bottomright,
    line5_set...
)

plot(pl_pc_site, pl_pc_cov, 
    fmt=:png, dpi=300, size=(660, 250),
)
# -

include("visualise_fig.jl")

["R0=14, α=0.05, pc=0.01, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240219_194520.jld2", "R0=14, α=0.05, pc=0.05, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240219_195211.jld2", "R0=14, α=0.05, pc=0.5, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240219_195921.jld2", "R0=14, α=0.05, pc=1.0, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240219_200546.jld2"] |> println
path_pc001 = "../dt_tmp_res/sens_ES_catchment_20240219_194520.jld2"
path_pc005 = "../dt_tmp_res/sens_ES_catchment_20240219_195211.jld2"
path_pc050 = "../dt_tmp_res/sens_ES_catchment_20240219_195921.jld2"
path_pc100 = "../dt_tmp_res/sens_ES_catchment_20240219_200546.jld2"
path_pcs_airport = [path_pc001, path_pc005, path_res2, path_pc050, path_pc100]

["R0=14, α=0.05, pc=0.01, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240219_204802.jld2", "R0=14, α=0.05, pc=0.05, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240219_210449.jld2", "R0=14, α=0.05, pc=0.5, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240219_212306.jld2", "R0=14, α=0.05, pc=1.0, n_sim=500, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240219_214020.jld2"] |> println
path_pc001 = "../dt_tmp_res/sens_ES_catchment_20240219_204802.jld2"
path_pc005 = "../dt_tmp_res/sens_ES_catchment_20240219_210449.jld2"
path_pc050 = "../dt_tmp_res/sens_ES_catchment_20240219_212306.jld2"
path_pc100 = "../dt_tmp_res/sens_ES_catchment_20240219_214020.jld2"
path_pcs_moz = [path_pc001, path_pc005, path_res3, path_pc050, path_pc100]

# +
pl_pc_site_airport = plot_sens_adaptor(
    path_pcs_airport, single_vis_pc_site!; 
    xlabel="Number of ES covered site",
    legendtitle="pc",
    legend = :bottomright,
    line5_set...
)
plot!(pl_pc_site_airport, left_margin=5Plots.mm, bottom_margin=3Plots.mm)

pl_pc_cov_airport = pl_pc_coverage = plot_sens_adaptor(
    path_pcs_airport, single_vis_pc_coverage!;
    xlabel="ES population coverage (%)", 
    legendtitle="pc",
    legend=:bottomright,
    line5_set...
)

pl_pc_site_moz = plot_sens_adaptor(
    path_pcs_moz, single_vis_pc_site!; 
    xlabel="Number of ES covered site",
    legendtitle="pc",
    legend = :bottomright,
    line5_set...
)
plot!(pl_pc_site_moz, left_margin=5Plots.mm, bottom_margin=3Plots.mm)

pl_pc_cov_moz = pl_pc_coverage = plot_sens_adaptor(
    path_pcs_moz, single_vis_pc_coverage!;
    xlabel="ES population coverage (%)", 
    legendtitle="pc",
    legend=:bottomright,
    line5_set...
)

plot(pl_pc_site_airport, pl_pc_cov_airport, 
    pl_pc_site_moz, pl_pc_cov_moz,
    layout= @layout[a b; c d],
    fmt=:png, dpi=300, size=(660, 250*2),
)
# -





