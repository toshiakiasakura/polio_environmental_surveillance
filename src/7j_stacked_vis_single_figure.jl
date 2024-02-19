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
using PyFormattedStrings

include("utils.jl")
include("model_meta_pop.jl")
include("visualise_fig.jl")
# -

# # Main figure 

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
    path_pcs_airport, single_vis_pc_site!; 
    xlabel="Number of ES covered site",
    legendtitle="pc",
    legend = :bottomright,
    line5_set...
)
plot!(pl_pc_site_moz, left_margin=5Plots.mm, bottom_margin=3Plots.mm)

pl_pc_cov_moz = pl_pc_coverage = plot_sens_adaptor(
    path_pcs_airport, single_vis_pc_coverage!;
    xlabel="ES population coverage (%)", 
    legendtitle="pc",
    legend=:bottomright,
    line5_set...
)

plot(pl_pc_site_airport, pl_pc_cov_airport, pl_pc_site_moz, pl_pc_cov_moz,
    layout= @layout[a b; c d],
    fmt=:png, dpi=300, size=(660, 250*2),
)
# -



# # Test

include("visualise_fig.jl")

tab = check_single_percentage(path_res1)
inds = [4, 5, 1, 3, 6, 2]
tab = tab[inds,:]
bin_labels = ["AFP only", "<-60 LT", "-60 ~ -1 LT", "0 ~ 59 LT", "≥60 LT", "ES only"]
x_grouped = @pipe [x for x in tab[:,:prop]] |> reshape(_, 1, 6)
colors = discretise_balance_color(bin_labels)
groupedbar(x_grouped,
    bar_position = :stack,
    bar_width=0.7,
    xticks=(1.0, ""),
    size=(300,600),
    legend=(1.1,0.5),
    right_margin=30Plots.mm,
    left_margin=10Plots.mm,
    labels=reshape(bin_labels, 1, 6),
    color=colors[:, end:-1:begin],
    ylabel="Probability of each pattern (%)",
    ylabelfontsize=12,
    ytickfontsize=12,
    fmt=:png, dpi=200,
    )

include("visualise_fig.jl")





pl_R0 = plot_sens_adaptor(
    path_Rs, single_vis_R0!; 
    legendtitle="R0",
    line5_set...
)
plot(pl_R0, size=(300,300), dpi=300)
#display(pl_R0)

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







# ## Sampling frequency





load(path_res1)

pl_α = plot_α_sens(path_freqs)

plot(pl_R0, pl_α)

include("visualise_fig.jl")

# +

grp_50 = fetch_early_det_50(path_res1)
pl = plot(grp_50[:, :ind_site], grp_50[:, :prop_50], xlim=[0,90])

# -

grp_50 = fetch_early_det_50(path_res1)
pl = plot(grp_50[:, :ind_site], grp_50[:, :prop_0], xlim=[0,90])
grp_50 = fetch_early_det_50(path_res2)
plot!(pl, grp_50[:, :ind_site], grp_50[:, :prop_0], xlim=[0,90])
grp_50 = fetch_early_det_50(path_res3)
plot!(pl, grp_50[:, :ind_site], grp_50[:, :prop_0], xlim=[0,90])





# # Figure check

path_res = "../dt_tmp_res/sens_ES_catchment_20240216_114507.jld2" # 200 simulations, International
path_res2 = "../dt_tmp_res/sens_ES_catchment_20240216_114741.jld2" 

path_trans = load(path_res)["path_trans"]

# Visualise.
paths = fetch_sim_paths(path_trans)
res = load(paths[1])["data"]
pc = res.pars.pc

res |> keys

res.pars

inc_prop = load(path_res)["inc_prop"]
sim_res = load(path_res)["sim_res"]



# # Result check
# single check will be removed after the long run results are obtained.

["R0=14, α=0.05, pc=0.25, n_sim=1000, pattern=population_size, ES_pattern=ES_population_size, path: ../dt_tmp_res/sens_ES_catchment_20240216_190000.jld2", "R0=14, α=0.05, pc=0.25, n_sim=1000, pattern=airport, ES_pattern=ES_population_size, path: ../dt_tmp_res/sens_ES_catchment_20240216_191737.jld2", "R0=14, α=0.05, pc=0.25, n_sim=1000, pattern=mozambique, ES_pattern=ES_population_size, path: ../dt_tmp_res/sens_ES_catchment_20240216_200927.jld2"] |> println

path_res1 = "../dt_tmp_res/sens_ES_catchment_20240216_190000.jld2"
path_res2 = "../dt_tmp_res/sens_ES_catchment_20240216_191737.jld2"
path_res3 = "../dt_tmp_res/sens_ES_catchment_20240216_200927.jld2"

kwds = Dict(:x_var=>"site", :xlim=>[0,90])
single_figure(path_res1; kwds...)
single_figure(path_res2; kwds...)
single_figure(path_res3; kwds...)

# ## R0 sensitivity

["R0=10, α=0.05, pc=0.25, n_sim=1000, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240216_184900.jld2", "R0=12, α=0.05, pc=0.25, n_sim=1000, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240216_191203.jld2", "R0=16, α=0.05, pc=0.25, n_sim=1000, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240216_195206.jld2"] |> println

path_R0_10 = "../dt_tmp_res/sens_ES_catchment_20240216_184900.jld2"
path_R0_14 = "../dt_tmp_res/sens_ES_catchment_20240216_191203.jld2"
path_R0_16 = "../dt_tmp_res/sens_ES_catchment_20240216_195206.jld2"

single_figure(path_R0_10; kwds...)
single_figure(path_R0_14; kwds...)
single_figure(path_R0_16; kwds...)

# ## α sensitivty

["R0=14, α=0.005, pc=0.25, n_sim=1000, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240216_202911.jld2", "R0=14, α=0.01, pc=0.25, n_sim=1000, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240216_205446.jld2", "R0=14, α=0.1, pc=0.25, n_sim=1000, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240216_211527.jld2", "R0=14, α=0.5, pc=0.25, n_sim=1000, pattern=population_size, path: ../dt_tmp_res/sens_ES_catchment_20240216_213551.jld2"] |> println

path_α_0005 = "../dt_tmp_res/sens_ES_catchment_20240216_202911.jld2"
path_α_001 = "../dt_tmp_res/sens_ES_catchment_20240216_205446.jld2"
#path_α_01 = " ../dt_tmp_res/sens_ES_catchment_20240216_211527.jld2"
path_α_05 = "../dt_tmp_res/sens_ES_catchment_20240216_213551.jld2"

single_figure(path_α_0005; kwds...)
single_figure(path_α_001; kwds...)
#single_figure(path_α_01; kwds...)
single_figure(path_α_05; kwds...)

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

path_res1 = "../dt_tmp_res/20231112_005011.ser" # 5000 samples
path_res2 = "../dt_tmp_res/20231114_211615.ser" # Importation is radiation
path_res3 = "../dt_tmp_res/20231112_043848.ser" # 5000 samples

path_res1_moz = "../dt_tmp_res/20231117_111513.ser"
path_res2_moz = "../dt_tmp_res/20231117_125438.ser"
path_res3_moz = "../dt_tmp_res/20231117_162723.ser"

# +
path_R10 = "../dt_tmp_res/20231117_101212.ser"
path_R12 = "../dt_tmp_res/20231117_122254.ser"
path_R16 = "../dt_tmp_res/20231117_152940.ser"

path_α0005 = "../dt_tmp_res/20231117_172755.ser" # 5000 simulations, R0=14.0, α=0.005,
path_α001 = "../dt_tmp_res/20231117_191759.ser" # 5000 simulations, R0=14.0, α=0.01,
path_α01= "../dt_tmp_res/20231117_215227.ser" # 5000 simulations, R0=14.0, α=0.1,
path_α05 = "../dt_tmp_res/20231118_013048.ser" # 5000 simulations, R0=14.0, α=0.5,

path_pc1 = "../dt_tmp_res/20231121_232432.ser"
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
println(paths[1])
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

include("model_meta_pop.jl")
include("visualise_fig.jl")

path_spatial = "../dt_tmp/spatial_params_agg230.ser"

tab = check_single_percentage(path_res1)
inds = [4, 5, 1, 3, 6, 2]
tab = tab[inds,:]
bin_labels = ["AFP only", "<-60 LT", "-60 ~ -1 LT", "0 ~ 59 LT", "≥60 LT", "ES only"]
x_grouped = @pipe [x for x in tab[:,:prop]] |> reshape(_, 1, 6)
colors = discretise_balance_color(bin_labels)
groupedbar(x_grouped,
    bar_position = :stack,
    bar_width=0.7,
    xticks=(1.0, ""),
    size=(300,600),
    legend=(1.1,0.5),
    right_margin=30Plots.mm,
    left_margin=10Plots.mm,
    labels=reshape(bin_labels, 1, 6),
    color=colors[:, end:-1:begin],
    ylabel="Probability of each pattern (%)",
    ylabelfontsize=12,
    ytickfontsize=12,
    fmt=:png, dpi=200,
    )

single_figure(path_spatial, path_res1; x_var="coverage")

single_figure(path_spatial, path_res1; x_var="site",
    xlim=[0,90]
)

# ## Multiple Figures

include("model_meta_pop.jl")
include("visualise_fig.jl")

three_scenario_results(path_spatial, path_res1, path_res2, path_res3; x_var="coverage")

three_scenario_results(path_spatial, path_res1, path_res2, path_res3;
    x_var="site", xlim=[0,100],
)





# ## Mozambique scenario

path_spatial_sorted = "../dt_tmp/spatial_params_agg230_moz_sorted.ser"

three_scenario_results(
    path_spatial_sorted,
    path_res1_moz,  path_res2_moz,  path_res3_moz;
    x_var="coverage", airport_order="mozambique")

include("visualise_fig.jl")

three_scenario_results(
    path_spatial_sorted,
    path_res1_moz, path_res2_moz, path_res3_moz;
    x_var="site", xlim=[0,100], airport_order="mozambique")

# ## Susceptible population

path_params, res_all1 = deserialize(path_res_unvac)
# ES population coverage
per_pop = cumsum(sp_pars_unvac.pop)/sum(sp_pars_unvac.pop)*100
sens_index = obtain_ES_sensitivity_index(sp_pars_unvac.pop, 0.01)
x_unvac = per_pop[sens_index]
nothing

pl1 = df_to_heatmap(res_all1[3], x_unvac, :ind_site; add_zero=true)
plot!(pl1,
    left_margin=5Plots.mm,
    xlabel="ES population coverage (%)",
    ylabel="Proportion (%)",
    title="Susceptible order",
)

# # Sampling frequency

# +
path_freq1 = "../dt_tmp_res/20231125_231613.ser" # 5000 samples
path_freq60 = "../dt_tmp_res/20231126_001339.ser" # 5000 samples

# TODO: fil lthe value here.
# -

single_figure(path_spatial, path_freq1; x_var="site",
    xlim=[0,90]
)

single_figure(path_spatial, path_freq60; x_var="site",
    xlim=[0,90]
)



# ## Sampling frequency and detection probabilities

include("model_meta_pop.jl")

sp_pars = deserialize(path_spatial)
per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
sens_index = obtain_ES_sensitivity_index(sp_pars.pop, 0.01)
nothing

#path_res1 = "../dt_tmp_res/20230719_175720.ser" # 2000 simulations, R0=14.0, α=0.05,
path_params, res_all1 = deserialize(path_res1)
df_res1, df_res2, df_res3 = res_all1
df_res = df_res2
x = df_res[: , :n_freq] |> unique
nothing

# +
x = df_res2[:, :n_freq] |> unique
xticks = [10, 20, 30, 40, 50]
pl1 = df_to_heatmap(df_res2, x, :n_freq; xticks=xticks)
plot!(pl1,
    left_margin=5Plots.mm,
    xlabel="Sampling frequency (day)",
    ylabel="Proportion (%)",
)
x = df_res1[:, :pop90] |> unique

xticks = [1, 3, 10, 30, 100, 300, 1000]
pl2 = df_to_heatmap(df_res1, x, :pop90; xticks=xticks)
plot!(pl2,
    left_margin=-5Plots.mm,
    xticks=(xticks, xticks),
    xlabel="ES detection sensitivity",
    xscale=:log10,
)

plot(pl1, pl2, size=(1200, 500),
    bottom_margin=10Plots.mm, fmt=:png,
)
# -

# ## R0 sensitivity

include("model_meta_pop.jl")

pls = []
for (i, path) in enumerate(path_lis_R0)
    path_params, res_all1 = deserialize(path)
    _, _, pars = deserialize(path_params * "/1.ser")
    R0 = round(pars.R0, digits=2)

    #pl = df_to_heatmap(res_all1[3], x_per_pop, :ind_site)
    pl = df_to_heatmap(res_all1[3], sens_index, :ind_site)
    ylabel =  (i % 2 == 1) ? "Proportion (%)" : ""
    #xlabel =  (i > 2) ? "ES population coverage (%)" : ""
    xlabel =  (i > 2) ? "Number of ES covered sites" : ""
    bottom_margin = (i > 2) ? 5Plots.mm : 5Plots.mm
    plot!(pl,
        left_margin=5Plots.mm,
        xlabel=xlabel,
        ylabel=ylabel,
        titlefontsize=12,
        bottom_margin=bottom_margin,
        xlim=[0,90],
    )
    annotate!(pl, (0.9, 0.9), text("R0=$(R0)", :white, :right), fontsize=12, fontcolor="white")
    push!(pls, pl)
end
plot(pls...,
    size=(1000,600), fmt=:png,
)

# ## α sensitivity

pls = []
for (i, path) in enumerate(path_lis_α)
    path_params, res_all1 = deserialize(path)
    _, _, pars = deserialize(path_params * "/1.ser")
    α = round(pars.α, digits=3)

    #pl = df_to_heatmap(res_all1[3], x_per_pop, :ind_site; add_zero=true,)
    pl = df_to_heatmap(res_all1[3], sens_index, :ind_site; add_zero=true)
    ylabel =  (i % 2 == 1) ? "Proportion (%)" : ""
    #xlabel =  (i > 2) ? "ES population coverage (%)" : ""
    #xlabel =  "ES population coverage (%)"
    xlabel =  "Number of ES covered sites"
    plot!(pl,
        left_margin=5Plots.mm,
        xlabel=xlabel,
        ylabel=ylabel,
        titlefontsize=12,
        xlim=[0,90],
    )
    annotate!(pl, (0.9, 0.9), text(f"α={α:.3f}", :white, :right), fontsize=12, fontcolor="white")
    #hline!(pl, [50], color=:black, linestyle=:dot, alpha=0.7, label=:none)
    #plot!(pl, [100,0], [0,100], lw=3, color=:black, ls=:dashdot, alpha=0.5)
    #hline!(pl, yticks, color=:white, alpha=0.5, ls=:dash, label=:none)
    #vline!(pl, yticks, color=:white, alpha=0.5, ls=:dash, label=:none)
    push!(pls, pl)
end
l = @layout [ a b; c d; e f]
plot(pls...,
    size=(1000,900),
    layout=l, fmt=:png
)

# ## pc sensitivity

single_figure(path_spatial, path_pc1; x_var="coverage",
)

single_figure(path_spatial, path_pc1; x_var="site",
    xlim=[0,90]
)






