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
using CategoricalArrays
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
using StatsPlots

include("util.jl")
include("model_meta_pop.jl")
include("visualise_fig.jl")
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






