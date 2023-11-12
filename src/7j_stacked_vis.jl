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

path_res1 = "../dt_tmp_res/20231112_005011.ser" # 5000 samples
path_res2 = "../dt_tmp_res/20231112_015655.ser" # 5000 samples
path_res3 = "../dt_tmp_res/20231112_043848.ser" # 5000 samples

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

per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
sens_index = obtain_ES_sensitivity_index(sp_pars.pop, 0.01)
x_per_pop =  per_pop[sens_index]
nothing

path_params, res_all1 = deserialize(path_res1)
df_res1, df_res2, df_res3 = res_all1
df_res = df_res3 # ES population coverage 
# ES population coverage
nothing

include("model_meta_pop.jl")

df_fil, bin_labels = create_lead_time_category(df_res)
nothing

pl = df_to_heatmap(df_res, x_per_pop, :ind_site; add_zero=true)
add_reverse_order_legend!(pl, bin_labels)
plot!(pl, 
    left_margin=5Plots.mm, right_margin=40Plots.mm, 
    xlabel="ES population coverage (%)",
    ylabel="Proportion (%)", 
    xlabelfontsize=14, ylabelfontsize=14,
)
plot!(pl, legend=(1.1, 0.9), dpi=300, fmt=:png)
display(pl)

pl = df_to_heatmap(df_res, sens_index, :ind_site; add_zero=true)
plot!(xlabel="Number of ES covered sites", ylabel="Proportion (%)", xlim=[0,100])

# ## Multiple Figures 

include("model_meta_pop.jl")

path_params1, res_all1 = deserialize(path_res1)
path_params2, res_all2 = deserialize(path_res2)
path_params3, res_all3 = deserialize(path_res3)
println(path_params1)
nothing

paths = fetch_sim_paths(path_params1)
res = deserialize(paths[1])
pc = res.pars.pc

pc=0.25

pl1 = df_to_heatmap(res_all1[3], x_per_pop, :ind_site; 
    add_zero=true, pc=pc)
plot!(pl1, 
    xlabel="ES population coverage (%)",
    ylabel="Proportion (%)",
    title="Population size scenario",
    left_margin=5Plots.mm, 
)
pl2 = df_to_heatmap(res_all2[3], x_per_pop, :ind_site; 
    add_zero=true, pc=pc)
plot!(pl2, 
    xlabel="ES population coverage (%)",
    #title="Airport scenario", 
    tmargin=50Plots.mm
)
airport_cov = per_pop[[11, 7, 62]]
for cov in airport_cov
    annotate!(cov, 100, text("↓", :bottom, 20, :black))
end
pl3 = df_to_heatmap(res_all3[3], x_per_pop, :ind_site; 
    add_zero=true, pc=pc)
plot!(pl3, 
    xlabel="ES population coverage (%)",
    #right_margin=40Plots.mm, 
    ylabel="Proportion (%)",
    title="Mozambique scenario",
)
pl4 = plot(showaxis = false, foreground_color_grid=:white,
    legend=(0.1, 0.9),
)
add_reverse_order_legend!(pl4, bin_labels)
pls = [pl1, pl2, pl3, pl4]
l = @layout [a b; c d]
pl = plot(pls..., 
    layout=l, dpi=300,
    bottom_margin=5Plots.mm,
    xtickfontsize=10, ytickfontsize=10,
)
plot!(pl, size=(1000, 700), fmt=:png)

include("model_meta_pop.jl")

# +
xlim = [0, 100]
xticks = [0, 20, 40, 60, 80, 100]
pl1 = df_to_heatmap(res_all1[3], sens_index, :ind_site; 
    add_zero=true, xticks=xticks)
plot!(pl1, 
    xlabel="Number of ES covered sites",
    ylabel="Proportion (%)",
    title="Population size scenario",
    left_margin=5Plots.mm, 
    xlim=xlim,
    
)
pl2 = df_to_heatmap(res_all2[3], sens_index, :ind_site; 
    add_zero=true, xticks=xticks)
plot!(pl2, 
    xlabel="Number of ES covered sites",
    #title="Airport scenario", 
    tmargin=50Plots.mm,
    xlim=xlim,
)
airport_cov = [11, 7, 62]
for cov in airport_cov
    annotate!(cov, 100, text("↓", :bottom, 20, :black))
end
pl3 = df_to_heatmap(res_all3[3], sens_index, :ind_site; 
    add_zero=true, xticks=xticks)
plot!(pl3, 
    xlabel="Number of ES covered sites",
    #right_margin=40Plots.mm, 
    ylabel="Proportion (%)",
    title="Mozambique scenario",
    xlim=xlim,
)
pl4 = plot(showaxis = false, foreground_color_grid=:white,
    legend=(0.1, 0.9),
)
add_reverse_order_legend!(pl4, bin_labels)
pls = [pl1, pl2, pl3, pl4]
l = @layout [a b; c d]
pl = plot(pls..., 
    layout=l, dpi=300,
    bottom_margin=5Plots.mm,
    xtickfontsize=10, ytickfontsize=10,
)
plot!(pl, size=(1000, 700), fmt=:png)
# -

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

# ## Sampling frequency and detection probabilities

include("model_meta_pop.jl")

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
    
    pl = df_to_heatmap(res_all1[3], x_per_pop, :ind_site)
    ylabel =  (i % 2 == 1) ? "Proportion (%)" : ""
    xlabel =  (i > 2) ? "ES population coverage (%)" : ""
    bottom_margin = (i > 2) ? 5Plots.mm : 5Plots.mm
    plot!(pl, 
        left_margin=5Plots.mm, 
        xlabel=xlabel,
        ylabel=ylabel,
        titlefontsize=12,
        bottom_margin=bottom_margin,
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
    
    pl = df_to_heatmap(res_all1[3], x_per_pop, :ind_site)
    ylabel =  (i % 2 == 1) ? "Proportion (%)" : ""
    xlabel =  (i > 2) ? "ES population coverage (%)" : ""
    plot!(pl, 
        left_margin=5Plots.mm, 
        xlabel=xlabel,
        ylabel=ylabel,
        titlefontsize=12,
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






