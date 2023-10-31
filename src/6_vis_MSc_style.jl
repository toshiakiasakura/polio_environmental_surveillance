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

# ## Baseline results

# +
path_params, res_all1 = deserialize(path_res1)
df, _, _ = res_all1
cond = df[: ,:pop90] .== 10
df = df[cond, :]
#first(df, 5)

#outcome_num_prop(df) |> display
#detect_pattern(df) |> countmap |> display

dfM = @pipe filter(x -> isnan(x["t_ES"]) == false, df)  
m_ES  = @pipe mean(dfM[:, "t_ES"]) |> round(_, digits=2)
dfM = @pipe filter(x -> isnan(x["t_AFP"]) == false, df)  
m_AFP = @pipe mean(dfM[:, "t_AFP"])  |> round(_, digits=2)
println("Mean ES: $m_ES days, Mean AFP: $m_AFP days")

dif = leadtime_diff(df)
leadtime_diff_statistics(dif)

df_base = copy(df)
nothing
# -


# ## International airports  

# +
path_params, res_all1 = deserialize(path_res2)
df, _, _ = res_all1
cond = df[: ,:pop90] .== 10
df = df[cond, :]
#first(df, 5)

#outcome_num_prop(df) |> display
#detect_pattern(df) |> countmap |> display

dfM = @pipe filter(x -> isnan(x["t_ES"]) == false, df)  
m_ES  = @pipe mean(dfM[:, "t_ES"]) |> round(_, digits=2)
dfM = @pipe filter(x -> isnan(x["t_AFP"]) == false, df)  
m_AFP = @pipe mean(dfM[:, "t_AFP"])  |> round(_, digits=2)
println("Mean ES: $m_ES days, Mean AFP: $m_AFP days")

dif = leadtime_diff(df)
leadtime_diff_statistics(dif)
df_inter = copy(df)
nothing
# -

# ## Figure 1.  

# +
days = res.pars.days
ts_ES = df_base[:, "t_ES"]
ts_AFP = df_base[:, "t_AFP"]
t_extinct = df_base[:, "t_extinct"]


# Baseline scenario 
pl1 = plot(
    xlabel="Day", ylabel="Cumulative probability \nof the first detection",
    legend=(0.6, 0.5),
)
annotate!((0.07, 0.95), "(A)")
cum_ES = cumulative_counts(ts_ES, days; prop=true)
cum_AFP = cumulative_counts(ts_AFP, days; prop=true)
cum_ES_cond = conditional_cumulative_prob(ts_ES, t_extinct, days)
cum_AFP_cond = conditional_cumulative_prob(ts_AFP, t_extinct, days)

plot!(pl1, 1:days, cum_ES, label="Prob. via ES", color=1)
plot!(pl1, 1:days, cum_ES_cond, label="Conditional Prob. via ES", color=1, linestyle=:dashdot)
plot!(pl1, 1:days, cum_AFP, label="Prob. via AFP surv.", color=2)
plot!(pl1, 1:days, cum_AFP_cond, label="Conditional Prob. via AFP surv.", color=2, linestyle=:dashdot)

dif = leadtime_diff(df_base)
x = [1 for i in 1:length(dif)]
pl2 = violin(x, dif, xticks=:none, ylabel="Lead time of ES (day)", legend=:none)
boxplot!(pl2, x, dif, fillalpha=0.75)
annotate!((0.15, 0.95), "(B)")

# International airports. 
days = res.pars.days
ts_ES = df_inter[:, "t_ES"]
ts_AFP = df_inter[:, "t_AFP"]
t_extinct = df_inter[:, "t_extinct"]


pl3 = plot(
    xlabel="Day", ylabel="Cumulative probability \nof the first detection",
    legend=(0.6, 0.5),
)
annotate!((0.07, 0.95), "(C)")
cum_ES = cumulative_counts(ts_ES, days; prop=true)
cum_AFP = cumulative_counts(ts_AFP, days; prop=true)
cum_ES_cond = conditional_cumulative_prob(ts_ES, t_extinct, days)
cum_AFP_cond = conditional_cumulative_prob(ts_AFP, t_extinct, days)

plot!(pl3, 1:days, cum_ES, label="Prob. via ES", color=1)
plot!(pl3, 1:days, cum_ES_cond, label="Conditional Prob. via ES", color=1, linestyle=:dashdot)
plot!(pl3, 1:days, cum_AFP, label="Prob. via AFP surv.", color=2)
plot!(pl3, 1:days, cum_AFP_cond, label="Conditional Prob. via AFP surv.", color=2, linestyle=:dashdot)

dif = leadtime_diff(df_inter)
x = [1 for i in 1:length(dif)]
pl4 = violin(x, dif, xticks=:none, ylabel="Lead time of ES (day)", legend=:none)
boxplot!(pl4, x, dif, fillalpha=0.75)
annotate!((0.15, 0.95), "(D)")
l = @layout [a{0.75w} b]
pl = plot(pl1, pl2, 
    fmt=:png, dpi=300, layout=l,
    size=(800,500), left_margin=5Plots.mm,
)


l = @layout [
    a{0.75w} b; 
    c{0.75w} d
]
pl = plot(pl1, pl2, pl3, pl4, 
    fmt=:png, dpi=300, layout=l,
    size=(900,700), left_margin=5Plots.mm,
)
display(pl)
savefig(pl, "../res/fig_baseline_cum_lead.png")
# -

# ## Visualise main figures

# +
#path_res1 = "../dt_tmp_res/20230719_175720.ser" # 2000 simulations, R0=14.0, α=0.05, 
path_params, res_all1 = deserialize(path_res1)

#path_res2 = "../dt_tmp_res/20230720_090053.ser" # 2000 simulations, R0=14.0, α=0.05, International airpot.
path_params, res_all2 = deserialize(path_res2)
println(path_params)
nothing
# -

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

# ### Check correlations 

per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
sens_index = obtain_ES_sensitivity_index(sp_pars.pop, 0.01)
df_res1, df_res2, df_res3 = res_all1
df_res = df_res3
df_diff = leadtime_diff_sensitivity(df_res, :ind_site)
df_diff[:, :per_pop] = per_pop[sens_index]
first(df_diff, 5)

pop = sp_pars.pop
unvac = sp_pars.unvac
per_pop = cumsum(pop)/sum(pop)*100
sens_index = obtain_ES_sensitivity_index(pop, 0.01)
w_pop_cov = mean((pop .- unvac)./pop, weights(pop))
println(w_pop_cov)
@pipe [12, 13, 14, 15, 16] .* (1-w_pop_cov) .|> round(_, digits=2) |> println 
nothing

# +
ind = abs.(per_pop.- 10) |> argmin
per = per_pop[ind]
println("Precent of population coverage :", per)
df_res_tmp = filter(x -> x.ind_site == ind, df_res) 
df_res_tmp = filter(x -> x.diff == x.diff, df_res_tmp)

xlabel = "Cumulative number of infections\nat 1st ES detection time"
ylabel = "Lead time of ES (day)"
xticks = [1, 10, 100, 1000, 10_000, 100_000]
pl1 = scatter(df_res_tmp[:, :R_inf_ES].+1, df_res_tmp[:, :diff],
    xlabel = xlabel, ylabel = ylabel, 
    xaxis=:log10, xticks=(xticks), 
    markersize=2.5, label=:none,
    #title="Population coverage $(per)",
)
annotate!(pl1, (0.1, 0.1), "(A)")

ind = abs.(per_pop.- 50) |> argmin
per = per_pop[ind]
println("Precent of population coverage :", per)
df_res_tmp = filter(x -> x.ind_site == ind, df_res) 
df_res_tmp = filter(x -> x.diff == x.diff, df_res_tmp)
xticks = [1, 10, 100, 1000]
pl2 = scatter(df_res_tmp[:, :R_inf_ES].+1, df_res_tmp[:, :diff],
    xlabel = xlabel, ylabel = ylabel, 
    xaxis=:log10,
    xticks=(xticks),
    markersize=2.5, label=:none,
    #title="Population coverage $(per)",
)
annotate!(pl2, (0.1, 0.1), "(B)")

pl = plot(pl1, pl2, size=(800,400), 
    left_margin=5Plots.mm, 
    bottom_margin=5Plots.mm, 
    fmt=:png,
    dpi=300,
)
display(pl)
savefig(pl, "../res/fig_cum_num_at_ES_det.png")
# -
# ## Susceptible population

path_params, res_all1 = deserialize(path_res_unvac)
nothing

# +
df_res1, df_res2, df_res3 = res_all1
xlabel = "ES population coverage (%)"

# ES population coverage
per_pop = cumsum(sp_pars_unvac.pop)/sum(sp_pars_unvac.pop)*100
sens_index = obtain_ES_sensitivity_index(sp_pars_unvac.pop, 0.01)
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
    xticks=[0,20,40,60,80,100],
    ylim=[-300, 300], legend=(0.9,0.3),
    )
    df_det = early_detect_prob_sensitivity(df_res, :ind_site) 

df_det = early_detect_prob_sensitivity(df_res, :ind_site) 
df_det[:, :per_pop] = per_pop[sens_index]
pl3 = vis_early_detect_prob(df_det, :per_pop,
    xlabel=xlabel, title="", legend=(0.7, 0.2)
)
plot!(pl3, [0, 100.0], [0,100.0], color=:black, label=:none, alpha=0.5, linestyle=:dash)
annotate!(pl1, (0.08, 0.95), "(A)")
annotate!(pl2, (0.08, 0.95), "(B)")
annotate!(pl3, (0.08, 0.95), "(C)")

pls = [pl1, pl2, pl3]
pl = plot(pls...,
    size=(400*3, 300*1),  # 400 * 2
    layout=(1,3),
    dpi=300,
    left_margin = 7Plots.mm, bottom_margin = 7Plots.mm,
)
savefig("../res/fig_es_population_cov_unvac.png")
display(pl)
#plot!(pl3, [0, 100.0], [0,100.0], color=:black, label=:none, linestyle=:dot)
# -

# ## ES catchment sensitivity
# R0 and sensitivity results 

pop = sp_pars.pop
unvac = sp_pars.unvac
per_pop = cumsum(pop)/sum(pop)*100
sens_index = obtain_ES_sensitivity_index(pop, 0.01)
w_pop_cov = mean((pop .- unvac)./pop, weights(pop))
nothing

paths = fetch_sim_paths(path_params)
n_sim = length(paths)

include("model_meta_pop.jl")

xlabel = "ES population coverage (%)"
ylabel_prob = "Detection pattern (%)"
ylabel_lead = "Lead time of ES (day)"
ylabel_early = "Probability of \nearly detection by ES"

# +
pls = []
for i in 1:6
    #legendtitle = i in [1,3,5] ? "AFP surv. or ES" : "50th lead tiem"
    push!(pls, plot(legendtitlefontsize=8))
end

colors = [1,2,3,4,5, 6]
for (i, path) in enumerate(path_lis_R0)
    path_params, res_all = deserialize(path)
    _, _, pars = deserialize(path_params * "/1.ser") 
    R0 = round(pars.R0, digits=2)
    Re0 = pars.R0 * (1-w_pop_cov)
    # ES catchment area
    label = f"{R0:.0f}" #, Re0 = {Re0:.2f}"
    df_res1, df_res2, df_res3 = res_all
    df_res = df_res3
    part_plot_detection_pattern(pls[1], df_res, :per_pop; 
        xlabel=xlabel,
        legend=(0.9,0.9),
        label=label,
        color=colors[i], legendtitle="R0",
    )
    part_plot_diff_sensitivity(pls[2], df_res, :per_pop; 
        xlabel=xlabel,
        legend=(0.9,0.4),
        label=label,
        ylim=[-200,200], legendtitle="R0",
    )
    part_plot_early_det(pls[3], df_res, :per_pop; 
        xlabel=xlabel,
        legend=(0.9,0.4),
        label=label, legendtitle="R0",
    )
end
plot!(pls[3], [0, 100.0], [0,100.0], color=:black, label=:none, alpha=0.5, linestyle=:dash)

for (i, path) in enumerate(path_lis_α)
    path_params, res_all = deserialize(path)
    _, _, pars = deserialize(path_params * "/1.ser") 
    α = round(pars.α, digits=3)
    label = f"{α:.3f}"
    
    # ES catchment area
    df_res1, df_res2, df_res3 = res_all
    df_res = df_res3
    part_plot_detection_pattern(pls[4], df_res, :per_pop; 
        xlabel=xlabel,
        label=label,
        color=colors[i], legendtitle="α",
    )
    part_plot_diff_sensitivity(pls[5], df_res, :per_pop; 
        xlabel=xlabel,
        label=label,
        legend=(0.9,0.4),
        ylim=[-200,200], legendtitle="α",
    )
    part_plot_early_det(pls[6], df_res, :per_pop; 
        xlabel=xlabel,
        legend=(0.9,0.4),
        label=label, legendtitle="α",
    )
end
plot!(pls[6], [0, 100.0], [0,100.0], color=:black, label=:none, alpha=0.5, linestyle=:dash)
annotate!(pls[1], (0.08, 0.95), "(A)")
annotate!(pls[2], (0.08, 0.95), "(B)")
annotate!(pls[3], (0.08, 0.95), "(C)")
annotate!(pls[4], (0.08, 0.95), "(D)")
annotate!(pls[5], (0.08, 0.95), "(E)")
annotate!(pls[6], (0.08, 0.95), "(F)")
pl = plot(pls..., layout=(2,3), size=(1200, 300*2),
    left_margin=5Plots.mm, bottom_margin=5Plots.mm,
    fmt=:png,
)
display(pl)
savefig(pl, "../res/fig_es_population_cov_sens.png")
# -

# ## Sampling frequency and detection probability

function vis_sampling_ES_sensitivity(res_all)
    df_res1, df_res2, df_res3 = res_all
    xlabel = "Sampling frequency (day)"
    
    # Sampling frequency
    df_res = df_res2
    tab = detection_pattern_sensitivity(df_res, :n_freq)
    pl1 = vis_detection_pattern(tab, n_sim, :n_freq; 
        xlabel=xlabel, title="",
    )
    
    df_diff = leadtime_diff_sensitivity(df_res, :n_freq)
    pl2 = vis_leadtime_diff_sensitivity(df_diff, :n_freq;  
        xlabel=xlabel, title="",
        #ylim=[nothing, 300], 
        legend=(0.9,0.3),
        )

    df_det = early_detect_prob_sensitivity(df_res, :n_freq) 
    pl3 = vis_early_detect_prob(df_det, :n_freq,
        xlabel=xlabel, title="",
    )
    
    # ES sensitivity
    xticks = (1, 3, 10, 30, 100, 300, 1000)
    
    xlabel = "ES detection sensitivity, Np90"
    df_res = df_res1
    tab = detection_pattern_sensitivity(df_res, :pop90)
    pl4 = vis_detection_pattern(tab, n_sim, :pop90; 
        xlabel=xlabel, title="",
        xscale=:log10, xticks=(xticks, xticks),
        legend=(0.6, 0.5),
    )
    
    df_diff = leadtime_diff_sensitivity(df_res, :pop90)
    pl5 = vis_leadtime_diff_sensitivity(df_diff, :pop90;  
        xlabel=xlabel, title="",
        xscale=:log10, xticks=(xticks, xticks),
        #ylim=[nothing, 350], 
        legend=(0.85,0.9),
        )

    df_det = early_detect_prob_sensitivity(df_res, :pop90) 
    pl6 = vis_early_detect_prob(df_det, :pop90,
        xscale=:log10, xticks=(xticks, xticks),
        xlabel=xlabel, title="",
    )
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
    savefig("../res/fig_sampling_ES_sensitivity.png")
    display(pl)
end

include("model_meta_pop.jl")
vis_sampling_ES_sensitivity(res_all1)

# ## Sampling frequency sensitivity

# +
pls = []
for i in 1:6
    #legendtitle = i in [1,3,5] ? "AFP surv. or ES" : "50th lead tiem"
    push!(pls, plot(legendtitlefontsize=8))
end

xlabel = "Sampling frequency (day)"
for (i, path) in enumerate(path_lis_R0)
    path_params, res_all = deserialize(path)
    _, _, pars = deserialize(path_params * "/1.ser") 
    R0 = round(pars.R0, digits=2)
    Re0 = pars.R0 * (1-w_pop_cov)
    # ES catchment area
    label = f"{R0:.0f}" #, Re0 = {Re0:.2f}"
    df_res1, df_res2, df_res3 = res_all
    df_res = df_res2
    part_plot_detection_pattern(pls[1], df_res, :n_freq; 
        xlabel=xlabel,
        legend=(0.9,0.95),
        label=label,
        color=colors[i], legendtitle="R0",
        xlim=[0,60],
    )
    part_plot_diff_sensitivity(pls[2], df_res, :n_freq; 
        xlabel=xlabel,
        legend=(0.9,0.6),
        label=label,
        #ylim=[-100,100], 
        legendtitle="R0",
        xlim=[0,60],
    )
    part_plot_early_det(pls[3], df_res, :n_freq; 
        xlabel=xlabel,
        legend=(0.9,0.9),
        label=label, legendtitle="R0",
        xlim=[0,60],
    )
end

for (i, path) in enumerate(path_lis_α)
    path_params, res_all = deserialize(path)
    _, _, pars = deserialize(path_params * "/1.ser") 
    α = round(pars.α, digits=3)
    label = f"{α:.3f}"
    
    # ES catchment area
    df_res1, df_res2, df_res3 = res_all
    df_res = df_res2
    part_plot_detection_pattern(pls[4], df_res, :n_freq; 
        legend=(0.8,0.95),
        xlabel=xlabel,
        label=label,
        color=colors[i], legendtitle="α",
        xlim=[0,60],
    )
    part_plot_diff_sensitivity(pls[5], df_res, :n_freq; 
        xlabel=xlabel,
        label=label,
        #ylim=[-100,100], 
        legendtitle="α",
        xlim=[0,60],
        legend=(0.8,0.4),
    )
    part_plot_early_det(pls[6], df_res, :n_freq; 
        xlabel=xlabel,
        label=label, legendtitle="α",
        xlim=[0,60],
        legend=(0.8,0.9),
    )
end
annotate!(pls[1], (0.08, 0.95), "(A)")
annotate!(pls[2], (0.08, 0.95), "(B)")
annotate!(pls[3], (0.08, 0.95), "(C)")
annotate!(pls[4], (0.08, 0.95), "(D)")
annotate!(pls[5], (0.08, 0.95), "(E)")
annotate!(pls[6], (0.08, 0.95), "(F)")
pl = plot(pls..., layout=(2,3), size=(1200, 300*2),
    left_margin=5Plots.mm, bottom_margin=5Plots.mm,
    fmt=:png,
)
display(pl)
savefig(pl, "../res/fig_sampling_sens.png")
# -

# ## ES sensitivity variability 

# +
pls = []
for i in 1:6
    #legendtitle = i in [1,3,5] ? "AFP surv. or ES" : "50th lead tiem"
    push!(pls, plot(legendtitlefontsize=8))
end

xticks = (1, 3, 10, 30, 100, 300, 1000)
xlabel = "ES detection sensitivity, Np90"
for (i, path) in enumerate(path_lis_R0)
    path_params, res_all = deserialize(path)
    _, _, pars = deserialize(path_params * "/1.ser") 
    R0 = round(pars.R0, digits=2)
    Re0 = pars.R0 * (1-w_pop_cov)
    # ES catchment area
    label = f"{R0:.0f}" #, Re0 = {Re0:.2f}"
    df_res1, df_res2, df_res3 = res_all
    df_res = df_res1
    part_plot_detection_pattern(pls[1], df_res, :pop90; 
        xlabel=xlabel,
        legend=(0.9,0.5),
        label=label,
        color=colors[i], legendtitle="R0",
        xscale=:log10, xticks=(xticks, xticks),
    )
    part_plot_diff_sensitivity(pls[2], df_res, :pop90; 
        xlabel=xlabel,
        legend=(0.9,0.6),
        label=label,
        #ylim=[-100,100], 
        legendtitle="R0",
        xscale=:log10, xticks=(xticks, xticks),
    )
    part_plot_early_det(pls[3], df_res, :pop90; 
        xlabel=xlabel,
        legend=(0.8,0.9),
        label=label, legendtitle="R0",
        xscale=:log10, xticks=(xticks, xticks),
    )
end

for (i, path) in enumerate(path_lis_α)
    path_params, res_all = deserialize(path)
    _, _, pars = deserialize(path_params * "/1.ser") 
    α = round(pars.α, digits=3)
    label = f"{α:.3f}"
    
    # ES catchment area
    df_res1, df_res2, df_res3 = res_all
    df_res = df_res1
    part_plot_detection_pattern(pls[4], df_res, :pop90; 
        xlabel=xlabel,
        label=label,
        color=colors[i], legendtitle="α",
        xscale=:log10, xticks=(xticks, xticks),
        legend=(0.8,0.5),
    )
    part_plot_diff_sensitivity(pls[5], df_res, :pop90; 
        xlabel=xlabel,
        label=label,
        #ylim=[-100,100], 
        legendtitle="α",
        xscale=:log10, xticks=(xticks, xticks),
    )
    part_plot_early_det(pls[6], df_res, :pop90; 
        xlabel=xlabel,
        label=label, legendtitle="α",
        legend=(0.8,0.9),
        xscale=:log10, xticks=(xticks, xticks),
    )
end
annotate!(pls[1], (0.08, 0.95), "(A)")
annotate!(pls[2], (0.08, 0.95), "(B)")
annotate!(pls[3], (0.08, 0.95), "(C)")
annotate!(pls[4], (0.08, 0.95), "(D)")
annotate!(pls[5], (0.08, 0.95), "(E)")
annotate!(pls[6], (0.08, 0.95), "(F)")
pl = plot(pls..., layout=(2,3), size=(1200, 300*2),
    left_margin=5Plots.mm, bottom_margin=5Plots.mm,
    fmt=:png,
)
display(pl)
savefig(pl, "../res/fig_es_sensitivity_sens.png")
# -

