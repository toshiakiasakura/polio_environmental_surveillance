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

path = "../dt_tmp/spatial_params_agg230.ser"
sp_pars = deserialize(path)
nothing

# ## Visualise all 

# +
path_res = "../dt_tmp_res/20230716_214255.ser"   # 100 simulations, n_sim=100, R0=14.0, α=0.05 , 84.87
path_res = "../dt_tmp_res/20230716_221956.ser" # 1000 simulations, R0=14.0, α=0.1,

path_res = "../dt_tmp_res/20230716_222655.ser"  # 1000 simulations, R0=13.0, α=0.05, 
path_res = "../dt_tmp_res/20230716_220459.ser"  # 1000 simulations, R0=14.0, α=0.05

# +

path_params, res_all = deserialize(path_res)
println(path_params)
df_res1, df_res2, df_res3 = res_all
nothing
# -

include("model_meta_pop.jl")
paths = fetch_sim_paths(path_params)
res = deserialize(paths[1])
n_sim = length(paths)
println("$n_sim simulations, R0=$(res.pars.R0), α=$(res.pars.α), ")
res.pars |> dump
vis_ES_sensitivity_res(res_all, n_sim, sp_pars)

# ### Check correlations 

per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
sens_index = obtain_ES_sensitivity_index(sp_pars.pop, 0.01)
df_res = df_res3
df_diff = leadtime_diff_sensitivity(df_res, :ind_site)
df_diff[:, :per_pop] = per_pop[sens_index]
first(df_diff, 5)

println("Precent of population coverage :", per_pop[5])
df_res_tmp = filter(x -> x.ind_site == 5, df_res) 
df_res_tmp = filter(x -> x.diff == x.diff, df_res_tmp)
pl1 = scatter(df_res_tmp[:, :R_final_num], df_res_tmp[:, :diff],
    xlabel = "Final AFP size",
    ylabel = "Laed time", 
)
pl2 = scatter(df_res_tmp[:, :R_final_site], df_res_tmp[:, :diff],
    xlabel = "Final site",
    ylabel = "Laed time", 
)
plot(pl1, pl2, size=(800,400), left_margin=5Plots.mm, bottom_margin=5Plots.mm, fmt=:png)

println("Precent of population coverage :", per_pop[90] )
df_res_tmp = filter(x -> x.ind_site == 90, df_res) 
df_res_tmp = filter(x -> x.diff == x.diff, df_res_tmp)
pl1 = scatter(df_res_tmp[:, :R_final_num], df_res_tmp[:, :diff],
    xlabel = "Final AFP size",
    ylabel = "Laed time", 
)
pl2 = scatter(df_res_tmp[:, :R_final_site], df_res_tmp[:, :diff],
    xlabel = "Final site",
    ylabel = "Laed time", 
)
plot(pl1, pl2, size=(800,400), left_margin=5Plots.mm, bottom_margin=5Plots.mm, fmt=:png)








# ## R0 and sensitivity results 

pop = sp_pars.pop
unvac = sp_pars.unvac
per_pop = cumsum(pop)/sum(pop)*100
sens_index = obtain_ES_sensitivity_index(pop, 0.01)
w_pop_cov = mean((pop .- unvac)./pop, weights(pop))
# TODO: adjust per_pop if not all sites are sampled.
nothing

paths = fetch_sim_paths(path_params)
n_sim = length(paths)

path_lis_R0 = [
    "../dt_tmp_res/20230716_220459.ser"  # 1000 simulations, R0=14.0, α=0.05
    "../dt_tmp_res/20230716_222655.ser"  # 1000 simulations, R0=13.0, α=0.05, 
]
#cmap = palette(:default)

# +
function part_plot_detection_pattern(pl, df_res::DataFrame, col; label="", kargs...)
    tab = nothing
    if col == :per_pop
        tab = detection_pattern_sensitivity(df_res, :ind_site)
        tab[:, :per_pop] = per_pop[sens_index]
    else
        tab = detection_pattern_sensitivity(df_res, col)
    end
    y = tab[:, :Both]/n_sim*100
    y_max = maximum(y)
    plot!(pl, tab[:, col], tab[:, :Both]/n_sim*100, 
        label=label,
        ylim=[0, y_max],
        ylabel=ylabel_prob; 
        kargs...
        )
end

function part_plot_diff_sensitivity(pl, df_res::DataFrame, col; label="",kargs...)
    df_diff = nothing
    if col == :per_pop
        df_diff = leadtime_diff_sensitivity(df_res, :ind_site)
        df_diff[:, :per_pop] = per_pop[sens_index]
    else
        df_diff = leadtime_diff_sensitivity(df_res, col)
    end
    plot!(pl, df_diff[:, col], df_diff[:, :q50],
        label=label,
        ylabel=ylabel_lead; 
        kargs...
    )
end

function part_plot_early_det(pl, df_res::DataFrame, col; label="", kargs...)
    df_det = nothing
    if col == :per_pop
        df_det = early_detect_prob_sensitivity(df_res, :ind_site)
        df_det[:, :per_pop] = per_pop[sens_index]
    else
        df_det= early_detect_prob_sensitivity(df_res, col)
    end
    y_max = maximum(df_det[:, :lead_30])
    plot!(pl, df_det[:, col], df_det[:, :lead_30],
        label=label,
        ylim=[0,y_max],
        ylabel=ylabel_early; 
        kargs...
    )
end

# +
pls = []
for i in 1:9
    #legendtitle = i in [1,3,5] ? "AFP surv. or ES" : "50th lead tiem"
    push!(pls, plot(legendtitlefontsize=8))
end

ylabel_prob = "Probability (%)"
ylabel_lead = "Lead time of ES (days)"
ylabel_early = "Probability of \nearly detection by ES"

for (i, path) in enumerate(path_lis_R0)
    path_params, res_all = deserialize(path)
    _, _, pars = deserialize(path_params * "/1.ser") 
    R0 = round(pars.R0, digits=2)
    Re0 = pars.R0 * (1-w_pop_cov)
    label = f"R0 = {R0:.2f}, Re0 = {Re0:.2f}"
    
    # ES catchment area
    df_res1, df_res2, df_res3 = res_all
    df_res = df_res3
    part_plot_detection_pattern(pls[1], df_res, :per_pop; 
        xlabel="Population coverage (%)",
        legend=(0.5,0.2),
        label=label,
        title="Detection via AFP surv. or ES",
    )
    part_plot_diff_sensitivity(pls[2], df_res, :per_pop; 
        xlabel="Population coverage (%)",
        legend=(0.5,0.2),
        label=label,
        title="Median of lead time",
    )
    part_plot_early_det(pls[3], df_res, :per_pop; 
        xlabel="Population coverage (%)",
        legend=(0.5,0.2),
        label=label,
        title="Early detection by ES \nwhen LT > 30 days or ES only",
    )
    
    # Sampling frequency
    df_res = df_res2
    part_plot_detection_pattern(pls[4], df_res, :n_freq; 
        xlabel="Sampling frequency (days)",
        label=label,
    )
    part_plot_diff_sensitivity(pls[5], df_res, :n_freq; 
        xlabel="Sampling frequency (days)",
        label=label,
    )
    part_plot_early_det(pls[6], df_res, :n_freq; 
        xlabel="Sampling frequency (days)",
        label=label,
    )
    
    # ES sensitivity
    xticks = (1, 3, 10, 30, 100, 300, 1000)
    df_res = df_res1
    part_plot_detection_pattern(pls[7], df_res, :pop90; 
        xlabel="ES sensitivity, Np90",
        xscale=:log10, xticks=(xticks, xticks),
        label=label,
    )
    part_plot_diff_sensitivity(pls[8], df_res, :pop90; 
        xlabel="ES sensitivity, Np90",
        xscale=:log10, xticks=(xticks, xticks),
        label=label,
    )
    part_plot_early_det(pls[9], df_res, :pop90; 
        xlabel="ES sensitivity, Np90",
        xscale=:log10, xticks=(xticks, xticks),
        label=label,
    )
end
plot(pls..., layout=(3,3), size=(1200, 300*3),
    left_margin=5Plots.mm, bottom_margin=5Plots.mm,
    fmt=:png,
)
# -
# # α sensitivity

path_lis_α = [
    "../dt_tmp_res/20230716_220459.ser"  # 1000 simulations, R0=14.0, α=0.05
    "../dt_tmp_res/20230716_221956.ser" # 1000 simulations, R0=14.0, α=0.1,
]
#cmap = palette(:default)



# +
pls = []
for i in 1:9
    #legendtitle = i in [1,3,5] ? "AFP surv. or ES" : "50th lead tiem"
    push!(pls, plot(legend_title="α", legendtitlefontsize=8))
end

ylabel_prob = "Probability (%)"
ylabel_lead = "Lead time of ES (days)"

for (i, path) in enumerate(path_lis_α)
    path_params, res_all = deserialize(path)
    _, _, pars = deserialize(path_params * "/1.ser") 
    α = round(pars.α, digits=2)
    label = f"{α:.2f}"
    
    # ES catchment area
    df_res1, df_res2, df_res3 = res_all
    df_res = df_res3
    part_plot_detection_pattern(pls[1], df_res, :per_pop; 
        xlabel="Population coverage (%)",
        label=label,
        title="Detection via AFP surv. or ES",
    )
    part_plot_diff_sensitivity(pls[2], df_res, :per_pop; 
        xlabel="Population coverage (%)",
        label=label,
        title="Median of lead time",
    )
    part_plot_early_det(pls[3], df_res, :per_pop; 
        xlabel="Population coverage (%)",
        label=label,
        title="Early detection by ES \nwhen LT > 30 days or ES only",
    )
    
    # Sampling frequency
    df_res = df_res2
    part_plot_detection_pattern(pls[4], df_res, :n_freq; 
        xlabel="Sampling frequency (days)",
        label=label,
    )
    part_plot_diff_sensitivity(pls[5], df_res, :n_freq; 
        xlabel="Sampling frequency (days)",
        label=label,
    )
    part_plot_early_det(pls[6], df_res, :n_freq; 
        xlabel="Sampling frequency (days)",
        label=label,
    )
    
    # ES sensitivity
    xticks = (1, 3, 10, 30, 100, 300, 1000)
    df_res = df_res1
    part_plot_detection_pattern(pls[7], df_res, :pop90; 
        xlabel="ES sensitivity, Np90",
        label=label,
        xscale=:log10, xticks=(xticks, xticks),
    )
    part_plot_diff_sensitivity(pls[8], df_res, :pop90; 
        xlabel="ES sensitivity, Np90",
        label=label,
        xscale=:log10, xticks=(xticks, xticks),
    )
    part_plot_early_det(pls[9], df_res, :pop90; 
        xlabel="ES sensitivity, Np90",
        label=label,
        xscale=:log10, xticks=(xticks, xticks),
    )
end
plot(pls..., layout=(3,3), size=(1200, 300*3),
    left_margin=5Plots.mm, bottom_margin=5Plots.mm,
    fmt=:png,
)
# -






# # Sensitivity analysis code for each parameter. 

# ## Sensitivity hazard 

df_res = sensitivity_ana_ES(path, par_AFP, par_ES, sensitivity_hazard)
nothing

tab = detection_pattern_sensitivity(df_res, :pop90)
xticks = (1, 3, 5, 10, 30, 50, 100, 300, 500, 1000)
xlabel = "Pop 90 for coefficient for hazard, g"
vis_detection_pattern(tab, n_sim, :pop90; xscale=:log10, xticks=(xticks, xticks), xlabel=xlabel)
df_diff = leadtime_diff_sensitivity(df_res, :pop90)
vis_leadtime_diff_sensitivity(df_diff, :pop90; xscale=:log10, xticks=(xticks, xticks))

# ## Sensitivity frequency

df_res = sensitivity_ana_ES(path, par_AFP, par_ES, sensitivity_frequency_sampling)
nothing

# +

tab = detection_pattern_sensitivity(df_res, :n_freq)
vis_detection_pattern(tab, n_sim, :n_freq; xlabel="Frequency of sampling (days)")
df_diff = leadtime_diff_sensitivity(df_res, :n_freq)
vis_leadtime_diff_sensitivity(df_diff, :n_freq; xlabel="Frequency of sampling (days)")
# -

# ### Sensitivity ES catchment_area 

per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
nothing

df_res = sensitivity_ana_ES(path, par_AFP, par_ES, sensitivity_ES_catchment_area)
nothing

tab = detection_pattern_sensitivity(df_res, :ind_site)
tab[:, :per_pop] = per_pop
vis_detection_pattern(tab, n_sim, :ind_site)
vis_detection_pattern(tab, n_sim, :per_pop; xlabel="Population coverage (%)")

df_diff = leadtime_diff_sensitivity(df_res, :ind_site)
vis_leadtime_diff_sensitivity(df_diff, :ind_site)
df_diff[:, :per_pop] = per_pop
vis_leadtime_diff_sensitivity(df_diff, :per_pop;  xlabel="Population coverage (%)")






