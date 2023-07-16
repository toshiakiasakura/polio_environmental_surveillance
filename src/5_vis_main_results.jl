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
using FreqTables
using LinearAlgebra
using DataFrames
using Dates
using Pipe
using Plots
using Rasters
using Serialization

include("util.jl")
include("model_meta_pop.jl")
# -

path = "../dt_tmp/spatial_params_agg230.ser"
sp_pars = deserialize(path)
nothing

# ## Visualise all 

path_res = "../dt_tmp_res/20230716_214255.ser"   # 100 simulations, n_sim=100, R0=14.0, α=0.05 , 84.87
path_params, res_all = deserialize(path_res)
df_res1, df_res2, df_res3 = res_all
nothing

paths = fetch_sim_paths(path_params)
res = deserialize(paths[1])
n_sim = length(paths)
res.pars |> dump
vis_ES_sensitivity_res(res_all, n_sim, sp_pars)







# ## R0 and sensitivity results 

using ColorSchemes
using PyFormattedStrings

per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
sens_index = obtain_ES_sensitivity_index(sp_pars.pop, 0.01)
# TODO: adjust per_pop if not all sites are sampled.
nothing

# +
path_res = "../dt_tmp_res/20230703_184604.ser" # 5000 simulations for R0=1.05
path_params, res_all = deserialize(path_res)
df_res1, df_res2, df_res3_1 = res_all

path_res = "../dt_tmp_res/20230703_192630.ser" # 5000 simulations for r0=1.1
path_params, res_all = deserialize(path_res)
df_res1, df_res2, df_res3_2 = res_all
nothing
# -

paths = fetch_sim_paths(path_params)
n_sim = length(paths)

path_lis_R0 = [
    "../dt_tmp_res/20230703_184604.ser", # 5000 simulations for R0=1.05
    "../dt_tmp_res/20230703_192630.ser", # 5000 simulations for r0=1.1
]
#cmap = palette(:default)



# +
pls = []
for i in 1:6
    legendtitle = i in [1,3,5] ? "AFP surv. or ES" : "50th lead tiem"
    push!(pls, plot(legendtitle=legendtitle, legendtitlefontsize=8))
end

ylabel_prob = "Probability (%)"
ylabel_lead = "Lead time of ES (days)"

for (i, path) in enumerate(path_lis)
    path_params, res_all = deserialize(path)
    _, _, pars = deserialize(path_params * "/1.ser") 
    R0 = round(pars.R0, digits=2)
    label = f"R0 = {R0:.2f}"
    function part_plot_detection_pattern(pl, df_res::DataFrame, col; kargs...)
        tab = nothing
        if col == :per_pop
            tab = detection_pattern_sensitivity(df_res, :ind_site)
            tab[:, :per_pop] = per_pop[sens_index]
        else
            tab = detection_pattern_sensitivity(df_res, col)
        end
        plot!(pl, tab[:, col], tab[:, :Both]/n_sim*100, 
            label=label,
            ylabel=ylabel_prob; 
            kargs...
            )
    end
    function part_plot_diff_sensitivity(pl, df_res::DataFrame, col; kargs...)
        df_diff = nothing
        if col == :per_pop
            df_diff = leadtime_diff_sensitivity(df_res, :ind_site)
            df_diff[:, :per_pop] = per_pop
        else
            df_diff = leadtime_diff_sensitivity(df_res, col)
        end
        plot!(pl, df_diff[:, col], df_diff[:, :q50],
            label=label,
            ylabel=ylabel_lead; 
            kargs...
        )
    end
    
    # ES catchment area
    df_res1, df_res2, df_res3 = res_all
    df_res = df_res3
    part_plot_detection_pattern(pls[1], df_res, :per_pop; 
        xlabel="Population coverage (%)")
    part_plot_diff_sensitivity(pls[2], df_res, :per_pop; 
        xlabel="Population coverage (%)")
    
    # Sampling frequency
    df_res = df_res2
    part_plot_detection_pattern(pls[3], df_res, :n_freq; 
        xlabel="Sampling frequency (days)")
    part_plot_diff_sensitivity(pls[4], df_res, :n_freq; 
        xlabel="Sampling frequency (days)")
    
    # ES sensitivity
    df_res = df_res1
    part_plot_detection_pattern(pls[5], df_res, :pop90; 
        xlabel="ES sensitivity, Np90")
    part_plot_diff_sensitivity(pls[6], df_res, :pop90; 
        xlabel="ES sensitivity, Np90")
end
plot(pls..., layout=(3,2), size=(800, 300*3),
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






