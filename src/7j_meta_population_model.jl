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

# # Transmission model 

sp_pars = deserialize("../dt_tmp_fix/pop_pi_agg1000_10unvac.ser")
keys(sp_pars)

n_site = length(sp_pars.pop)
nothing

pars = SEIRMetaModelParams(
    R0=1.10,
    N0=sp_pars.pop,
    π_mat=sp_pars.π_mat,
    n_site=n_site,
    α=0.050,
    days=365*3,
    )
pars |> dump

# ## Run simulations 

include("model_meta_pop.jl")

Random.seed!(47)
pl = plot()
htmaps = []
for _ in 1:6
    rec, outcome, pars = run_sim(pars; rec_flag=true)
    #println(outcome.R_final_num, ", ", outcome.R_final_site)
    ht = heatmap_meta_pop(rec.Ia)
    push!(htmaps, ht)
end
plot(htmaps..., size=(800,800), layout=(3,2), title="R0 = $(pars.R0)", fmt=:png)

Random.seed!(48)
run_and_save_sim(pars; n_sim=5000)

# ## Baseline result for Surveillance part 

# +
path = "../dt_tmp/20230630_162931" # 100 simulations
#path = "../dt_tmp/20230630_163728" # 1000 simulations

path = "../dt_tmp/20230702_213018" # 100 simulations
path = "../dt_tmp/20230702_213608" # 5000 simulations for R0=1.05
#path = "../dt_tmp/20230702_221538" # 5000 simulation for R0=1.10
# -

n_sim = fetch_sim_paths(path) |> length |> println
paths = fetch_sim_paths(path)
res = deserialize(paths[1])
res.pars |> dump

# +
par_AFP = AFPSurParams()
# Hazard rate
p = 0.9
g = - log(1-p)/10
# Catchment area
pop = sp_pars.pop
sum(pop[1:15])/sum(pop)*100
area = fill(0, n_site)
area[1:15] .= 1

par_ES = ESParams(g=g, area=area)
nothing
# -

dump(par_AFP)
dump(par_ES)

include("model_meta_pop.jl")

sim_res = collect_summary_statistics(path, par_AFP, par_ES)
nothing

df = DataFrame(sim_res)
n_sim = size(df)[1]
#h1 = histogram(df[:, :R_final_num], bins=20, title="Final size")
h2 = histogram(df[:, :R_final_site], bins=20, title="Site at least one infection\n# of sim =$n_sim, R0=$(res.pars.R0)", 
    xlabel="Number of sites", ylabel="Count", left_margin=8Plots.mm)
h3 = histogram(df[:, :R_final_AFP], bins=20, title="Final size of AFP", xlabel="AFP cases")
plot(h2, h3,  size=(800,400), fmt=:png)

x = df[:, :R_final_AFP]
histogram(x[x.>0], fmt=:png, xlabel="AFP cases", title="Final size of AFP, excluding 0")

# +
df = DataFrame(sim_res)
#first(df, 5)

outcome_num_prop(df) |> display
detect_pattern(df) |> countmap |> display

dfM = @pipe filter(x -> isnan(x["t_ES"]) == false, df)  
m_ES  = @pipe mean(dfM[:, "t_ES"]) |> round(_, digits=2)
dfM = @pipe filter(x -> isnan(x["t_AFP"]) == false, df)  
m_AFP = @pipe mean(dfM[:, "t_AFP"])  |> round(_, digits=2)
println("Mean ES: $m_ES days, Mean AFP: $m_AFP days")

vis_cumulative_prob(df, pars.days; title="sim =$n_sim, R0=$(res.pars.R0)")
# -

diff = leadtime_diff(df)
leadtime_diff_statistics(diff)

# ## Sensitivity all

include("model_meta_pop.jl")

#path = "../dt_tmp/20230702_213018" # 100 simulations
#path = "../dt_tmp/20230702_213608" # 5000 simulations for R0=1.05
#path = "../dt_tmp/20230702_221538" # 5000 simulation for R0=1.10
n_sim = fetch_sim_paths(path) |> length

df_res1, df_res2, df_res3 = sensitivity_ana_all(path, par_AFP, par_ES)
nothing

df_res = df_res1
tab = detection_pattern_sensitivity(df_res, :pop90)
xticks = (1, 3, 5, 10, 30, 50, 100, 300, 500, 1000)
xlabel = "Pop 90 for coefficient for hazard, g"
vis_detection_pattern(tab, n_sim, :pop90; xscale=:log10, xticks=(xticks, xticks), xlabel=xlabel)
df_diff = leadtime_diff_sensitivity(df_res, :pop90)
vis_leadtime_diff_sensitivity(df_diff, :pop90; xscale=:log10, xticks=(xticks, xticks))

df_res = df_res2
tab = detection_pattern_sensitivity(df_res, :n_freq)
vis_detection_pattern(tab, n_sim, :n_freq; xlabel="Frequency of sampling (days)")
df_diff = leadtime_diff_sensitivity(df_res, :n_freq)
vis_leadtime_diff_sensitivity(df_diff, :n_freq; xlabel="Frequency of sampling (days)")

# +
per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
nothing

df_res = df_res3
tab = detection_pattern_sensitivity(df_res, :ind_site)
tab[:, :per_pop] = per_pop
vis_detection_pattern(tab, n_sim, :ind_site)
vis_detection_pattern(tab, n_sim, :per_pop; xlabel="Population coverage (%)")

df_diff = leadtime_diff_sensitivity(df_res, :ind_site)
vis_leadtime_diff_sensitivity(df_diff, :ind_site)
df_diff[:, :per_pop] = per_pop
vis_leadtime_diff_sensitivity(df_diff, :per_pop;  xlabel="Population coverage (%)")
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






