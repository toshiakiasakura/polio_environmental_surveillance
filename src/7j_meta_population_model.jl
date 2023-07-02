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

sp_pars = deserialize("../dt_tmp_fix/pop_pi_agg1000_10unvac.ser")
keys(sp_pars)

n_site = length(sp_pars.pop)
nothing

pars = SEIRMetaModelParams(
    R0=2.0,
    N0=sp_pars.pop,
    π_mat=sp_pars.π_mat,
    n_site=n_site,
    α=0.050
    )
pars |> dump

include("model_meta_pop.jl")

# ## Run simulations 

Random.seed!(45)
@time rec, outcome, pars = run_sim(pars; rec_flag=true)
println(outcome.R_final_num, ", ", outcome.R_final_site)
heatmap_meta_pop(rec.Ia)

@pipe rec.H_AFP |> sum

@time res = run_sim(pars; rec_flag=true)
println("$t_AFP, $(outcome.t_extinct), $t_ES")

# +

run_and_save_sim(pars; n_sim=100)
# -

# ## Surveillance part 

include("model_meta_pop.jl")

path = "../dt_tmp/20230630_162931" # 100 simulations
#path = "../dt_tmp/20230630_163728" # 1000 simulations

p = 0.9
λ0 = - log(1-p)/10

par_AFP = AFPSurParams()
area = fill(0, n_site)
area[1:2] .= 1
par_ES = ESParams(λ0=λ0, area=area)
nothing

sim_res = collect_summary_statistics(path, par_AFP, par_ES)
nothing





df = DataFrame(sim_res)
#h1 = histogram(df[:, :R_final_num], bins=20, title="Final size")
h2 = histogram(df[:, :R_final_site], bins=20, title="Site at least one infection")
h3 = histogram(df[:, :R_final_AFP], bins=20, title="Final size of AFP")
plot(h2, h3,  size=(800,400))

df = DataFrame(sim_res)
#first(df, 5)
outcome_num_prop(df) |> display
detect_pattern(df) |> countmap |> display
vis_cumulative_prob(df, pars.days)

diff = leadtime_diff(df)
leadtime_diff_statistics(diff)

# ## Sensitivity analysis for ES catchment area

include("model_meta_pop.jl")

per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
nothing

df_res = sensitivity_catchment_area(path, par_AFP, par_ES)
nothing

n_sim = fetch_sim_paths(path) |> length

tab = detection_pattern_sensitivity(df_res, :ind_site)
tab[:, :per_pop] = per_pop
vis_detection_pattern(tab, n_sim, :ind_site)
vis_detection_pattern(tab, n_sim, :per_pop)

include("model_meta_pop.jl")

df_diff = leadtime_diff_sensitivity(df_res, :ind_site)
vis_leadtime_diff_sensitivity(df_diff, :ind_site)
df_diff[:, :per_pop] = per_pop
vis_leadtime_diff_sensitivity(df_diff, :per_pop)


