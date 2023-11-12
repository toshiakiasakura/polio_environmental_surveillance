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
using SparseArrays

include("util.jl")
include("geo_ana.jl")
include("model_meta_pop.jl")
include("main_model_meta_pop.jl")
# -

path_spatial = "../dt_tmp/spatial_params_agg230.ser"

# ## International airport location, index 

# +
sp_pars = deserialize(path_spatial)
df = sp_pars.df
# Tambo International
lat = -26.12825796201514
lon = 28.242074092511
dist = harversine_dist.(lat, lon, df[:, :lat], df[:, :lon]) 
println("argmin(dist): $(argmin(dist)), dist: $(dist[argmin(dist)])")

# Cape Town 
lat =  -33.970502228847884
lon = 18.600228711334545
dist = harversine_dist.(lat, lon, df[:, :lat], df[:, :lon]) 
println("argmin(dist): $(argmin(dist)), dist: $(dist[argmin(dist)])")

# King Shaka 
lat =  -29.608764960536764
lon = 31.115368797913593
dist = harversine_dist.(lat, lon, df[:, :lat], df[:, :lon]) 
println("argmin(dist): $(argmin(dist)), dist: $(dist[argmin(dist)])")


# -

# # Transmission model 

n_sim = 5000

R0 = 14.0
α = 0.05
pc = 0.25

imp_ws = CSV.read("../data/imp_ws.csv", DataFrame, header=false)[:,1]
nothing

include("model_meta_pop.jl")
include("main_model_meta_pop.jl")

function run_all_patterns(;
        R0=1.0, α=0.05, pc=1.0, 
        imp_ws=[1.0], 
        n_sim=100,
        path_spatial=""
    )
    summary_info = []
    for pattern in ["population_size", "airport", "mozambique"]
        path_trans = run_transmission_model(;
            R0=R0, α=α, pc=pc, imp_ws=imp_ws,
            n_sim=n_sim,
            path_spatial_params=path_spatial, 
            pattern=pattern
        )
        par_AFP, par_ES = set_par_AFP_ES(path_spatial_params=path_spatial, pc=pc)
        path_save = sensitivity_all_summary(par_AFP, par_ES; path_trans=path_trans)
        validate_meta_population_model(path_trans, path_spatial)
        remove_all_leaving_one(path_trans)
        rec = "R0=$(R0), α=$(α), pc=$(pc), n_sim=$(n_sim), pattern=$(pattern), path: $(path_save)"
        push!(summary_info, rec)
    end
    println(rec)
end

run_all_patterns(
    R0=R0, α=α, pc=pc, imp_ws=imp_ws, 
    n_sim=n_sim,
    path_spatial=path_spatial,
    )

# ## Debug field

# +
path_trans = run_transmission_model(;
    R0=R0, α=α, pc=pc, imp_ws=imp_ws,
    n_sim=n_sim,
    path_spatial_params=path_spatial, 
    pattern="population_size"
)
par_AFP, par_ES = set_par_AFP_ES(path_spatial_params=path_spatial, pc=pc)
#baseline_results(par_AFP, par_ES; 
#    path_trans=path_trans
#)

sensitivity_all_summary(par_AFP, par_ES; path_trans=path_trans)
validate_meta_population_model(path_trans, path_spatial)
remove_all_leaving_one(path_trans)
# -

path_trans = run_transmission_model(;
    R0=R0, α=α, pc=pc, n_sim=n_sim,
    path_spatial_params=path_spatial, 
    pattern="airport"
)
par_AFP, par_ES = set_par_AFP_ES(path_spatial_params=path_spatial)
sensitivity_all_summary(par_AFP, par_ES; path_trans=path_trans)
validate_meta_population_model(path_trans, path_spatial)
remove_all_leaving_one(path_trans)

path_trans = run_transmission_model(;
    R0=R0, α=α, pc=pc, imp_ws=imp_ws,
    n_sim=n_sim,
    path_spatial_params=path_spatial, 
    pattern="mozambique"
)
par_AFP, par_ES = set_par_AFP_ES(path_spatial_params=path_spatial)
sensitivity_all_summary(par_AFP, par_ES; path_trans=path_trans)
validate_meta_population_model(path_trans, path_spatial)
remove_all_leaving_one(path_trans)












