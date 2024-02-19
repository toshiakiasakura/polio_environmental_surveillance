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
include("geo_ana.jl")
include("model_meta_pop.jl")
include("main_model_meta_pop.jl")

# ## Data preparation

# +
imp_ws_moz = CSV.read("../data/imp_ws_moz.csv", DataFrame, header=false)[:,1]
imp_ws_airport = CSV.read("../data/imp_ws_airport.csv", DataFrame, header=false)[:,1]

imp_ws_moz_sorted = CSV.read("../data/imp_ws_moz_sorted.csv", DataFrame, header=false)[:,1]
imp_ws_airport_sorted = CSV.read("../data/imp_ws_airport_sorted.csv", DataFrame, header=false)[:,1]
nothing
# -

# # Test, single model run 

R0 = 14.0
imp_ws = [1.0]

path_trans = run_transmission_model(
    R0=14.0, α=0.05, pc=0.25, imp_ws=[1.0],
    n_sim=10, pattern="population_size"
)

par_AFP, par_ES = set_par_AFP_ES(;pc=0.25, pattern="population_size")
nothing

path_objs = fetch_sim_paths(path_trans)
path = path_objs[1]
res = load(path)["data"]
res |> typeof

path_save = save_sensitivity_ES_catchment_area(
    par_AFP, par_ES, path_trans; 
    inc_prop=0.05
)

path_save = save_sensitivity_sampling_frequency(
    par_AFP, par_ES, path_trans; 
)

remove_all_leaving_one(path_trans)

function post_processing(path_trans, path_spatial, pc)
    par_AFP, par_ES = set_par_AFP_ES(path_spatial_params=path_spatial, pc=pc)
    path_save = sensitivity_all_summary(par_AFP, par_ES; path_trans=path_trans)
    #validate_meta_population_model(path_trans, path_spatial)
    remove_all_leaving_one(path_trans)
    return path_save
end





    remove_all_leaving_one(path_trans)
