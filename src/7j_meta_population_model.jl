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

p = 0.9
λ0 = - log(1-p)/10

sp_params = deserialize("../dt_tmp_fix/pop_pi_agg1000_10unvac.ser")
keys(sp_params)

n_site = length(sp_params.pop)
ES_area = fill(0, n_site)
ES_area[1:4] .= 1
nothing

params = SEIRMetaModelParams(
    R0=3.0,
    λ0=λ0,
    N0=sp_params.pop,
    π_mat=sp_params.π_mat,
    n_site=n_site,
    ES_area=ES_area,
    )
params |> dump

# ### TODO list.
# TODO: Update `update_model`  
# TODO: Change AFP surveillance and ES surveillance appropriately  
# TODO: Write the meta-SEIR model including the full model.   
# TODO: Run the model.  












