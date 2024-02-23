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
include("test_utils.jl")

using Test

# # Data preparation

path_f = glob("../data_pop/zaf_f_*.tif")
path_m = glob("../data_pop/zaf_m_*.tif")
path_f |> length |> println
path_m |> length |> println
@test length(path_f) == 18
@test length(path_m) == 18

println(path_f)




