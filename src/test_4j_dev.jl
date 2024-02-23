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

using Test

# ## Data preparation

# +
imp_ws_moz = CSV.read("../data/imp_ws_moz.csv", DataFrame, header=false)[:,1]
imp_ws_airport = CSV.read("../data/imp_ws_airport.csv", DataFrame, header=false)[:,1]

imp_ws_moz_sorted = CSV.read("../data/imp_ws_moz_sorted.csv", DataFrame, header=false)[:,1]
imp_ws_airport_sorted = CSV.read("../data/imp_ws_airport_sorted.csv", DataFrame, header=false)[:,1]
nothing
# -

# ## Model notation
#

sp_pars = read_spatial_params_file("ES_population_size")
nothing

π_mat = sp_pars.π_mat
# Check a digonal element of π is 0.0.
@test all(diag(π_mat) .== 0.0)
nothing

π_tmp = [0 0.1 0.5; 0 0 0.5; 0 0 0 ]
I = [10, 20, 30]
@test sum(π_tmp' .* I', dims=2)[:, 1] == [0, 1, 15]

# # Test a single site model with a meta-population model. 

w_pop_cov = 0.9196
Nc = 100_000
N_unvac = round(Nc * (1-w_pop_cov) ) |> Int64
pars = SEIRMetaModelParams(
    R0 = 14.0, α=0.05, pc=0.25,
    Nc = [Nc, Nc], N_unvac=[N_unvac, N_unvac],
    n_site=2,
    π_mat = [0 0; 0 0],
    imp_ws=[1, 0]
)


Random.seed!(48)
path_trans = run_and_save_sim(pars; 
    n_sim=10, pattern="population_size"
)

path_lis = glob(path_trans * "/*")

pl = plot()
for (i, path) in enumerate(path_lis)
    @unpack rec = load(path)["data"]
    x = length(rec.Ic[1, :])
    n_tot = sum(rec.Ic[1, :])
    plot!(pl, 1:x, rec.Ic[1,:]/pars.pc, label="y$(i), $(n_tot)")
end
pl

pl = plot()
for (i, path) in enumerate(path_lis)
    @unpack rec = load(path)["data"]
    x = length(rec.Ic[1, :])
    n_tot = sum(rec.Z_A5_6)
    plot!(pl, 1:x, rec.Z_A5_6, label="$(i), $(n_tot)")
end
pl

remove_all_transmission_results(path_trans)


