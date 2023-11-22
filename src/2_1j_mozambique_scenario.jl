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
using CSV
using ColorSchemes
using DataFrames
using Dates
using GeoDataFrames
import GeoDataFrames as GDF
using Pipe
using Plots
using Rasters
using Shapefile

include("util.jl")
include("geo_ana.jl")
include("model_meta_pop.jl")
# -

# ## Read WorldPop data

path1 = "../data/zaf_f_5_2020_constrained.tif"
path2 = "../data/zaf_m_5_2020_constrained.tif"
zaf_5_agg = merge_two_map_data(path1, path2)
## Draw population map
#plot(zaf_5_agg, fmt=:png) |> display
#df_zaf = raster_to_df(zaf_5_agg)
#cut_validate_raster_dataframe!(df_zaf, cut_off_pop)
nothing

df_zaf = CSV.read("../res/table_pop_top.csv", DataFrame)
nothing

cut_off_pop = 100

path3 = "../data/moz_f_5_2020_constrained.tif"
path4 = "../data/moz_m_5_2020_constrained.tif"

moz_5_agg = merge_two_map_data(path3, path4)
plot(moz_5_agg, fmt=:png) |> display

df_moz = raster_to_df(moz_5_agg)
cut_validate_raster_dataframe!(df_moz, cut_off_pop)

histogram(log10.(df_zaf[:, :value])) |> display
histogram(log10.(df_moz[:, :value])) |> display

# ## Calculate probability of pi with two districts 

# +

cols = ["lon", "lat", "value", "cnt"]
df_zaf[:, :cnt] .= "zaf"
df_moz[:, :cnt] .= "moz"
sort!(df_zaf, "value", rev=true) # to ensure the order is the same as other files.
sort!(df_moz, "value", rev=true)

df_mer = vcat(df_zaf[:, cols], df_moz)
pop = df_mer[:, :value]
mat = df_mer[:, [:lat, :lon]] |> Matrix

@time π_mat = calculate_probability_of_pi(pop, mat)
nothing
# -

n_zaf = size(df_zaf)[1]
n_moz = size(df_moz)[1]
println("# of sites in zaf: $(n_zaf), in moz: $(n_moz)")
ind_zaf = 1:n_zaf
ind_moz = (n_zaf + 1):(n_zaf + n_moz)

π_inter = π_mat[ind_moz, ind_zaf]
pop_moz = df_moz[:, :value]
prod = π_inter .* pop_moz
imp_ws = sum(prod, dims=1)[1,:] ./ sum(prod)
nothing

println("Imp weights vector length : $(length(imp_ws))")

histogram(log10.(imp_ws))

CSV.write("../data/imp_ws_moz.csv",  Tables.table(imp_ws), writeheader=false)

imp_ws_read = CSV.read("../data/imp_ws_moz.csv", DataFrame, header=false)[:,1]
nothing

dfM = copy(df_zaf)
dfM[:, :imp_ws] = imp_ws
zaf_map = copy(zaf_5_agg)
zaf_map[:] .= 0
for r in eachrow(dfM)
    zaf_map[At(r.lon), At(r.lat)] = log10.(r.imp_ws)
end
pl = plot(zaf_map,
    xlim=[15,35], ylim=[-35, -21],
    axis=nothing, border=:none,
    dpi=100, fmt=:png, size=(800,500),
)
add_zaf_borders!(pl)

# ## Calculate radiation model weighted Airport scenario

df_zaf = CSV.read("../res/table_pop_top.csv", DataFrame)
nothing

pop = df_zaf[:, "value"] 
mat = df_zaf[:, [:lat, :lon]] |> Matrix
π_mat = calculate_probability_of_pi(pop, mat)
nothing

ind = [11, 7, 62] 
weight = [4_342_611, 1_156_996, 188_243]
π_mat[diagind(π_mat)] .= 1/6
imp_prob = sum(π_mat[ind, :] .*weight, dims=1)[1, :]
imp_prob = imp_prob/sum(imp_prob)
CSV.write("../data/imp_ws_airport.csv",  Tables.table(imp_prob), writeheader=false)

# ## Prepare Mozambique risk based simulation 

imp_ws_moz= CSV.read("../data/imp_ws_moz.csv", DataFrame, header=false)[:,1]
imp_ws_airport = CSV.read("../data/imp_ws_airport.csv", DataFrame, header=false)[:,1]
nothing

# +
sort_inds = sortperm(imp_ws_moz, rev=true)
imp_ws_moz_sorted = imp_ws_moz[sort_inds]
CSV.write("../data/imp_ws_moz_sorted.csv",  Tables.table(imp_ws_moz_sorted), writeheader=false)

imp_ws_airport_sorted = imp_ws_airport[sort_inds]
CSV.write("../data/imp_ws_airport_sorted.csv",  Tables.table(imp_ws_airport_sorted), writeheader=false)

df_zaf_sort = df_zaf[sort_inds, :]
nothing
# -

pop = df_zaf_sort[:, :value]
unvac = df_zaf_sort[:, :unvac]
mat = df_zaf_sort[:, [:lat, :lon]] |> Matrix
@time π_mat_moz = calculate_probability_of_pi(pop, mat)
nothing

spatial_p = (pop=pop, unvac=unvac, π_mat=π_mat_moz, df=df_zaf_sort)
path = "../dt_tmp/spatial_params_agg230_moz_sorted.ser"
println(path)
serialize(path, spatial_p)







