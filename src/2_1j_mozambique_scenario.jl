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

# ## Read WorldPop data

path = "../data_pop/zaf_merged_0_4.tif"
zaf_0_4 = read(Raster(path))
path = "../data_pop/moz_merged_0_4.tif"
moz_0_4 = read(Raster(path))
nothing

df_zaf = CSV.read("../res/table_pop_top.csv", DataFrame)
nothing

plot(moz_0_4, fmt=:png) |> display

cut_off_pop = 100
df_moz = raster_to_df(moz_0_4)
cut_validate_raster_dataframe!(df_moz, cut_off_pop)

histogram(log10.(df_zaf[:, :value]), title="zaf") |> display
histogram(log10.(df_moz[:, :value]), title="moz") |> display

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

histogram(log10.(imp_ws), title="Importation prob. from Mozambique")

# a file for importation weights.
CSV.write("../data/imp_ws_moz.csv",  Tables.table(imp_ws), writeheader=false)

imp_ws_read = CSV.read("../data/imp_ws_moz.csv", DataFrame, header=false)[:,1]
nothing

dfM = copy(df_zaf)
dfM[:, :imp_ws] = imp_ws
zaf_map = copy(zaf_0_4)
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

include("utils.jl")

df_zaf = CSV.read("../res/table_pop_top.csv", DataFrame)
nothing

sp_pars = read_spatial_params_file("ES_population_size")
nothing

pop = sp_pars.df[:, "value"]
mat = sp_pars.df[:, [:lat, :lon]] |> Matrix
π_mat = sp_pars.π_mat
df = sp_pars.df
nothing

# +
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



# +
ind = [10, 7, 61]
weight = [4_342_611, 1_156_996, 188_243] # Travel volume
π_mat[diagind(π_mat)] .= 1/6
imp_prob = sum(π_mat[ind, :] .*weight, dims=1)[1, :]
imp_prob = imp_prob/sum(imp_prob)

# Importation probability in the airport introduction scenario. 
CSV.write("../data/imp_ws_airport.csv",  Tables.table(imp_prob), writeheader=false)
# -

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
path = "../dt_tmp/spatial_params_agg230_moz_sorted.jld2"
println(path)
jldsave(path; data=spatial_p)







