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
# -

# # Vaccination coverage data

# +
path = "../data/zaf_OPV_HEXA_vaccine_coverage_2020.csv"
df_vac = CSV.read(path, DataFrame)
VE = 0.63
VE1 = 1 - (1- VE)
VE2 = 1 - (1 - VE)^2
VE3 = 1 - (1 - VE)^3
VE4 = 1 - (1 - VE)^4

CV4 = df_vac.HEXA4
CV3 = df_vac.HEXA3
CV2 = df_vac.HEXA2 
CV1 = df_vac.HEXA1 

dif3 = [i < 0 ? 0.0 : i for i in CV3 .- CV4]
dif2 = [i < 0 ? 0.0 : i for i in CV2 .- CV3]
dif1 = [i < 0 ? 0.0 : i for i in CV1 .- CV2]
EVP = CV4.*VE4 .+ dif3 .* VE3 .+ dif2 .* VE2 .+ dif1 .* VE1
df_vac[!, "EVP"] = EVP
nothing
# -

describe(df_vac[:, "EVP"])

# +
#path_shape = "../data/geoBoundaries-ZAF-ADM2.shp"
path_shape = "../dt_geoBoundaries-ZAF-ADM2-all/geoBoundaries-ZAF-ADM2.shp"
df_geo = GDF.read(path_shape)
df_geo = innerjoin(df_geo, df_vac[:, ["shapeName", "EVP"]], on="shapeName")

l = @layout [a{0.95w} b]
pl = plot(axis=nothing, border=:none, dpi=300)
lower = 0.75
amp = reverse(ColorSchemes.amp)
for r in eachrow(df_geo)
    sc = (r.EVP - lower)/(1 - lower)
    plot!(pl, r.geometry, color=amp[sc])
end
p2 = heatmap(rand(2,2), clims=(lower,1), framestyle=:none, 
    c=cgrad(amp), cbar=true, lims=(-1,0),
)
pl = plot(pl, p2, layout=l, right_margin=10Plots.mm, fmt=:png)
annotate!((0.1, 0.95), "(B)")
display(pl)
# savefig(pl, "../res/fig_vaccine_coverage.png")
# -

df_save = copy(df_vac)
cols = ["shapeName", "sample_size",  "OPV0", "OPV1",
 "HEXA1",  "HEXA2",  "HEXA3",  "HEXA4",  "EVP", ]
covs = ["OPV0", "OPV1", "HEXA1",  "HEXA2",  "HEXA3",  "HEXA4",  "EVP", ]
df_save = df_save[:, cols]
df_save[:, covs] = df_save[:, covs] .* 100
df_save[:, "EVP"] = round.(df_save[:, "EVP"], digits=3)
CSV.write("../res/table_vaccine_coverage.csv", df_save)

# ## Read WorldPop data

path = "../data/zaf_f_5_2020_constrained.tif"
f_5_zaf = read(Raster(path))
plot(f_5_zaf)

path = "../data/zaf_m_5_2020_constrained.tif"
m_5_zaf = read(Raster(path))
m_5_zaf = replace_missing(m_5_zaf, 0)
plot(m_5_zaf)

# #### Edit Point

agg_scale = 230 # Latitude diff: 21.31 km, Longitude diff: 19.74 km
# agg_scale = 110 # Latitude diff: 10.19 km, Longitude diff: 9.44 km
m_5_zaf_agg = Rasters.aggregate(sum, m_5_zaf, agg_scale; skipmissingval=true)
nothing

# Print basic info.
nr, nc = size(m_5_zaf_agg)
println("row: $nr, col: $nc, grid num: $(nr*nc)")
get_gridsize(m_5_zaf_agg)

# Draw population map
plot(m_5_zaf_agg, fmt=:png) |> display

include("util.jl")

# Check population size change.
m_5_zaf |> sum |> println
m_5_zaf_agg |> sum |> println

# +
cut_off_pop = 100
df_ras = raster_to_df(m_5_zaf_agg)

println("Original size: ",size(df_ras))
df_ras[!, :value] = @pipe df_ras[:, :value] .|> round(_; digits=0)
filter!(x -> x.value > 0.0, df_ras)
n_bf = df_ras[:, :value] |> sum
println("Total population size before removing: ", n_bf)
println("Before removing: ", size(df_ras))
filter!(x -> x.value > cut_off_pop, df_ras)
n_af = df_ras[:, :value] |> sum
println("Total population size after removing: ", n_af)
println("After removing: ", size(df_ras))
nothing
# -

println("% or removal: ", (n_bf - n_af)/n_bf * 100)

bins = [2 + 0.25*i for i in 0:12]
pl = histogram(log10.(df_ras[:, :value]), 
    bins=bins,
    ylabel="Number of locations",
    xlabel="log10(Population size)",
    legend=:none, dpi=300, fmt=:png
)
savefig(pl, "../res/fig_pop_histogram.png")
display(pl)

zaf_map = copy(m_5_zaf_agg)
cond = zaf_map .< cut_off_pop
zaf_map[cond] .= 0
pl = plot(zaf_map, 
    xlim=[15,35], ylim=[-35, -21],
    axis=nothing, border=:none,
    right_margins=8Plots.mm,
    dpi=300, fmt=:png
)
annotate!((0.1, 0.95), "(A)")
add_zaf_borders!(pl)
display(pl)
# savefig(pl, "../res/fig_pop_map.png")

# ## Relate vaccination coverage data to population data. 

df_ras  # population with coordinates
df_geo # vaccination coverage by district with polygon data.
nothing

dist = []
for r in eachrow(df_ras)
    point = ArchGDAL.createpoint(r.lon, r.lat)
    flag = false
    for g in eachrow(df_geo)
        pol = g.geometry
        if  ArchGDAL.within(point, pol)
            push!(dist, g.shapeName)
            flag = true
            break
        end
    end
    if flag == false
        push!(dist, "not determined")
    end
end
df_ras[!, "shapeName"] = dist
nothing

(df_ras[:, "shapeName"] .== "not determined") |> sum

df_fil = filter(x -> x.shapeName == "not determined", df_ras)
pl = plot()
scatter!(pl, df_fil[:, :lon], df_fil[:, :lat], )
add_zaf_borders!(pl)

function calculate_minimum_dist_given_polygon(lon, lat, pol)
    pol_sim = ArchGDAL.simplify(pol, 0.05)
    boundary = ArchGDAL.boundary(pol_sim)
    n = ArchGDAL.ngeom(boundary)
    dist = Inf
    for i in 1:n
        lon2, lat2, z = ArchGDAL.getpoint(boundary, i - 1)
        dist_tmp = harversine_dist(lat, lon, lat2, lon2)
        dist = minimum([dist, dist_tmp])
    end
    return dist
end

for (i,r) in enumerate(eachrow(df_ras))
    if r.shapeName == "not determined"
        dist = Inf
        dist_name = "not determined"
        for r_geo in eachrow(df_geo)
            dist_tmp = calculate_minimum_dist_given_polygon(r.lon, r.lat, r_geo.geometry)
            if dist > dist_tmp
                dist = dist_tmp
                dist_name = r_geo.shapeName
            end
        end
        df_ras[i, "shapeName"] = dist_name
    end
end

(df_ras[:, "shapeName"] .== "not determined") |> sum |> println
df_mer = leftjoin(df_ras, df_vac[:, ["shapeName", "EVP"]], on="shapeName")
nothing

df_mer[!, :unvac] = @pipe (df_mer[:, :value] .* (1 .- df_mer[:, :EVP])) .|> round(_, digits=0)
sort!(df_mer, "value", rev=true)
df_save = copy(df_mer)
df_save[:, :EVP] = df_save[:, :EVP] .* 100
CSV.write("../res/table_pop_top.csv", df_save)
first(df_mer, 15)

# Weighted vaccine coverage. 
prop = df_mer[:, :value]./sum(df_mer[:, :value])
weighted_EVP = (df_mer[:, :EVP] .* prop) |> sum

# ## Calculate the probability of mobilisations

include("util.jl")

# +
pop = df_mer[:, "value"] 
unvac = @pipe df_mer[:, :unvac] .|> round(_, digits=0) .|> Int64
mat = df_mer[:, [:lat, :lon]] |> Matrix
π_mat = calculate_probability_of_pi(pop, mat)

spatial_p = (pop=pop, unvac=unvac, π_mat=π_mat, df=df_mer)
path = "../dt_tmp/spatial_params_agg$(agg_scale).ser"
println(path)
serialize(path, spatial_p)
# -



# ## Visualisation

include("util.jl")

function visualise_probs(ind::Int64, title::String)
    zaf_map = copy(m_5_zaf_agg)
    zaf_map[:] .= 0
    n = size(df_mer)[1]
    for i in 1:n
        x = df_mer[i, :lon]
        y = df_mer[i, :lat]
        zaf_map[ At(x), At(y)] = log10(π_mat[i, ind])
    end
    pl = plot()
    plot!(zaf_map,
        xlim=[15,35], ylim=[-35, -21],
        axis=nothing, border=:none,
        logscale=:log10, 
        #colorbar_title="log10(πij)",
        title=title,
    )
    add_zaf_borders!(pl)
    return pl
end

# +
ind = 1
pl1 = histogram(log10.(π_mat[:, ind]),
    legend=:none, 
    xlabel = "log10(πi1)", ylabel="Number of locations",
)
annotate!(pl1, (0.1, 0.95), "(A)")
dist = df_mer[ind, "shapeName"]
pl2 = visualise_probs(ind, "1st populous location\n$(dist)")
annotate!(pl2, (0.1, 0.95), "(B)")

ind = 2
pl3 = histogram(log10.(π_mat[:, ind]),
    legend=:none, 
    xlabel = "log10(πi2)", ylabel="Number of locations",
)
annotate!(pl3, (0.1, 0.95), "(C)")
dist = df_mer[ind, "shapeName"]
pl4 = visualise_probs(ind, "2nd populous location\n$(dist)")
annotate!(pl4, (0.1, 0.95), "(D)")

ind = 3
pl5 = histogram(log10.(π_mat[:, ind]),
    legend=:none, 
    xlabel = "log10(πi3)", ylabel="Number of locations",
)
annotate!(pl5, (0.1, 0.95), "(E)")
dist = df_mer[ind, "shapeName"]
pl6 = visualise_probs(ind, "3rd populous location\n$(dist)")
annotate!(pl6, (0.1, 0.95), "(F)")

l = @layout [
    a{0.3w} b
    c{0.3w} d
    e{0.3w} f
]
pl = plot(pl1, pl2, pl3, pl4, pl5, pl6,
    dpi=300, size=(800, 350 * 3), 
    layout=l, left_margin=5Plots.mm,
)
display(pl)
savefig(pl, "../res/fig_pi_map.png")
# -






