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

using ArchGDAL
using Rasters, Dates
using DataFrames
using Pipe
using Plots

# ## Read WorldPop data

path = "../data/zaf_f_5_2020_constrained.tif"
f_5_zaf = read(Raster(path))
plot(f_5_zaf)

path = "../data/zaf_m_5_2020_constrained.tif"
m_5_zaf = read(Raster(path))
m_5_zaf = replace_missing(m_5_zaf, 0)
plot(m_5_zaf)

m_5_zaf_agg = Rasters.aggregate(sum, m_5_zaf, 200; skipmissingval=true)
plot(m_5_zaf_agg, fmt=:png) |> display
nr, nc = size(m_5_zaf_agg)
nr*nc

(37.2912 - 16.4579)/129 * 111 # 200 aggregation, 19211 grids -> eventually around 1500-2000
# > 17.92 km 
(37.2912 - 16.4579)/51 * 111 # 500 aggregation, 3009 grids -> eventually around 500-800 
# > 45.34 km 

m_5_zaf_agg = Rasters.aggregate(sum, m_5_zaf, 1000; skipmissingval=true)

m_5_zaf |> sum |> println
m_5_zaf_agg |> sum |> println

size(m_5_zaf_agg)

plot(m_5_zaf_agg, fmt=:png)

coords = []
lon = lookup(m_5_zaf_agg, X)
lat = lookup(m_5_zaf_agg, Y)
for lo in lon, la in lat
    push!(coords, (lon=lo, lat=la, value=m_5_zaf_agg[At(lo), At(la)]))
end
df = DataFrame(coords)
first(df, 5) |> display
df |> size

df[!, :value] = @pipe df[:, :value] .|> round(_; digits=0)
cond = df[:, :value] .>= 100
df[cond, :value] |> size

histogram(log10.(df[cond, :value]), ylabel=:count,xlabel="log10(pop_size)")

# ## Basic usage

lon, lat = X(25:1:30), Y(25:1:30)
ti = Ti(DateTime(2001):Month(1):DateTime(2002))
ras = Raster(rand(lon, lat, ti))

lon = lookup(ras, X) # if X is longitude
lat = lookup(ras, Y) # if Y is latitude

ras[Ti(8)]

using Rasters, RasterDataSources, ArchGDAL, Plots

A = Raster(WorldClim{BioClim}, 5)
madagascar = view(A, X(43.25 .. 50.48), Y(-25.61 .. -12.04)) # Note the space between .. -12
plot(madagascar)



madagascar

lon = lookup(madagascar, X)

ismissing.(madagascar) |> sum

plot(madagascar .* 5)

# ## Modelling 

lon = lookup(madagascar, X) # if X is longitude
lat = lookup(madagascar, Y) # if Y is latitude

lo = lon[20]
la = lat[20]
madagascar[At(lo), At(la)]

coords = []
for lo in lon, la in lat
    push!(coords, (lon=lo, lat=la, value=madagascar[At(lo), At(la)]))
end
df = DataFrame(coords)
first(df, 5)

# This is matrix for corresponding values.
mat = unstack(df, :lon, :lat, :value)[:, 2:end] |> Matrix
nothing

# ## Aggregate 

plot(madagascar)

using Rasters: Center
ag = Rasters.aggregate(sum, madagascar, (X(6), Y(4)); skipmissingval=true)
size(madagascar) |> println
size(ag) |> println
plot(ag)

ag[1,1] |> typeof

replace_missing(madagascar , 0) |> sum |> println
replace_missing(ag, 0) |> sum |> println

# ### Data manipulation 

ag_new = ag
ag_new[ag .< 700] .= ag[1,1]
plot(ag_new)


