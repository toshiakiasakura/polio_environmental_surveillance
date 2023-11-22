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

# # Vaccination coverage data

df_vac = obtain_EVP()
nothing

describe(df_vac[:, "EVP"])

path_shape = "../dt_geoBoundaries-ZAF-ADM2-all/geoBoundaries-ZAF-ADM2.shp"
df_geo = GDF.read(path_shape)
df_geo = innerjoin(df_geo, df_vac[:, ["shapeName", "EVP"]], on="shapeName")
nothing

create_vaccination_coverage_map(df_geo)

save_vaccination_coverage_data(df_vac)

g = 0.23025
I = 1:1000
w = 1 .- exp.(-0.97.* g.*I)
plot(log10.(I), w*100, xticks=([0,1,2,3],[1,10,100,1000]), label=false,
    xlabel="Number of infections, I(t)", ylabel="Probability of detection, wi,t",
    fmt=:png, fontsize=22
)

# ## Read WorldPop data

path1 = "../data/zaf_f_5_2020_constrained.tif"
path2 = "../data/zaf_m_5_2020_constrained.tif"

pl = plot(framestyle=:none, fmt=:png, figsize=(800,600))
add_zaf_borders!(pl)



zaf_5_agg = merge_two_map_data(path1, path2)
nothing

# #### Edit Point

# Draw population map
plot(zaf_5_agg, fmt=:png) |> display

cut_off_pop = 100
df_ras = raster_to_df(zaf_5_agg)
cut_validate_raster_dataframe!(df_ras, cut_off_pop)

bins = [2 + 0.25*i for i in 0:12]
pl = histogram(log10.(df_ras[:, :value]), 
    bins=bins,
    ylabel="Number of locations",
    xlabel="log10(Population size)",
    legend=:none, dpi=300, fmt=:png
)
savefig(pl, "../res/fig_pop_histogram.png")
display(pl)

create_map_after_cutting(zaf_5_agg, cut_off_pop)

# ## Relate vaccination coverage data to population data. 

include("geo_ana.jl")

df_ras  # population with coordinates
df_geo # vaccination coverage by district with polygon data.
nothing

df_geo |> typeof

df_ras = relate_df_ras_to_district_info(df_ras, df_geo)
nothing

n = (df_ras[:, "shapeName"] .== "not determined") |> sum 
println("Check unclassified points are zero: ", n)

# df_mer: df_ras + EVP and unvaccinated population info.
df_mer = leftjoin(df_ras, df_vac[:, ["shapeName", "EVP"]], on="shapeName")
df_mer[!, :unvac] = @pipe (df_mer[:, :value] .* (1 .- df_mer[:, :EVP])) .|> round(_, digits=0)
sort!(df_mer, "value", rev=true)
pop = df_mer[:, :value]
df_mer[!, :cum_per] = cumsum(pop)./sum(pop).*100 
nothing

# Save top30 population district data.
df_save = copy(df_mer)
df_save[!, :cum_prop] = @pipe cumsum(df_save[:, :value])/sum(df_save[:, :value]).*100 .|> round(_, digits=1)
df_save[:, :EVP] = df_save[:, :EVP] .* 100
CSV.write("../res/table_pop_top.csv", df_save)
first(df_mer, 30)

# Weighted vaccine coverage. 
prop = df_mer[:, :value]./sum(df_mer[:, :value])
weighted_EVP = (df_mer[:, :EVP] .* prop) |> sum

# ### Visualise the unimmunised population

visualise_unimmune(df_mer, zaf_5_agg, "")

# ### Calculate population and ES coverage by district

# +
national_ES_cov = 8.59
df_prop = copy(df_mer)
df_prop[:, :flag_pc100] = (df_prop[:, :cum_per] * 1.0) .< national_ES_cov
df_prop[:, :flag_pc50] = (df_prop[:, :cum_per] * 0.5) .< national_ES_cov
df_prop[:, :flag_pc25] = (df_prop[:, :cum_per] * 0.25) .< national_ES_cov
df_prop[:, :flag_pc20] = (df_prop[:, :cum_per] * 0.20) .< national_ES_cov

prop_f = (x,y; pc=1.0) -> sum(x .* y) ./ sum(x) .* 100 .* pc
dfM = @pipe df_prop |> groupby(_, :shapeName) |> 
    DataFrames.combine(_, :value =>( x -> sum(x) ) => :tot_pop, 
        [:value, :flag_pc100] => prop_f => :prop_pc100,
        [:value, :flag_pc50] =>( (x,y) -> prop_f(x,y; pc=0.5) )=> :prop_pc50,
        [:value, :flag_pc25] =>( (x,y) -> prop_f(x,y; pc=0.25) )=> :prop_pc25,
        #[:value, :flag_pc20] =>( (x,y) -> prop_f(x,y; pc=0.2) )=> :prop_pc20,
    )
@pipe dfM |> sort(_, :prop_pc25, rev=true) |> first(_, 20)
# -

# ### Visualise where is the main points 
# TODO: Move this part to another code file, and create the gif picture. 

cum_per = df_mer[:, :cum_per]
"$(round(cum_per[1],digits=1))"

# +
# Prepare data
cum_per = df_mer[:, :cum_per]
pop = df_mer[:, :value]
sens_index = obtain_ES_sensitivity_index(pop, 0.01)
zaf_map = copy(zaf_5_agg)
zaf_map[:] .= 0

# Create multiple files.
l_ind =  1
for (i, r_ind) in enumerate(sens_index)
    ind = l_ind:r_ind
    for r in eachrow(df_mer[ind,:])
        zaf_map[ At(r.lon), At(r.lat)] = 100
    end
    l_ind = r_ind + 1
    #annotate!(pl, r.lon, r.lat, text(r.cum_per, :left, 4))
    pl = plot(zaf_map,
        xlim=[15,35], ylim=[-35, -21],
        color=:blue,
        colorbar=false,
        axis=nothing, border=:none,
        dpi=100, fmt=:png, size=(1200,1200),
    )
    per = round(cum_per[r_ind],digits=1)
    annotate!(20, -20, text("%Population size: $(per)%\nNumber of sites: $(r_ind)", :left, 30))
    add_zaf_borders!(pl)
    file = lpad(i, 4, "0")
    savefig(pl, "../dt_tmp/gif_map/$(file).png")
end
# Make gif in python
# -

sum(sens_index .<= 33) 

n = (df_mer[:, :cum_per] .<= 8.59/0.25) |> sum
println("pc 25, maximum index: ", sum(sens_index .<= n), ", Number of sites: ", n)
n = (df_mer[:, :cum_per] .<= 8.59/0.30) |> sum
println("pc 30, maximum index: ", sum(sens_index .<= n), ", Number of sites: ", n)
n = (df_mer[:, :cum_per] .<= 8.59/0.5) |> sum
println("pc 50, maximum index: ", sum(sens_index .<= n), ", Number of sites: ", n)

# ## Calculate the probability of mobilisations

include("util.jl")

# +
df_mer = sort(df_mer, "value", rev=true)
pop = df_mer[:, "value"] 
unvac = @pipe df_mer[:, :unvac] .|> round(_, digits=0) .|> Int64
mat = df_mer[:, [:lat, :lon]] |> Matrix
π_mat = calculate_probability_of_pi(pop, mat)

spatial_p = (pop=pop, unvac=unvac, π_mat=π_mat, df=df_mer)
path = "../dt_tmp/spatial_params_agg$(agg_scale).ser"
println(path)
serialize(path, spatial_p)
# -
pop_rev = sort(pop, rev=false)
prop = pop_rev./sum(pop_rev)
cum = cumsum(prop)
scatter(log10.(pop_rev), cum, title="How population_size scenairio works")

histogram(log10.(pop))

# ## Visualisation

include("util.jl")

function visualise_probs(df_mer, ind::Int64, title::String)
    zaf_map = copy(zaf_5_agg)
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
pl1 = histogram(log10.(π_mat[ind, :]),
    legend=:none, 
    xlabel = "log10(π1i)", ylabel="Number of locations",
)
annotate!(pl1, (0.1, 0.95), "(A)")
dist = df_mer[ind, "shapeName"]
pl2 = visualise_probs(df_mer, ind, "1st populous location\n$(dist)")
annotate!(pl2, (0.1, 0.95), "(B)")

ind = 2
pl3 = histogram(log10.(π_mat[ind, :]),
    legend=:none, 
    xlabel = "log10(π2i)", ylabel="Number of locations",
)
annotate!(pl3, (0.1, 0.95), "(C)")
dist = df_mer[ind, "shapeName"]
pl4 = visualise_probs(df_mer, ind, "2nd populous location\n$(dist)")
annotate!(pl4, (0.1, 0.95), "(D)")

ind = 3 
pl5 = histogram(log10.(π_mat[ind, :]),
    legend=:none, 
    xlabel = "log10(π3i)", ylabel="Number of locations",
)
annotate!(pl5, (0.1, 0.95), "(E)")
dist = df_mer[ind, "shapeName"]
pl6 = visualise_probs(df_mer, ind, "3rd populous location\n$(dist)")
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
# ## Unvaccinated base 

# +
df_mer2 = sort(df_mer, "unvac"; rev=true)
pop = df_mer2[:, "value"] 
unvac = @pipe df_mer2[:, :unvac] .|> round(_, digits=0) .|> Int64
mat = df_mer2[:, [:lat, :lon]] |> Matrix
π_mat = calculate_probability_of_pi(pop, mat)

spatial_p = (pop=pop, unvac=unvac, π_mat=π_mat, df=df_mer2)
path = "../dt_tmp/spatial_params_agg$(agg_scale)_unvac.ser"
println(path)
serialize(path, spatial_p)
# -
first(df_mer2, 10)

# +
ind = 1
pl1 = histogram(log10.(π_mat[:, ind]),
    legend=:none, 
    xlabel = "log10(πi1)", ylabel="Number of locations",
)
annotate!(pl1, (0.1, 0.95), "(A)")
dist = df_mer2[ind, "shapeName"]
pl2 = visualise_probs(df_mer2, ind, "1st populous location\n$(dist)")
annotate!(pl2, (0.1, 0.95), "(B)")

ind = 2
pl3 = histogram(log10.(π_mat[:, ind]),
    legend=:none, 
    xlabel = "log10(πi2)", ylabel="Number of locations",
)
annotate!(pl3, (0.1, 0.95), "(C)")
dist = df_mer2[ind, "shapeName"]
pl4 = visualise_probs(df_mer2, ind, "2nd populous location\n$(dist)")
annotate!(pl4, (0.1, 0.95), "(D)")

ind = 3
pl5 = histogram(log10.(π_mat[:, ind]),
    legend=:none, 
    xlabel = "log10(πi3)", ylabel="Number of locations",
)
annotate!(pl5, (0.1, 0.95), "(E)")
dist = df_mer2[ind, "shapeName"]
pl6 = visualise_probs(df_mer2, ind, "3rd populous location\n$(dist)")
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
savefig(pl, "../res/fig_pi_map_unvac.png")
# -






