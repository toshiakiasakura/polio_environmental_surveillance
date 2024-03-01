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

# # Create population map data

path1 = "../data_pop/zaf_f_0_2020_constrained.tif"
path_lis = [
    "../data_pop/zaf_f_1_2020_constrained.tif",
    "../data_pop/zaf_m_0_2020_constrained.tif",
    "../data_pop/zaf_m_1_2020_constrained.tif",
]
map_sum = read_agg_map(path1)
for path in path_lis
    map_sum += read_agg_map(path)
end
write("../data_pop/zaf_merged_0_4.tif", map_sum, force=true)

path1 = "../data_pop/moz_f_0_2020_constrained.tif"
path_lis = [
    "../data_pop/moz_f_1_2020_constrained.tif",
    "../data_pop/moz_m_0_2020_constrained.tif",
    "../data_pop/moz_m_1_2020_constrained.tif",
]
map_sum = read_agg_map(path1)
for path in path_lis
    map_sum += read_agg_map(path)
end
write("../data_pop/moz_merged_0_4.tif", map_sum, force=true)

path_f = glob("../data_pop/zaf_f_*.tif")
path_m = glob("../data_pop/zaf_m_*.tif")
path_f |> length |> println
path_m |> length |> println

map_sum = read_agg_map(path_f[1])
@showprogress for path in vcat(path_f[2:end], path_m)
    map_sum += read_agg_map(path)
end
write("../data_pop/zaf_merged_whole_pop.tif",
    map_sum, force=true)

# # Vaccination coverage data

df_vac = obtain_EVP()
nothing

describe(df_vac[:, "EVP"])

path_shape = "../dt_geoBoundaries-ZAF-ADM2-all/geoBoundaries-ZAF-ADM2.shp"
df_geo = GDF.read(path_shape)
df_geo = innerjoin(df_geo, df_vac[:, ["shapeName", "EVP"]], on="shapeName")
nothing

include("geo_ana.jl")

# Save information
create_vaccination_coverage_map(df_geo)
save_vaccination_coverage_data(df_vac) # For the publication purpose.

# ## Read and check WorldPop data

path_zaf_0_4 = "../data_pop/zaf_merged_0_4.tif"
path_zaf_whole = "../data_pop/zaf_merged_whole_pop.tif"

zaf_0_4 = read(Raster(path_zaf_0_4))
zaf_whole = read(Raster(path_zaf_whole))
nothing

print_basic_map_info(zaf_0_4)

print_basic_map_info(zaf_whole)

pl = plot( )
add_zaf_borders!(pl)
plot!(pl, zaf_0_4,
    framestyle=:none, fmt=:png,
    xlim=[15, 35], ylim=[-35, -21], border=:none,
    title="", right_margins=12Plots.mm
)

pl = plot( )
add_zaf_borders!(pl)
plot!(pl, zaf_whole,
    framestyle=:none, fmt=:png,
    xlim=[15, 35], ylim=[-35, -21], border=:none,
    title="", right_margins=12Plots.mm
)

cut_off_pop = 100
df_ras = raster_to_df(zaf_0_4)
cut_validate_raster_dataframe!(df_ras, cut_off_pop)

visualise_population(df_ras, :value, zaf_0_4, "")

bins = [2 + 0.25*i for i in 0:12]
pl = histogram(log10.(df_ras[:, :value]),
    bins=bins,
    ylabel="Number of locations",
    xlabel="log10(Population size)",
    legend=:none, dpi=300, fmt=:png
)
savefig(pl, "../res/fig_pop_histogram.png")
display(pl)

# ## Relate vaccination coverage data to population data.

include("geo_ana.jl")

# +
# df_ras: population with coordinates after filtering.
zaf_0_4 = read(Raster(path_zaf_0_4))
df_ras = raster_to_df(zaf_0_4)
cut_validate_raster_dataframe!(df_ras, cut_off_pop)

# Merge the whole population value.
zaf_whole = read(Raster(path_zaf_whole))
df_ras_whole = @pipe raster_to_df(zaf_whole) |>
    rename(_, :value => :value_whole)
df_ras = leftjoin(df_ras, df_ras_whole, on=[:lon, :lat])

df_geo # vaccination coverage by district with polygon data.
nothing
# -

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

# Save top 30 population district data (with Total population).
df_save = copy(df_mer)
df_save[!, :cum_prop] = @pipe cumsum(df_save[:, :value])/sum(df_save[:, :value]).*100 .|> round(_, digits=1)
df_save[:, :EVP] = df_save[:, :EVP] .* 100
CSV.write("../res/table_pop_top.csv", df_save)
first(df_mer, 5)

# ### Visualise the unimmunised population

include("geo_ana.jl")

pl1 = visualise_population(df_mer, :value, zaf_0_4, "")
annotate!(pl1, 17, -22, text("A", 26, :left, :black))
pl2 = visualise_population(df_mer, :value_whole, zaf_0_4, "")
annotate!(pl2, 17, -22, text("B", 26, :left, :black))
plot(pl1, pl2, dpi=300, fmt=:png, size=(1200, 400))

visualise_population(df_mer, :unvac, zaf_0_4, "")

# ### Calculate population and ES coverage by district

# +
national_ES_cov = 11.3 # old value, 8.59
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
tab = @pipe dfM |> sort(_, :prop_pc25, rev=true) |> first(_, 20)
CSV.write("../res/table_ES_pop.csv", tab)
tab
# -

# ## Calculate the probability of mobilisations

agg_scale = 230

# +
sort!(df_mer, "value", rev=true)
pop = df_mer[:, "value"]
unvac = @pipe df_mer[:, :unvac] .|> round(_, digits=0) .|> Int64
mat = df_mer[:, [:lat, :lon]] |> Matrix
π_mat = calculate_probability_of_pi(pop, mat)

spatial_p = (pop=pop, unvac=unvac, π_mat=π_mat, df=df_mer)
path = "../dt_tmp/spatial_params_agg$(agg_scale).jld2"
println(path)
jldsave(path; data=spatial_p)
# -
pop_rev = sort(pop, rev=false)
prop = pop_rev./sum(pop_rev)
cum = cumsum(prop)
scatter(log10.(pop_rev), cum, title="How population_size scenairio works",
    xlabel="Order", ylabel="Cumulative proportion (%)"
)

# ## Unvaccinated base

# +
df_mer2 = sort(df_mer, "unvac"; rev=true)
pop = df_mer2[:, "value"]
unvac = @pipe df_mer2[:, :unvac] .|> round(_, digits=0) .|> Int64
mat = df_mer2[:, [:lat, :lon]] |> Matrix
π_mat = calculate_probability_of_pi(pop, mat)

spatial_p = (pop=pop, unvac=unvac, π_mat=π_mat, df=df_mer2)
path = "../dt_tmp/spatial_params_agg$(agg_scale)_unvac.jld2"
println(path)
jldsave(path; data=spatial_p)
# -

