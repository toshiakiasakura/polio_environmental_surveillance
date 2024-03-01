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

path = "../dt_tmp/spatial_params_agg230.jld2"
spatial_p = load(path)["data"]
spatial_p |> keys |> println
df = spatial_p.df
π_mat = spatial_p.π_mat
nothing

path_zaf_0_4 = "../data_pop/zaf_merged_0_4.tif"
zaf_0_4 = read(Raster(path_zaf_0_4))
nothing

# # Probability of importation for each scenario

function visualise_importation_prob(df_mer; 
        title::String = "", colorbar = true,
        colorbar_title = "",
    )
    zaf_map = copy(zaf_0_4)
    zaf_map[:] .= 0
    n = size(df_mer)[1]
    for i in 1:n
        x = df_mer[i, :lon]
        y = df_mer[i, :lat]
        zaf_map[ At(x), At(y)] = log10.(df_mer[i, :imp_prob])
    end
    pl = plot()
    plot!(zaf_map,
        xlim=[15,35], ylim=[-35, -21],
        axis=nothing, border=:none,
        logscale=:log10,
        #colorbar_title="log10(πij)",
        title=title,
        right_margin=5Plots.mm,
        top_margin = 0Plots.mm,
        colorbar = colorbar,
        colorbar_title = colorbar_title,
        color = :matter,
        clims = (-6.5, -0.5),
        fmt=:png, dpi=200,
    )
    add_zaf_borders!(pl)
    return pl
end

sp_pars = read_spatial_params_file("ES_population_size")
imp_pop = normalize(sp_pars.pop, 1)
imp_ws_moz = CSV.read("../data/imp_ws_moz.csv", DataFrame, header=false)[:,1]
imp_ws_airport = CSV.read("../data/imp_ws_airport.csv", DataFrame, header=false)[:,1]
nothing

log10.(imp_pop) |> describe
log10.(imp_ws_airport) |> describe
log10.(imp_ws_moz) |> describe

# +
add_annotation!(pl, text_) =  annotate!(pl, 16, -22, text(text_, :black, :left, 18))

df[:, :imp_prob] = imp_pop
pl_pop = visualise_importation_prob(df; colorbar=false)
add_annotation!(pl_pop, "A")

df[:, :imp_prob] = imp_ws_airport
pl_airport = visualise_importation_prob(df; colorbar=false)
add_annotation!(pl_airport, "B")

df[:, :imp_prob] = imp_ws_moz
pl_moz = visualise_importation_prob(df; colorbar=false)
add_annotation!(pl_moz, "C")
layout = @layout [a b; c c]
plot(pl_pop, pl_airport, pl_moz, layout=layout, 
    fmt=:png, dpi=300, size=(1200, 800))
# -

df[:, :imp_prob] = imp_pop
pl_pop = visualise_importation_prob(df; colorbar=true, colorbar_title="log10(Probability of importation)")
plot!(pl_pop, size=(800,600))

# ## For figure 1 

df[:, :imp_prob] = imp_pop
pl_pop = visualise_importation_prob(df; colorbar=false)
plot!(pl_pop, right_margin=-25Plots.mm)
nothing

# +
# International airport information
ind = [10, 7, 61]
travel_volume = [4_342_611, 1_156_996, 188_243]
#weight = normalize((travel_volume).^(1/2) , 1) # Travel volume
weight = [0.7, 0.5, 0.3]
cm = palette(:tab10)

zaf_map = copy(zaf_0_4)
zaf_map[:] .= 0
pl_airport = plot()
plot!(pl_airport, zaf_map,
    xlim=[15,35], ylim=[-35, -21],
    axis=nothing, border=:none,
    logscale=:log10,
    #colorbar_title="log10(πij)",
    title="",
    right_margin=-10Plots.mm,
    top_margin = 0Plots.mm,
    colorbar = false,
    color = :matter,
    clims = (-6.5, -0.5),
    fmt=:png, dpi=200,
)
add_zaf_borders!(pl_airport)
# Plot the internatinoal airport positions
scatter!(pl_airport, df[ind, :lon], df[ind, :lat], 
    label=false, color=cm[1:3], 
    markershape=:circle,
    markersize=weight*50, alpha=1,
)

nothing
# -

include("geo_ana.jl")

pl_moz = plot()
plot!(pl_moz, zaf_map,
    axis=nothing, border=:none,
    title="",
    xlim=[15,45], ylim=[-35, -10],
)
add_zaf_ADM0_borders!(pl_moz)
add_moz_borders!(pl_moz)
scatter!(pl_moz, 
    [24, 33, 40, 33], 
    [-30, -15, -13, -22], 
    label=false, color=cm[4],
    markershape=:rect,
    markersize=[12, 8, 8, 8],
)
nothing

plot(pl_pop, pl_airport, pl_moz, 
    fmt=:png, size=(1800, 600), dpi=150,
    layout=@layout[a b c],
)







# ## Prepare ES layout for figure 1

function plot_ES_covered_sites(
    ES_pattern, r_ind
)
    spatial_p = read_spatial_params_file(ES_pattern)
    df = spatial_p.df
    
    # Prepare data
    zaf_map = copy(zaf_0_4)
    zaf_map[:] .= 0
    
    ind = 1:r_ind
    for r in eachrow(df[ind,:])
        zaf_map[ At(r.lon), At(r.lat)] = 100
    end
    pl = plot(zaf_map,
        xlim=[15,35], ylim=[-35, -21],
        color=:blue,
        colorbar=false,
        title="",
        axis=nothing, border=:none,
        dpi=100, fmt=:png, size=(1200,1200),
    )
    add_zaf_borders!(pl)
    pl
end

pl_pop = plot_ES_covered_sites("ES_population_size", 50)
plot!(pl_pop, right_margin=-20Plots.mm)
pl_moz = plot_ES_covered_sites("ES_mozambique_imp_risk", 50)
plot(pl_pop, pl_moz, fmt=:png, dpi=150, size=(1200,600))

# ## Prepare multiple figures with ES coverage sites for gif

# Prepare data
cum_per = df[:, :cum_per]
pop = df[:, :value]
sens_index = obtain_ES_sensitivity_index(pop, 0.01)
zaf_map = copy(zaf_0_4)
zaf_map[:] .= 0
nothing

# Create multiple files.
l_ind =  1
for (i, r_ind) in enumerate(sens_index)
    ind = l_ind:r_ind
    for r in eachrow(df[ind,:])
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

# +
# These index are for gif creation process in python.
# Use the maximum index for the inclusion criteria.

# sum(sens_index .<= n) is used to specify r_ind. 
national_ES_cov = 11.3 # for old value, 8.59
n = (df[:, :cum_per] .<= national_ES_cov/0.25) |> sum
println("pc 25, maximum index: ", sum(sens_index .<= n), ", Number of sites: ", n)
n = (df[:, :cum_per] .<= national_ES_cov/0.30) |> sum
println("pc 30, maximum index: ", sum(sens_index .<= n), ", Number of sites: ", n)
n = (df[:, :cum_per] .<= national_ES_cov/0.5) |> sum
println("pc 50, maximum index: ", sum(sens_index .<= n), ", Number of sites: ", n)
# -

# ## Visualise the propability of mobilisation from 3 top populous areas

include("utils.jl")

function visualise_probs(df_mer, ind::Int64, title::String)
    zaf_map = copy(zaf_0_4)
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
        color=:matter,
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
annotate!(pl1, (0.1, 0.95), "A")
dist = df[ind, "shapeName"]
pl2 = visualise_probs(df, ind, "1st populous location\n$(dist)")
annotate!(pl2, (0.1, 0.95), "B")

ind = 2
pl3 = histogram(log10.(π_mat[ind, :]),
    legend=:none,
    xlabel = "log10(π2i)", ylabel="Number of locations",
)
annotate!(pl3, (0.1, 0.95), "C")
dist = df[ind, "shapeName"]
pl4 = visualise_probs(df, ind, "2nd populous location\n$(dist)")
annotate!(pl4, (0.1, 0.95), "D")

ind = 3
pl5 = histogram(log10.(π_mat[ind, :]),
    legend=:none,
    xlabel = "log10(π3i)", ylabel="Number of locations",
)
annotate!(pl5, (0.1, 0.95), "E")
dist = df[ind, "shapeName"]
pl6 = visualise_probs(df, ind, "3rd populous location\n$(dist)")
annotate!(pl6, (0.1, 0.95), "F")

l = @layout [
    a{0.3w} b
    c{0.3w} d
    e{0.3w} f
]
pl = plot(pl1, pl2, pl3, pl4, pl5, pl6,
    dpi=300, size=(800, 350 * 3), fmt=:png,
    layout=l, left_margin=5Plots.mm,
)
display(pl)
savefig(pl, "../res/fig_pi_map.png")
# -
# ## Unvaccinated base

path = "../dt_tmp/spatial_params_agg230_unvac.jld2"
spatial_p |> keys |> println
df = spatial_p.df
π_mat = spatial_p.π_mat
nothing

# +
ind = 1
pl1 = histogram(log10.(π_mat[:, ind]),
    legend=:none,
    xlabel = "log10(πi1)", ylabel="Number of locations",
)
annotate!(pl1, (0.1, 0.95), "(A)")
dist = df[ind, "shapeName"]
pl2 = visualise_probs(df, ind, "1st populous location\n$(dist)")
annotate!(pl2, (0.1, 0.95), "(B)")

ind = 2
pl3 = histogram(log10.(π_mat[:, ind]),
    legend=:none,
    xlabel = "log10(πi2)", ylabel="Number of locations",
)
annotate!(pl3, (0.1, 0.95), "(C)")
dist = df[ind, "shapeName"]
pl4 = visualise_probs(df, ind, "2nd populous location\n$(dist)")
annotate!(pl4, (0.1, 0.95), "(D)")

ind = 3
pl5 = histogram(log10.(π_mat[:, ind]),
    legend=:none,
    xlabel = "log10(πi3)", ylabel="Number of locations",
)
annotate!(pl5, (0.1, 0.95), "(E)")
dist = df[ind, "shapeName"]
pl6 = visualise_probs(df, ind, "3rd populous location\n$(dist)")
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






