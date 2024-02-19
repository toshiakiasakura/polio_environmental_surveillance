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

# ## Prepare multiple figures with ES coverage sites for gif

# +
# Prepare data
cum_per = df[:, :cum_per]
pop = df[:, :value]
sens_index = obtain_ES_sensitivity_index(pop, 0.01)
zaf_map = copy(zaf_0_4)
zaf_map[:] .= 0

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
# -

sum(sens_index .<= 33)

# These index are for gif creation process in python.
# Use the maximum index for the inclusion criteria.
national_ES_cov = 8.59
n = (df[:, :cum_per] .<= national_ES_cov/0.25) |> sum
println("pc 25, maximum index: ", sum(sens_index .<= n), ", Number of sites: ", n)
n = (df[:, :cum_per] .<= national_ES_cov/0.30) |> sum
println("pc 30, maximum index: ", sum(sens_index .<= n), ", Number of sites: ", n)
n = (df[:, :cum_per] .<= national_ES_cov/0.5) |> sum
println("pc 50, maximum index: ", sum(sens_index .<= n), ", Number of sites: ", n)

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
dist = df[ind, "shapeName"]
pl2 = visualise_probs(df, ind, "1st populous location\n$(dist)")
annotate!(pl2, (0.1, 0.95), "(B)")

ind = 2
pl3 = histogram(log10.(π_mat[ind, :]),
    legend=:none,
    xlabel = "log10(π2i)", ylabel="Number of locations",
)
annotate!(pl3, (0.1, 0.95), "(C)")
dist = df[ind, "shapeName"]
pl4 = visualise_probs(df, ind, "2nd populous location\n$(dist)")
annotate!(pl4, (0.1, 0.95), "(D)")

ind = 3
pl5 = histogram(log10.(π_mat[ind, :]),
    legend=:none,
    xlabel = "log10(π3i)", ylabel="Number of locations",
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






