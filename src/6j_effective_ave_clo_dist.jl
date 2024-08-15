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
include("visualise_fig.jl")

# # Prepare the data

path_res1 = "../dt_tmp_hpc/sens_ES_catchment_20240302_050027059.jld2"
path_res2 = "../dt_tmp_hpc/sens_ES_catchment_20240302_091756993.jld2"
path_res3 = "../dt_tmp_hpc/sens_ES_catchment_20240302_190345473.jld2"
path_res1_moz = "../dt_tmp_hpc/sens_ES_catchment_20240303_020804257.jld2"
path_res2_moz = "../dt_tmp_hpc/sens_ES_catchment_20240303_064157726.jld2"
path_res3_moz = "../dt_tmp_hpc/sens_ES_catchment_20240303_160629569.jld2"

spatial_p = read_spatial_params_file("ES_population_size")
df_pop = spatial_p.df
pop = spatial_p.pop
nothing

# # Distance and early detection prob.

include("visualise_fig.jl")

p_imp_pop = pop/sum(pop) .|> Float64
imp_ws_airport = CSV.read("../data/imp_ws_airport.csv", DataFrame, header=false)[:,1]
imp_ws_moz = CSV.read("../data/imp_ws_moz.csv", DataFrame, header=false)[:,1]
nothing

# +
df50_pop = calculate_distance_and_prop_50(df_pop, p_imp_pop, path_res1)
df50_air = calculate_distance_and_prop_50(df_pop, imp_ws_airport, path_res2)
df50_moz = calculate_distance_and_prop_50(df_pop, imp_ws_moz, path_res3)

pl_w = plot(
    xlabel="Average minimum distance \nto ES covered sites (km)",
    ylabel="Early detection probability (%)", fmt=:png, 
    )
marker = (:circle,3)
plot!(pl_w, df50_pop[:, :dist], df50_pop[:, :prop_50],  label="Population size scenario", marker=marker)
plot!(pl_w, df50_air[:, :dist], df50_air[:, :prop_50],  label="Airport scenario", marker=marker)
plot!(pl_w, df50_moz[:, :dist], df50_moz[:, :prop_50],  label="Mozambique scenario", marker=marker)
display(pl_w)
# println("When the leap is observed aounrd index 5 to 6, it covers the iLembe district")
# -

# ## Outbreak potential and distance

include("utils.jl")

n_obs = 10
R0 = 14
# Calculate the effective reproduction number for each grid.
Reff = R0.* (1 .- df_pop[:, :EVP])
nothing

# +
pl = histogram(Reff, bins=0.0:0.2:4.0, ylim=[0, Inf], label="",
    xlabel="Effective reproduction number at patch i",
    ylabel="Counts", fmt=:png, dpi=200
)
p_twin = twinx(pl)

R = 0.1:0.1:4.0
p_out = observe_more_than_y_cases_in_final_dist.(10, R)
plot!(p_twin, R, p_out, color="red", ylim=[0,1],
    ylabel="Probability of an outbreak \nwith ≥10 infections",
    legend=(0.7,0.5),
    label="",#label="More than 5 cases",
)

# -

# ## Visualise the relationship between the early detection probability and weighted minimum distance to ES-covered patches

function visualise_weighted_ave_and_early_detection(
    spatial_p::NamedTuple, imp_ws_airport, imp_ws_moz;
    label1, label2, label3,
    path_res1, path_res2, path_res3,
    ylabel = "Simulated early detection probability (%)",
) 
    df = spatial_p.df
    pop = spatial_p.pop

    # Calculate the effective reproduction number for each patch.
    n_obs = 10
    R0 = 14
    Reff = R0.* (1 .- df[:, :EVP])
    p_out = observe_more_than_y_cases_in_final_dist.(n_obs, Reff)

    imp_ws_pop_out = normalize(pop .* p_out, 1)
    imp_ws_air_out = normalize(imp_ws_airport .* p_out, 1)
    imp_ws_moz_out = normalize(imp_ws_moz .* p_out, 1)
    df50_pop = calculate_distance_and_prop_50(df, imp_ws_pop_out, path_res1)
    df50_air = calculate_distance_and_prop_50(df, imp_ws_air_out, path_res2)
    df50_moz = calculate_distance_and_prop_50(df, imp_ws_moz_out, path_res3)

    yticks = [0, 20, 40, 60, 80, 100]
    pl_w_out = plot(
        xlabel="Weighted average minimum distance to \nES-covered patches (km)",
        ylabel=ylabel, 
        ylim=[-5,100], 
        yticks=(yticks, yticks),
        fmt=:png,
        tickfontsize=8,
        legendfontsize=9,
        foreground_color_legend=nothing,
        background_color_legend=nothing,
        )
    marker = (:circle,2.5)
    plot!(pl_w_out, df50_pop[:, :dist], df50_pop[:, :prop_50],  label=label1, marker=marker, lw=2)
    plot!(pl_w_out, df50_air[:, :dist], df50_air[:, :prop_50],  label=label2, marker=marker, lw=2)
    plot!(pl_w_out, df50_moz[:, :dist], df50_moz[:, :prop_50],  label=label3, marker=marker, lw=2)

    return pl_w_out, (df50_pop, df50_air, df50_moz)
end

function add_ind_site_annotations!(pl, df; inds, x_adj, y_adj, color)
    df_fil = filter(x -> x[:ind_site] in inds, df)
    for (i,r) in enumerate(eachrow(df_fil))
        annotate!(pl, 
            r[:dist] + x_adj[i], 
            r[:prop_50] .+ y_adj[i], 
            text(r[:ind_site], color, :center, 10))
    end
end

ES_pattern = "ES_population_size"
spatial_p = read_spatial_params_file(ES_pattern)
imp_ws_airport = read_imp_ws_data("airport", ES_pattern)
imp_ws_moz = read_imp_ws_data("mozambique", ES_pattern)
pl_ES_pop, dfs = visualise_weighted_ave_and_early_detection(spatial_p, imp_ws_airport, imp_ws_moz;
    label1="IMP-POP/ES-POP", label2="IMP-AIR/ES-POP", label3="IMP-LBC/ES-POP",
    path_res1=path_res1, path_res2=path_res2, path_res3=path_res3,
)
nothing

blue, orange, green = palette(:default)[1:3]
df50_pop, df50_air, df50_moz = dfs
add_ind_site_annotations!(pl_ES_pop, df50_pop; 
    inds = [1, 30, 57, 154], x_adj = [0, 0, 15, 20], y_adj = [5, 4, 3, 3],
    color=blue
)
add_ind_site_annotations!(pl_ES_pop, df50_air;
    inds  = [1, 5, 9],  x_adj = [0, 9, -8], y_adj = [5, 0, 0],
    color = orange,
)
add_ind_site_annotations!(pl_ES_pop, df50_moz;
    inds  = [1, 45, 49, 74],  x_adj = [0, 0, 0, 12], y_adj = [5, -4, -4, 0],
    color = green,
)

# +
ES_pattern = "ES_mozambique_imp_risk"
spatial_p = read_spatial_params_file(ES_pattern)
spatial_p.df[:, :EVP] ./= 100
imp_ws_airport_ES_moz = read_imp_ws_data("airport", ES_pattern)
imp_ws_moz_ES_moz = read_imp_ws_data("mozambique", ES_pattern)

pl_ES_moz, dfs = visualise_weighted_ave_and_early_detection(spatial_p, imp_ws_airport_ES_moz, imp_ws_moz_ES_moz;
    label1="IMP-POP/ES-LBC", label2="IMP-AIR/ES-LBC", label3="IMP-LBC/ES-LBC",
    path_res1=path_res1_moz, path_res2=path_res2_moz, path_res3=path_res3_moz,
    ylabel="",
)
nothing
# -

df50_pop, df50_air, df50_moz = dfs
add_ind_site_annotations!(pl_ES_moz, df50_pop; 
    inds = [1, 23, 152], x_adj = [0, 0, 0], y_adj = [4, -4, 4],
    color=blue
)
add_ind_site_annotations!(pl_ES_moz, df50_air;
    inds  = [1, 23, 152],  x_adj = [0, 20, -8], y_adj = [4, 0, 4],
    color = orange,
)
add_ind_site_annotations!(pl_ES_moz, df50_moz;
    inds  = [1, 11, 28],  x_adj = [0, 0, -10], y_adj = [-5, -4, -3],
    color = green,
)

# +
# Figure adjustment
annotate!(pl_ES_pop, (0.2, 0.85), text("A", :black, :left, :bottom, 24))
xticks = (0, 50, 100, 150, 200, 250, 300, 350)
plot!(pl_ES_pop, left_margin=5Plots.mm, bottom_margin=6Plots.mm, 
    xticks=(xticks, xticks),
)

annotate!(pl_ES_moz, (0.2, 0.85), text("B", :black, :left, :bottom, 24))

plot(pl_ES_pop, pl_ES_moz,
    fmt=:png, dpi=300, size=(1000,400)
)
# -










