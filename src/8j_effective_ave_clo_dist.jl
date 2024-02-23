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

path_res1 = "../dt_tmp_res/sens_ES_catchment_20240219_145929.jld2"
path_res2 = "../dt_tmp_res/sens_ES_catchment_20240219_150715.jld2"
path_res3 = "../dt_tmp_res/sens_ES_catchment_20240219_152555.jld2"

path_res1_moz = "../dt_tmp_res/sens_ES_catchment_20240217_123152.jld2"
path_res2_moz = "../dt_tmp_res/sens_ES_catchment_20240217_124140.jld2"
path_res3_moz = "../dt_tmp_res/sens_ES_catchment_20240217_130352.jld2"

df_zaf = CSV.read("../res/table_pop_top.csv", DataFrame)
nothing

# # Distance and early detection prob.

include("visualise_fig.jl")

pop = df_zaf[:, :value]
p_imp_pop = pop/sum(pop)
imp_ws_airport = CSV.read("../data/imp_ws_airport.csv", DataFrame, header=false)[:,1]
imp_ws_moz = CSV.read("../data/imp_ws_moz.csv", DataFrame, header=false)[:,1]
nothing

# +
d_sens_pop, prop_50_pop = calculate_distance_and_prop_50(df_zaf, p_imp_pop, path_res1)
d_sens_air, prop_50_air = calculate_distance_and_prop_50(df_zaf, imp_ws_airport, path_res2)
d_sens_moz, prop_50_moz = calculate_distance_and_prop_50(df_zaf, imp_ws_moz, path_res3)

pl_w = plot(
    xlabel="Average minimum distance \nto ES covered sites (km)",
    ylabel="Early detection probability (%)", fmt=:png,
    )
marker = (:circle,3)
plot!(pl_w, d_sens_pop, prop_50_pop[:, :prop_50],  label="Population size scenario", marker=marker)
plot!(pl_w, d_sens_air, prop_50_air[:, :prop_50],  label="Airport scenario", marker=marker)
plot!(pl_w, d_sens_moz, prop_50_moz[:, :prop_50],  label="Mozambique scenario", marker=marker)
#display(pl_w)
println("When the leap is observed aounrd index 5 to 6, it covers the iLembe district")
# -

# ## Outbreak potential and distance

include("utils.jl")

n_obs = 10
R0 = 14
# Calculate the effective reproduction number for each grid.
Reff = R0.* (1 .- df_zaf[:, :EVP]./100)
nothing

# +
pl = histogram(Reff, bins=0.0:0.2:4.0, ylim=[0, Inf], label="",
    xlabel="Reproduction number",
    ylabel="Counts", fmt=:png,
)
p_twin = twinx(pl)

R = 0.1:0.1:4.0
p_out = observe_more_than_y_cases_in_final_dist.(10, R)
plot!(p_twin, R, p_out, color="red", ylim=[0,1],
    ylabel="Probability of an outbreak \nwith 10 or more infections",
    legend=(0.7,0.5),
    label="",#label="More than 5 cases",
)

# -

# Set the observe final distribution.
p_out = observe_more_than_y_cases_in_final_dist.(n_obs, Reff)
nothing

# For population scenario
p_imp_pop_out = pop.*p_out/sum(pop.*p_out)
p_imp_air = imp_ws_airport.*p_out/sum(imp_ws_airport.*p_out)
p_imp_moz = imp_ws_moz.*p_out/sum(imp_ws_moz.*p_out)
nothing

# +
d_sens_pop, prop_50_pop = calculate_distance_and_prop_50(df_zaf, p_imp_pop_out, path_res1)
d_sens_air, prop_50_air = calculate_distance_and_prop_50(df_zaf, p_imp_air, path_res2)
d_sens_moz, prop_50_moz = calculate_distance_and_prop_50(df_zaf, p_imp_moz, path_res3)

pl_w_out = plot(
    xlabel="Average minimum distance to \nES covered sites (km)",
    ylabel="Early detection prob. (%)", fmt=:png,
    )
marker = (:circle,3)
plot!(pl_w_out, d_sens_pop, prop_50_pop[:, :prop_50],  label="Population size scenario", marker=marker)
plot!(pl_w_out, d_sens_air, prop_50_air[:, :prop_50],  label="Airport scenario", marker=marker)
plot!(pl_w_out, d_sens_moz, prop_50_moz[:, :prop_50],  label="Mozambique scenario", marker=marker)
println("When the leap is observed aounrd index 5 to 6, it covers the iLembe district")

# Figure adjustment
annotate!(pl_w, 30, 90, text("A", :black, :left, 24))
plot!(pl_w, left_margin=5Plots.mm, bottom_margin=6Plots.mm)

annotate!(pl_w_out, 30, 90, text("B", :black, :left, 24))

plot(pl_w, pl_w_out,
    fmt=:png, dpi=300, size=(1000,400)
)
# -






