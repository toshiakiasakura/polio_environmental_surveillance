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
using Optim
using Pipe
using Plots
using Rasters
using Shapefile

include("utils.jl")
include("geo_ana.jl")
include("model_meta_pop.jl")
include("visualise_fig.jl")
# -

df_zaf = CSV.read("../res/table_pop_top.csv", DataFrame)
nothing

# # Calculate the d,ave

pop = df_zaf[:, :value]
mat = df_zaf[:, [:lat, :lon]] |> Matrix
sens_index = obtain_ES_sensitivity_index(pop, 0.01)
p_imp = pop/sum(pop)
d_sens = calculate_d_ave_over_index(df_zaf, p_imp)
nothing

pl = plot(sens_index, d_sens, xlabel="# of ES covered sites", ylabel="d, ave", fmt=:png,
    title="Relationship between # of ES covered sites and distance",
)
scatter!(pl, sens_index, d_sens)

# +
path_res1 = "../dt_tmp_res/20231112_005011.ser" # 5000 samples
prop_50 = fetch_early_det_50(path_res1)

plot(sens_index, prop_50, xlim=[0,100])
# -

# # Distance and early detection prob.

include("visualise_fig.jl")

path_res1 = "../dt_tmp_res/20231112_005011.ser" # 5000 samples
path_res2 = "../dt_tmp_res/20231114_211615.ser" # Importation is radiation
path_res3 = "../dt_tmp_res/20231112_043848.ser" # 5000 samples

path_res1_moz = "../dt_tmp_res/20231117_111513.ser"
path_res2_moz = "../dt_tmp_res/20231117_125438.ser"
path_res3_moz = "../dt_tmp_res/20231117_162723.ser"

path_spatial = "../dt_tmp/spatial_params_agg230.ser"

# +
sp_pars = deserialize(path_spatial)
per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
sens_index = obtain_ES_sensitivity_index(sp_pars.pop, 0.01)
x_per_pop =  per_pop[sens_index]

pop = df_zaf[:, :value]
nothing
# -

p_imp_pop = pop/sum(pop)
imp_ws_airport = CSV.read("../data/imp_ws_airport.csv", DataFrame, header=false)[:,1]
imp_ws_moz = CSV.read("../data/imp_ws_moz.csv", DataFrame, header=false)[:,1]
nothing

# +
# For old airport pattern.
#n_sites = pop |> length
#ind_imp, inbound = [11, 7, 62], [4_342_611, 1_156_996, 188_243]
#inbound = inbound/sum(inbound)
#p_imp_air = fill(0.0, n_sites)
#for (ind, p) in zip(ind_imp, inbound)
#    p_imp_air[ind] = p
#end
# -

d_sens_pop, prop_50_pop = calculate_distance_and_prop_50(df_zaf, p_imp_pop, path_res1)
d_sens_air, prop_50_air = calculate_distance_and_prop_50(df_zaf, imp_ws_airport, path_res2)
d_sens_moz, prop_50_moz = calculate_distance_and_prop_50(df_zaf, imp_ws_moz, path_res3)
nothing

pl = plot(
    xlabel="Average minimum distance to ES covered sites (km)",
    ylabel="Early detection prob. (%)", fmt=:png,
    )
marker = (:circle,3)
plot!(pl, d_sens_pop, prop_50_pop,  label="Population scenario", marker=marker)
plot!(pl, d_sens_air, prop_50_air,  label="Airport scenario", marker=marker)
plot!(pl, d_sens_moz, prop_50_moz,  label="Mozambique scenario", marker=marker)
display(pl)
println("When the leap is observed aounrd index 5 to 6, it covers the iLembe district")

# ## Outbreak potential and distance

include("utils.jl")

n_obs = 10

R0 = 14
Reff = R0.* (1 .- df_zaf[:, :EVP]./100)
nothing

# +
pl = histogram(Reff, bins=0.0:0.2:4.0, ylim=[0, Inf], label="",
    xlabel="Reproduction number",
    ylabel="Counts", fmt=:png,
)
p_twin = twinx(pl)

R = 0.1:0.1:4.0
p_out = observe_more_than_y_cases_in_final_dist.(5, R)
plot!(p_twin, R, p_out, label="More than 5 cases", color="red", ylim=[0,1],
    ylabel="Probability", legend=(0.7,0.5)
)
p_out = observe_more_than_y_cases_in_final_dist.(10, R)
plot!(p_twin, R, p_out, label="More than 10 cases", color="green", ylim=[0,1])
p_out = observe_more_than_y_cases_in_final_dist.(30, R)
plot!(p_twin, R, p_out, label="More than 30 cases", color="yellow", ylim=[0,1])
# -

# Set the observe final distribution.
p_out = observe_more_than_y_cases_in_final_dist.(n_obs, Reff)
nothing

# For population scenario
p_imp_pop_out = pop.*p_out/sum(pop.*p_out)
p_imp_air = imp_ws_airport.*p_out/sum(imp_ws_airport.*p_out)
p_imp_moz = imp_ws_moz.*p_out/sum(imp_ws_moz.*p_out)
nothing

d_sens_pop, prop_50_pop = calculate_distance_and_prop_50(df_zaf, p_imp_pop_out, path_res1)
d_sens_air, prop_50_air = calculate_distance_and_prop_50(df_zaf, p_imp_air, path_res2)
d_sens_moz, prop_50_moz = calculate_distance_and_prop_50(df_zaf, p_imp_moz, path_res3)

pl = plot(
    xlabel="Average minimum distance to ES covered sites (km)",
    ylabel="Early detection prob. (%)", fmt=:png,
    )
marker = (:circle,3)
plot!(pl, d_sens_pop, prop_50_pop,  label="Population scenario", marker=marker)
plot!(pl, d_sens_air, prop_50_air,  label="Airport scenario", marker=marker)
plot!(pl, d_sens_moz, prop_50_moz,  label="Mozambique scenario", marker=marker)
display(pl)
println("When the leap is observed aounrd index 5 to 6, it covers the iLembe district")

# ## This figure is difficult to be interpreted!

# +
pop = df_zaf[:, :value]
per_pop = cumsum(pop)/sum(pop)*100
sens_index = obtain_ES_sensitivity_index(pop, 0.01)
per_pop_ind = per_pop[sens_index]
pl = plot(
    xlabel="ES covered population",
    ylabel="100 - Early detection prob (%)"
)
p_twin = twinx(pl,
)
plot!(pl, per_pop_ind, 100 .- prop_50_pop, label="", marker=:circle)
plot!(p_twin, per_pop_ind, d_sens_pop, label="Population scenario", linestyle=:dash,
    ylabel="Average minimum distance"
)

plot!(pl, per_pop_ind, 100 .- prop_50_air, label="", marker=:circle)
plot!(p_twin, per_pop_ind, d_sens_air, label="Airport scenario", linestyle=:dash)

plot!(pl, per_pop_ind, 100 .- prop_50_moz, label="", marker=:circle)
plot!(p_twin, per_pop_ind, d_sens_moz, label="Mozambique scenario", linestyle=:dash)
# -
# ## Mozambique scenairo

path_spatial_sorted = "../dt_tmp/spatial_params_agg230_moz_sorted.ser"
sp_pars = deserialize(path_spatial_sorted)
df_zaf_moz = sp_pars.df
pop = sp_pars.pop
nothing

p_imp_pop = pop/sum(pop)
imp_ws_airport = CSV.read("../data/imp_ws_airport_sorted.csv", DataFrame, header=false)[:,1]
imp_ws_moz = CSV.read("../data/imp_ws_moz_sorted.csv", DataFrame, header=false)[:,1]
nothing

R0 = 14
Reff = R0.* (1 .- df_zaf_moz[:, :EVP]./100)
p_out = observe_more_than_y_cases_in_final_dist.(n_obs, Reff)
nothing

# For population scenario
p_imp_pop_out = pop.*p_out/sum(pop.*p_out)
p_imp_air = imp_ws_airport.*p_out/sum(p_imp_air.*p_out)
p_imp_moz = imp_ws_moz.*p_out/sum(imp_ws_moz.*p_out)
nothing

d_sens_pop, prop_50_pop = calculate_distance_and_prop_50(df_zaf_moz, p_imp_pop_out, path_res1_moz)
d_sens_air, prop_50_air = calculate_distance_and_prop_50(df_zaf_moz, p_imp_air, path_res2_moz)
d_sens_moz, prop_50_moz = calculate_distance_and_prop_50(df_zaf_moz, p_imp_moz, path_res3_moz)
nothing

pl = plot(
    xlabel="Average minimum distance to ES covered sites (km)",
    ylabel="Early detection prob. (%)", fmt=:png,
    )
marker = (:circle,3)
plot!(pl, d_sens_pop, prop_50_pop,  label="Population scenario", marker=marker)
plot!(pl, d_sens_air, prop_50_air,  label="Airport scenario", marker=marker)
plot!(pl, d_sens_moz, prop_50_moz,  label="Mozambique scenario", marker=marker)
display(pl)
println("When the leap is observed aounrd index 5 to 6, it covers the iLembe district")






# ### Borel-Tanner distribution

x = 1:20
pl = plot()
R = 0.2
y = borel_tanner_dist.(x, R)
plot!(pl, x, y, label="R=$(R)")
R = 0.5
y = borel_tanner_dist.(x, R)
plot!(pl, x, y, label="R=$(R)")

x = 1:100_000
borel_tanner_dist.(x, 0.2) |> sum |> println
borel_tanner_dist.(x, 0.8) |> sum |> println
borel_tanner_dist.(x, 0.9) |> sum |> println
borel_tanner_dist.(x, 1.0) |> sum |> println
borel_tanner_dist.(x, 2.0) |> sum |> println

# Extinction probability for the offspring distribution
R = 10.0
function f(q, R)
    #if q[1] < 0; return Inf; end
    return abs(exp(-R*(1-q[1])) - q[1])
end
qs = []
Rs = [1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0]
for R in Rs
    res = optimize(x-> f(x, R), [0.4], LBFGS())
    push!(qs, res.minimizer[1])
end

println("Extinction probability is matched with the sum of the Borel-Tanner dist: ",qs[3])





# ## Calculate probability of pi with two districts

# +
df_zaf[:, :cnt] .= "zaf"
df_moz[:, :cnt] .= "moz"
sort!(df_zaf, "value", rev=true) # to ensure the order is the same as other files.
sort!(df_moz, "value", rev=true)

df_mer = vcat(df_zaf, df_moz)
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

CSV.write("../data/imp_ws.csv",  Tables.table(imp_ws), writeheader=false)

imp_ws_read = CSV.read("../data/imp_ws.csv", DataFrame, header=false)[:,1]
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











