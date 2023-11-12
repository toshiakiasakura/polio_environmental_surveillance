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
using FreqTables
using LinearAlgebra
using DataFrames
using Dates
using Pipe
using Plots
using Rasters
using Serialization
using SparseArrays

include("util.jl")
include("model_meta_pop.jl")
# -

# # Transmission model 

#path = "../dt_tmp/spatial_params_agg110.ser"
path = "../dt_tmp/spatial_params_agg230.ser"
#path = "../dt_tmp/spatial_params_agg230_unvac.ser"
sp_pars = deserialize(path)
@unpack pop, unvac, π_mat = sp_pars
nothing

n_site = length(sp_pars.pop)

pars = SEIRMetaModelParams(
    R0=14.0,
    α=0.05,
    N_tot=sp_pars.pop,
    N_unvac=sp_pars.unvac,
    π_mat=sp_pars.π_mat,
    n_site=n_site,
    days=365*3,
    )
pars |> dump

w_pop_cov = mean((pop .- unvac)./pop, weights(pop))
println(w_pop_cov)
Re0 = pars.R0 * (1-w_pop_cov)

# +
#(1-w_pop_cov).* [12,13,14,15,16]
# -

# ## International airport location, index 

# +
df = sp_pars.df
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


# -

# ## Run simulations 

include("model_meta_pop.jl")

# +
pl = plot()
htmaps = []

seed_num = [48, 1, 21, 60, 65, 66] # For R0=14.
for i in 1:6
    Random.seed!(seed_num[i])
    rec, outcome, pars = run_sim(pars; rec_flag=true)
    #println(outcome.R_final_num, ", ", outcome.R_final_site)
    ht = heatmap_meta_pop(rec.I[1:100,:] |> Matrix)
    push!(htmaps, ht)
end
plot(htmaps..., size=(800,800), layout=(3,2), title="", 
    fmt=:png, dpi=300
)
# -

# ### Long run simulations 

n_sim=50

Random.seed!(48)
dump(pars)
path = @time run_and_save_sim(pars; n_sim=n_sim)

# ## Check 

paths = fetch_sim_paths(path)

# ## Baseline result for Surveillance part 

# +
# path = "../dt_tmp/20230716_211227" # n_sim=100, R0=14.0, α=0.05 , 84.87
# path = "../dt_tmp/20230716_212937" # n_sim=1000, R0=14.0, α=0.05 , 929.581375 seconds
# path = "../dt_tmp/20230716_211326" # n_sim=1000, R0=14.0, α=0.10 , 822.106013 seconds
# -

n_sim = fetch_sim_paths(path) |> length |> println
paths = fetch_sim_paths(path)
res = deserialize(paths[1])
res.pars |> dump

# +
par_AFP = AFPSurParams()
# Hazard rate
p = 0.9
g = - log(1-p)/10
cov_rate = 0.13 # Previously set to be 0.5
# Catchment area
pop = sp_pars.pop
cum_prop = cumsum(pop/sum(pop))
cov50 = abs.(cum_prop .- cov_rate) |> argmin
println("Cum index:", cov50)
sum(pop[1:cov50])/sum(pop)*100 |> println
area = fill(0, n_site)
area[1:cov50] .= 1

par_ES = ESParams(g=g, area=area)
nothing
# -

dump(par_AFP)
dump(par_ES)

include("model_meta_pop.jl")

sim_res = @time collect_summary_statistics(path, par_AFP, par_ES)
nothing

# +
df = DataFrame(sim_res)
n_sim = size(df)[1]
#h1 = histogram(df[:, :R_final_num], bins=20, title="Final size")
h2 = histogram(df[:, :R_final_site], bins=20, 
    xlabel="Number of sites having >=1 infections", ylabel="Density", 
    norm=true, legend=:none
)
h3 = histogram(df[:, :R_final_AFP], bins=20, 
    xlabel="Number of final AFP cases", 
    norm=true, legend=:none,
    ylabel="Density"
)

max_f = maximum(df[: ,:R_final_AFP])
df_fil = filter(x -> x.R_final_AFP != 0, df)
histogram!(h3, df_fil[:, :R_final_AFP], bins=20, legend=:none,
    inset = bbox(0.5Plots.w, 0.1Plots.h, 0.5Plots.w, 0.5Plots.h), 
    subplot=2,
    xlabel="Number of final AFP cases excluding 0",
    ylabel="Density",
    norm=true,
)
l = @layout [a{0.4w} b]
plot(h2, h3,  size=(900,400), fmt=:png, dpi=300,
    left_margin=5Plots.mm, 
    bottom_margin=5Plots.mm, 
    layout=l, 
)

# +
days = res.pars.days
ts_ES = df[:, "t_ES"]
ts_AFP = df[:, "t_AFP"]
t_extinct = df[:, "t_extinct"]


pl1 = plot(
    xlabel="Day", ylabel="Cumulative probability of first detection",
    legend=(0.6, 0.7),
)
annotate!((0.07, 0.95), "(A)")
cum_ES = cumulative_counts(ts_ES, days; prop=true)
cum_AFP = cumulative_counts(ts_AFP, days; prop=true)
cum_ES_cond = conditional_cumulative_prob(ts_ES, t_extinct, days)
cum_AFP_cond = conditional_cumulative_prob(ts_AFP, t_extinct, days)

plot!(pl1, 1:days, cum_ES, label="Prob. via ES", color=1)
plot!(pl1, 1:days, cum_ES_cond, label="Conditional Prob. via ES", color=1, linestyle=:dashdot)
plot!(pl1, 1:days, cum_AFP, label="Prob. via AFP surv.", color=2)
plot!(pl1, 1:days, cum_AFP_cond, label="Conditional Prob. via AFP surv.", color=2, linestyle=:dashdot)

dif = leadtime_diff(df)
x = [1 for i in 1:length(dif)]
pl2 = violin(x, dif, xticks=:none, ylabel="Lead time of ES (day)", legend=:none)
boxplot!(pl2, x, dif, fillalpha=0.75)
annotate!((0.15, 0.95), "(B)")
l = @layout [a{0.75w} b]
pl = plot(pl1, pl2, 
    fmt=:png, dpi=300, layout=l,
    size=(800,500), left_margin=5Plots.mm,
)
display(pl)
savefig(pl, "../res/fig_baseline_cum_lead.png")

# +
#first(df, 5)

outcome_num_prop(df) |> display
detect_pattern(df) |> countmap |> display

dfM = @pipe filter(x -> isnan(x["t_ES"]) == false, df)  
m_ES  = @pipe mean(dfM[:, "t_ES"]) |> round(_, digits=2)
dfM = @pipe filter(x -> isnan(x["t_AFP"]) == false, df)  
m_AFP = @pipe mean(dfM[:, "t_AFP"])  |> round(_, digits=2)
println("Mean ES: $m_ES days, Mean AFP: $m_AFP days")

#vis_cumulative_prob(df, pars.days; title="sim =$n_sim, R0=$(res.pars.R0)")
# -

dif = leadtime_diff(df)
leadtime_diff_statistics(dif)

# ## Sensitivity all

include("model_meta_pop.jl")

#path = "../dt_tmp/20230702_213018" # 100 simulations
#path = "../dt_tmp/20230702_213608" # 5000 simulations for R0=1.05
#path = "../dt_tmp/20230702_221538" # 5000 simulation for R0=1.10
n_sim = fetch_sim_paths(path) |> length

# +
#path_objs = fetch_sim_paths(path)
#res = deserialize(path_objs[1])
# -

res_all = sensitivity_ana_all(path, par_AFP, par_ES)
nothing

now_str = get_today_time()
path_save = "../dt_tmp_res/$(now_str).ser"
res_all_path = (path=path, res_all=res_all)
println(path_save)
serialize(path_save, res_all_path)










