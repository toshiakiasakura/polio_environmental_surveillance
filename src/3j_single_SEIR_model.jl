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
using CSV
using DataFrames
using Distributions
using ShiftedArrays
using StatsPlots
using Optim

include("model.jl")
# -

path = "../dt_tmp/spatial_params_agg110.ser"
spatial_p = deserialize(path)
#spatial_p = (pop=pop, unvac=unvac, π_mat=π_mat, df=df_mer)
@unpack pop, unvac, π_mat = spatial_p
nothing

w_pop_cov = mean((pop .- unvac)./pop, weights(pop))
println(w_pop_cov)
N_tot = 100_000
N_unvac = round(N_tot * (1- w_pop_cov) )
g = - log(1-0.9)/10

pars = SEIRModelParams(
    N_tot=N_tot, N_unvac=N_unvac, 
    days=365*3, g=g, 
    R0=14
)
dump(pars)

Re0 = pars.R0 * (1-w_pop_cov)

# ## Check AFP incubation time 
# Note: TO check this model change A as [1000, 0, 0, 0, 0, 0] within `initialize_model` function.

params_inc = SEIRModelParams(
    N_tot=N_tot, N_unvac=N_unvac, 
    days=365*3, g=g, 
    R0=0, I0_init=1
)

include("model.jl")

pl = plot(xlim=[0,40], ylim=[0, 0.1], xticks=[0, 7, 15, 21, 25, 33, 39], fmt=:png,
    xlabel="day", ylabel="Prob",
)
for _ in 1:100
    rec, outcome = run_sim(params_inc; rec_flag=true)
    plot!(pl, 1:rec.days, rec.Z_A5_6/1000, label="New AFP", legend=:none)
end
display(pl)

# # Check single run

Random.seed!(5)
rec, outcome = run_sim(pars; rec_flag=true)
println(outcome)
pl = plot(xlim=[0,365*3], ylim=[0, 200])
plot!(pl, 1:rec.days, rec.E, label="E")
plot!(pl, 1:rec.days, rec.I, label="I")
plot!(pl, 1:rec.days, rec.Z_A5_6*100, label="New AFP")
plot!(pl, 1:rec.days, rec.R, label="R")

outcomes = []
N_sim = 10000
for i in 1:N_sim
    rec, outcome = run_sim(pars; rec_flag=false)
    push!(outcomes, outcome)
end
df = DataFrame(outcomes)
first(df, 5)

outcome_num_prop(df)

# +
days = pars.days
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
savefig(pl, "../res/fig_single_prob_lead.png")
# -

println("R0 = $(pars.R0), Re0 = $((1-w_pop_cov)*pars.R0)")

dif = leadtime_diff(df)
leadtime_diff_statistics(dif)

# +
println("Need of inclusion of other vaccinated individuals to the model.")
max_f = maximum(df[: ,:R_final_AFP])
pl2 = histogram(df[:, :R_final_AFP], bins=0:max_f, legend=:none, 
    xlabel="Final AFP cases",
    ylabel="Probability",
    norm=true,
    fmt=:png, dpi=300
    
)
df_fil = filter(x -> x.R_final_AFP != 0, df)
histogram!(pl2, df_fil[:, :R_final_AFP], bins=1:max_f, legend=:none,
    inset = bbox(0.4Plots.w, 0.1Plots.h, 0.5Plots.w, 0.5Plots.h), 
    subplot=2,
    xlabel="Final AFP cases excluding 0",
    ylabel="Probability",
    norm=true,
)
# -




