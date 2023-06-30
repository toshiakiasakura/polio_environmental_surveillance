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

params = SEIRModelParams(N0=10_000, days=360, λ0=0.2302, R0=2)
dump(params)

rec, outcome = run_sim(params; rec_flag=true)
println(outcome)
pl = plot(xlim=[0,360], ylim=[0,10000])
plot!(pl, 1:rec.days, rec.E, label="E")
plot!(pl, 1:rec.days, rec.Ia, label="Ia")
plot!(pl, 1:rec.days, rec.I_AFP*100, label="I_AFP per")
plot!(pl, 1:rec.days, rec.R, label="R")

outcomes = []
N_sim = 1000
for i in 1:N_sim
    rec, outcome = run_sim(params; rec_flag=false)
    push!(outcomes, outcome)
end
df = DataFrame(outcomes)
first(df, 5)

outcome_num_prop(df)

detect_pattern(df)

# +
days = params.days
ts_ES = df[:, "t_ES"]
ts_AFP = df[:, "t_AFP"]
t_extinct = df[:, "t_extinct"]
cum_ES = cumulative_counts(ts_ES, params.days; prop=true)
cum_AFP = cumulative_counts(ts_AFP, params.days; prop=true)
plot(1:params.days, cum_ES, label="ES", fmt=:png)
plot!(1:params.days, cum_AFP, label="AFP")

cum_ES = conditional_cumulative_prob(ts_ES, t_extinct, params.days)
cum_AFP = conditional_cumulative_prob(ts_AFP, t_extinct, params.days)
plot!(1:params.days, cum_ES, label="Pc, ES")
plot!(1:params.days, cum_AFP, label="Pc, AFP")


# -

diff = leadtime_diff(df)
leadtime_diff_statistics(diff)

println("Need of inclusion of other vaccinated individuals to the model.")
histogram(df[:, :R_final], bins=100)




