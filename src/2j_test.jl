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
using Accessors
using CSV
using DataFrames
using Distributions
using Interpolations
using Parameters
using Plots
using ProtoStructs
using ProgressMeter
using Random
using StatsBase
using StatsPlots

include("util.jl")
include("model.jl")
# -

# ## Virus shedding

p = BaseParams(N0=10_000)
dump(p)

# +
virus = CSV.read("../data/Andrew2023.csv", DataFrame)
days = virus[:, :day] .|> Float64
vals = virus[:, :log_conc_shedding]

gs = relative_virus_shedding(days, vals)
#ys = ys/sum(ys
plot(days, vals, marker=:circle, label="Observed values",
     ylabel="Shedding conc., log(CID_50)",
     xticks=[i for i in 0:7:100],
)
xs = gs |> keys |> collect
ys = gs |> values |> collect
inds = sortperm(xs)

# -

# TODO: calculate hazard value analytically.
#λ0 = 76
λ0 = 1
ts = 1:p.days
I_new = vcat([1] , [0 for _ in 2:p.days])
λ = [hazard_function(I_new, gs, t; λ0=λ0)  for t in ts]
ω = 1 .- exp.(-λ)
pl = plot(
    title="1 infection at time 1", 
    ylim=(0,1.0), 
    #yticks=((1.0, 0.8)),
    yticks=([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5]),
    fmt=:png,
)
plot!(pl, ts, I_new, label="I_new", xlabel="day", ylabel="ω(t)")
plot!(pl, ts, λ, label="Virus shedding")
plot!(pl, ts, ω, label="Virus shedding")
display(pl)



# ## Simulations 

include("util.jl")
include("model.jl")

p = BaseParams()
dump(p)

# +
Random.seed!(10)
@time res = run_sim(p)
d_afp1, d_wt1 = fetch_first_date(p, res)

pl = plot(xlim=[0,150], xlabel=:days, ylabel="Count", fmt=:png)
plot!(pl, 1:p.days, res.Ia/1000, label="Ia /1000")
plot!(pl, 1:p.days, res.I_AFP, label="I_AFP")
plot!(pl, 1:p.days, res.I_new/100, label="I_new / 100")
plot!(pl, 1:p.days, res.H_new, label="H_new")
#plot!(pl, 1:p.days, res.R/1000, label="R /1000")
# -

n_sim = 1000
ds_afp1 = []
ds_wt1 = []
for i in 1:n_sim
    for _ in 1:1e4
        res = run_sim(p)
        n = res.I_new |> sum
        if n > 10; break; end
    end
    d_afp1, d_wt1 = fetch_first_date(p, res)
    ds_afp1 = vcat(ds_afp1, d_afp1)
    ds_wt1 = vcat(ds_wt1, d_wt1)
end

surv_afp1 = [0 for _ in 1:p.days]
surv_wt1 = [0 for _ in 1:p.days]
for (da, dw) in zip(ds_afp1, ds_wt1)
    surv_afp1[da] += 1
    surv_wt1[dw] += 1
end
surv_afp1 = 1 .- cumsum(surv_afp1)/n_sim
surv_wt1 = 1 .- cumsum(surv_wt1)/n_sim
pl = plot(xticks=0:30:150, fmt=:png, 
    xlim=[0,150], 
)
plot!(pl, 1:p.days, surv_afp1, label="AFP surveillance")
plot!(pl, 1:p.days, surv_wt1, label="Environmental surveillance")

lead_wt1 = ds_afp1 .- ds_wt1 
describe(lead_wt1) |> println
pl = plot(size=(200,400), fmt=:png)
df = DataFrame(Dict(:x=>lead_wt1))
@df df violin!(pl, :x)

n_delay = (lead_wt1 .< 0) |> sum 
n_delay/n_sim*100

# ## Interpolate.jl 

A_x1 = 1:.1:10
A_x2 = 1:.5:20
f(x1, x2) = log(x1+x2)
A = [f(x1,x2) for x1 in A_x1, x2 in A_x2]
itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
sitp = Interpolations.scale(itp, A_x1, A_x2)
sitp(5., 10.) # exactly log(5 + 10)
sitp(5.6, 7.1) # approximately log(5.6 + 7.1)

A = rand(20)
A_x = 1.0:2.0:40.0
nodes = (A_x,)
itp = interpolate(nodes, A, Gridded(Linear()))
x = [3, 3.5, 4, 4.5, 5.5, 6.0,  8]
plot(nodes, A, marker=:circle)
plot!(x, [itp(xx) for xx in x], marker=:circle)

A_x = 1.:2.:40.
A = [log(x) for x in A_x]
itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
sitp = Interpolations.scale(itp, A_x)
sitp(3.) # exactly log(3.)
sitp(3.5) 
plot(A_x, A, marker=:circle)
x = [3, 3.5, 4, 4.5, 5.5, 6.0,  8]
plot!(x, [sitp(xx) for xx in x], marker=:circle)


