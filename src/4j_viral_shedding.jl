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

using CSV
using DataFrames
using Distributions
using ShiftedArrays
using StatsPlots
using Optim
#using Turing

df = CSV.read("../data/Radboud_2013_Pt.csv", DataFrame)
plot(df[:, "day"], df[:, "Pt"], marker=:circle)

lis = append!([2], 8:16)
df_pt = df[lis,:]
scale = df_pt[:, :Pt] .* 7
df_pt[:,"pdf"] = df_pt[:, :Pt]/sum(scale)
df_pt[!, :day] = df_pt[:, :day] .|> Float64
plot(df_pt[:, :day], df_pt[:,:pdf], title="pdf of proportion of shedding virus over time")

pt = df_pt[:, :pdf]
days = df_pt[:, :day]
function KLd(x)
    α = x[1]
    θ = x[2]
    qt = Distributions.pdf.(Gamma(α, θ), days)
    d = pt .* (log.(pt) .- log.(qt))
    return sum(d)
end

res = optimize(KLd, [0.0,0.0], [Inf, Inf], [1.0, 10.0])

α, θ = Optim.minimizer(res)
println("$α, $θ")
x = 0:63
y = pdf.(Gamma(α, θ), x) 
pl = plot(title="pdf of proportion of shedding virus over time")
plot!(pl, days, pt, label="Observed", marker=:circle)
plot!(pl, x, y, label="Fitted")
display(pl)

# +
pt = df_pt[:, :pdf]
days = df_pt[:, :day]
function conv(t, γ2; γ1 = 1/4.0)
    return γ1*γ2/(γ1 - γ2)*(exp(-γ2*t) - exp(-γ1*t))
end
function target(x)
    γ2 = x[1]
    qt = conv.(days, γ2)
    d = pt .* (log.(pt) .- log.(qt))
    return sum(d)
end

res = optimize(target, [0], [1/4.0], [1/8])

# +

γ2 = Optim.minimizer(res)
println("$γ2 , $(1/γ2)")
# -

x = 0:63
y_conv = conv.(x, γ2)
y_gamma = pdf.(Gamma(α, θ), x) 
pl = plot(title="pdf of proportion of shedding virus over time\nγ1=1/4, γ2=1/15.02",
    xlabel="day", ylabel="PDF",
    fmt=:png,
)
plot!(pl, days, pt, label="Observed", marker=:circle)
plot!(pl, x, y_gamma, label="Gamma fitted")
plot!(pl, x, y_conv, label="SEIR model")
display(pl)

# # Hazard rate

p = 0.9
λ0 = - log(1-p)/10
println(λ0)
I = 1:100
ωt = 1 .- exp.(-λ0.*I)
plot(I, ωt, ylim=[0,1.0])




