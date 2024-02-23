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

# ## Log normal distribution

function f(x, data)
    lognorm = LogNormal(x[1], x[2])
    cdfs = cdf.(lognorm, data)
    return sum(abs.(cdfs .- q))
end

q = [0.5, 0.8]
data_x = [3.7, 13]
data_10_high = data_x./10
data_5_high = data_x./5
data_10_low = data_x.*10
data_5_low = data_x.*5


quantile(Normal(), 0.8)

get_f_res(data, inits) = optimize(x-> f(x, data), inits) |> Optim.minimizer
f_res = get_f_res(data_x, [0.0, 1.0])
f_res_10_high = get_f_res(data_10_high, [1.0, 1.0])
f_res_5_high = get_f_res(data_5_high, [1.0, 1.0])
f_res_5_low = get_f_res(data_5_low, [1.0, 1.0])
f_res_10_low = get_f_res(data_10_low, [2.4, 1.0])

println(f_res)
println(f_res_10_high)
println(f_res_5_high)
println(f_res_5_low)
println(f_res_10_low)

function draw_curve(pl, data, f_res, color, label)
    lognorm = LogNormal(f_res[1], f_res[2])
    scatter!(pl, data, q, ylim=[0,1], xlim=[0,30], color=color, label=label)
    xs = 0:0.1:30
    ys = cdf.(lognorm, xs)
    plot!(xs, ys, label="", color=color)
end

# +
pl = plot(
    xlabel="Infectious individuals per 100,000 population",
    ylabel="Probability of detection",
    fmt=:png, dpi=200,
)
cm = palette(:tab10)

draw_curve(pl, data_x, f_res, cm[1], "Baseline")
draw_curve(pl, data_10_high, f_res_10_high, cm[2], "10 times higher")
draw_curve(pl, data_5_high, f_res_5_high, cm[3], "5 times higher")
draw_curve(pl, data_5_low, f_res_5_low, cm[4], "5 times higher")
draw_curve(pl, data_10_low, f_res_10_low, cm[5], "10 times lower")
# -








