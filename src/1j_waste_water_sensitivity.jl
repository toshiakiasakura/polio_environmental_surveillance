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

using DataFrames
using Distributions
using GLM
using Optim
using Plots

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


get_f_res(data, inits) = optimize(x-> f(x, data), inits) |> Optim.minimizer
f_res = get_f_res(data_x, [0.0, 1.0])
f_res_10_high = get_f_res(data_10_high, [1.0, 1.0])
f_res_5_high = get_f_res(data_5_high, [1.0, 1.0])
f_res_5_low = get_f_res(data_5_low, [1.0, 1.0])
f_res_10_low = get_f_res(data_10_low, [2.0, 1.0])

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

pl = plot(
    xlabel="Incidence (/100,000)",
    ylabel="Probability of detection",
)
draw_curve(pl, data_x, f_res, "blue", "Baseline")
draw_curve(pl, data_10_high, f_res_10_high, "orange", "10 times higher")
draw_curve(pl, data_5_high, f_res_5_high, "orange", "5 times higher")
draw_curve(pl, data_5_low, f_res_5_low, "green", "5 times higher")
draw_curve(pl, data_10_low, f_res_10_low, "green", "10 times lower")





# # Probit function
# Since around very low values, it does not reach 0. Discard....

function probit_func(x; a=1, b=1)
    norm = Normal()
    return cdf(norm, a + b*x)
end


data1 = DataFrame(X=[0, 4, 9]./100_000, Y=[0, 0.5, 0.8])
data2 = DataFrame(X=[1.6, 3.6]./100_000, Y=[0.5, 0.8])
data3 = DataFrame(X=[40, 90]./100_000, Y=[0.5, 0.8])

probit = glm(@formula(Y ~ X), data1, Binomial(), ProbitLink())
a1, b1 = round.(coef(probit); digits=3)
probit = glm(@formula(Y ~ X), data2, Binomial(), ProbitLink())
a2, b2 = round.(coef(probit); digits=3)
probit = glm(@formula(Y ~ X), data3, Binomial(), ProbitLink())
a3, b3 = round.(coef(probit); digits=3)

display([a1, b1])
display([a2, b2])
display([a3, b3])

x = (0.1:0.1:100) ./ 100_000
y1 = probit_func.(x; a=a1, b=b1)
y2 = probit_func.(x; a=a2, b=b2)
y3 = probit_func.(x; a=a3, b=b3)
nothing

scatter(data1.X, data1.Y, ylim=[0,1])
scatter!(data2.X, data2.Y)
scatter!(data3.X, data3.Y)
plot!(x,y1)
plot!(x,y2)
plot!(x,y3)




