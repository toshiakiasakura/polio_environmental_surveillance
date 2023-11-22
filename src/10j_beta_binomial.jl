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

using Distributions
using Plots

#var = n*p*(1-p)*(s+n)/(s+1)
log_var = log(n) + log(p) + log(1-p) + log(s+n) - log(s+1)
exp(log_var)

var(binom)*k_var

log_var = log(n) + log(p) + log(1-p) + log(k_var)
exp(log_var)

beta_binom

n = 11000
p = 1 - exp(-1.0/n)
ρ = 1/5000
s = 1/ρ - 1 
beta_binom = BetaBinomial(n, p*s, (1-p)*s)
rnds = rand(beta_binom, 1000)
histogram(rnds) |> display
rnds = rand(Binomial(n, p), 1000)
histogram(rnds) |> display

# +
R0_list = [i*0.1 for i in 1:30]

n_obs = 5 # more than this value. 
for N in [10, 1000, 5_000] in
    p_5_betabin = []
    p_5_bin = []
    ρ =  1/5000
    for R0 in R0_list
        n = N - 1
        p = 1 - exp(-R0/n)
        s = 1/ρ - 1
        beta_binom = BetaBinomial(n, p*s, (1-p)*s)
        binom = Binomial(n, p)
        x = 0:n
        y = pdf.(beta_binom, x) 
        y_binom = pdf.(binom, x) 
        push!(p_5_betabin, 1 - sum(y[1:(n_obs+1)]))
        push!(p_5_bin, 1 - sum(y_binom[1:(n_obs+1)]))
    end
    pl = plot(title="Prob. of observing more than $(n_obs) cases \ngiven $(n+1) population with ρ=$(ρ)",
        fmt=:png,
    )
    plot!(pl, R0_list, p_5_betabin, xlabel="R0", ylabel="Prob.", label="BetaBinom")
    plot!(pl, R0_list, p_5_bin, label="Binom")
    pl |> display
end

# -


