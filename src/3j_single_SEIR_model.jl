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
using ShiftedArrays
using Optim

include("utils.jl")
include("model.jl")
# -

path = "../dt_tmp/spatial_params_agg230.jld2"
spatial_p = load(path)["data"]
@unpack pop, unvac, π_mat = spatial_p
nothing

pop_child = spatial_p.df[:,:value] |> sum
pop_whole = spatial_p.df[:,:value_whole] |> sum
prop_child = Float64(pop_child / pop_whole)

w_pop_cov = mean((pop .- unvac)./pop, weights(pop))
println(w_pop_cov)
N_tot = 100_000
N_unvac = round(N_tot * (1- w_pop_cov) ) |> Int64
Pop_whole = N_tot/prop_child
nothing

base_kwds = Dict(
    :N_tot=>N_tot, :N_unvac=>N_unvac,
    :R0=>14, :pc=>0.25, :Pop_whole=>Pop_whole
)

pars = SEIRModelParams(;base_kwds...) 
dump(pars)

Re0 = pars.R0 * (1-w_pop_cov)

# ## Check AFP incubation time 

pars_inc = SEIRModelParams(I0_init=0 ;base_kwds...)
dump(pars_inc)

pl = plot(xlim=[0,40], #ylim=[0, 0.1], 
    xticks=[0, 7, 15, 21, 25, 33, 39], fmt=:png,
    xlabel="day", ylabel="Prob",
)
dist = fill(0, pars_inc.days)
for i in 1:100
    rec, outcome = run_sim(pars_inc; rec_flag=true, pattern="check_incubation")
    dist += rec.Z_A5_6
    plot!(pl, 1:rec.days, rec.Z_A5_6/1000, label="New AFP", legend=:none, alpha=0.1)
    if i == 100
        plot!(pl, 1:rec.days, dist/sum(dist), label="New AFP", legend=:none)
    end
end
display(pl)

# # Check single run

Random.seed!(32)
rec, outcome = run_sim(pars; rec_flag=true)
println(outcome)
pl = plot(xlim=[0,365*3], ylim=[0, 200])
plot!(pl, 1:rec.days, rec.E, label="E")
plot!(pl, 1:rec.days, rec.I, label="I")
plot!(pl, 1:rec.days, rec.Z_A5_6*100, label="New AFP")
plot!(pl, 1:rec.days, rec.R, label="R")

function multiple_run(pars_)
    # Simulation part
    outcomes = []
    N_sim = 10000
    for i in 1:N_sim
        rec, outcome = run_sim(pars_; rec_flag=false)
        push!(outcomes, outcome)
    end
    df = DataFrame(outcomes)
    return df
end

# +
#function multiple_run_and_sim(pars)
#    df = multiple_run(pars)
#    first(df, 5) |> display 
#    # visualisation part
#    outcome_num_prop(df) |> display
#    
#    # Visualisation
#    days = pars.days
#    ts_ES = df[:, "t_ES"]
#    ts_AFP = df[:, "t_AFP"]
#    t_extinct = df[:, "t_extinct"]
#    cum_ES = cumulative_counts(ts_ES, days; prop=true)
#    cum_AFP = cumulative_counts(ts_AFP, days; prop=true)
#    cum_ES_cond = conditional_cumulative_prob(ts_ES, t_extinct, days)
#    cum_AFP_cond = conditional_cumulative_prob(ts_AFP, t_extinct, days)
#
#    pl1 = plot(
#        xlabel="Day", ylabel="Cumulative probability of first detection",
#        legend=(0.6, 0.7),
#    )
#    annotate!((0.07, 0.95), "(A)")
#
#    plot!(pl1, 1:days, cum_ES, label="Prob. via ES", color=1)
#    plot!(pl1, 1:days, cum_ES_cond, label="Conditional Prob. via ES", color=1, linestyle=:dashdot)
#    plot!(pl1, 1:days, cum_AFP, label="Prob. via AFP surv.", color=2)
#    plot!(pl1, 1:days, cum_AFP_cond, label="Conditional Prob. via AFP surv.", color=2, linestyle=:dashdot)
#
#    dif = leadtime_diff(df)
#    x = [1 for i in 1:length(dif)]
#    pl2 = violin(x, dif, xticks=:none, ylabel="Lead time of ES (day)", legend=:none)
#    boxplot!(pl2, x, dif, fillalpha=0.75)
#    annotate!((0.15, 0.95), "(B)")
#    l = @layout [a{0.75w} b]
#    pl = plot(pl1, pl2, 
#        fmt=:png, dpi=300, layout=l,
#        size=(800,500), left_margin=5Plots.mm,
#    )
#    display(pl)
#    savefig(pl, "../res/fig_single_prob_lead.png")
#    
#    # Discritive results.
#    println("Proportion of early detection", ((dif.> 0) |> sum)/length(dif)*100 )
#    
#    dif = leadtime_diff(df)
#    leadtime_diff_statistics(dif)
#    
#    # Final size distribution
#    println("Need of inclusion of other vaccinated individuals to the model.")
#    max_f = maximum(df[: ,:R_final_AFP])
#    pl2 = histogram(df[:, :R_final_AFP], bins=0:max_f, legend=:none, 
#        xlabel="Final AFP cases",
#        ylabel="Probability",
#        norm=true,
#        fmt=:png, dpi=300
#
#    )
#    df_fil = filter(x -> x.R_final_AFP != 0, df)
#    histogram!(pl2, df_fil[:, :R_final_AFP], bins=1:max_f, legend=:none,
#        inset = bbox(0.4Plots.w, 0.1Plots.h, 0.5Plots.w, 0.5Plots.h), 
#        subplot=2,
#        xlabel="Final AFP cases excluding 0",
#        ylabel="Probability",
#        norm=true,
#    )
#end

# +
#pars = SEIRModelParams(;base_kwds...)
#pars.pc = 1.0
#multiple_run_and_sim(pars)

# +
#pars_pc25 = SEIRModelParams(;base_kwds...)
#multiple_run_and_sim(pars_pc25)
# -

# # Effect of pc

function plot_cum!(pl, df; label="", color=:blue)
    cum_ES = cumulative_counts(df[:, "t_ES"], days; prop=true)
    cum_AFP = cumulative_counts(df[:, "t_AFP"], days; prop=true)
    plot!(pl, xlabel="Day", ylabel="Cumulative probability of first detection",
        fmt=:png
    )
    plot!(pl, 1:days, cum_ES, label="ES $(label)", color=color)
    plot!(pl, 1:days, cum_AFP, label="AFP $(label)", color=color, ls=:dash)
end

pars_pc1 = SEIRModelParams(; base_kwds...)
pars_pc1.pc = 1.0
pars_pc001 = SEIRModelParams(; base_kwds...)
pars_pc001.pc = 0.01
pars_pc025 = SEIRModelParams(; base_kwds...)
days = pars_pc025.days

# +
pl = plot()
df = multiple_run(pars_pc1)
plot_cum!(pl, df; label="pc100", color=:blue)

df = multiple_run(pars_pc025)
plot_cum!(pl, df; label="pc25", color=:red)

df = multiple_run(pars_pc001)
plot_cum!(pl, df; label="pc001", color=:green)
# -

# # Effect of sensitivity 

# +
pars_pc25 = SEIRModelParams(; base_kwds...)
pars_sens10_high = SEIRModelParams(
    ES_μ=-0.994, ES_σ=1.493
    ;base_kwds...)

pars_sens10_low = SEIRModelParams(
    ES_μ=3.611, ES_σ=1.493
    ;base_kwds...)

# +
pl = plot()
df = multiple_run(pars_pc25)
plot_cum!(pl, df; label="pc25", color=:blue)

df = multiple_run(pars_sens10_high)
plot_cum!(pl, df; label="10 high sens", color=:orange)

df = multiple_run(pars_sens10_low)
plot_cum!(pl, df; label="10 olower sens", color=:red)
# -

# # Effect of population size 

include("model.jl")

# +
N_tot = 100_000
pars_pc25 = SEIRModelParams(; base_kwds...)

N_tot = 10_000
pars_10000 = SEIRModelParams(
    N_tot=N_tot, N_unvac=round(N_tot * (1- w_pop_cov) ), 
    Pop_whole=N_tot/prop_child, 
    days=365*3, 
    R0=14, pc=0.25,
)
N_tot = 5_000
pars_5000 = SEIRModelParams(
    N_tot=N_tot, N_unvac=round(N_tot * (1- w_pop_cov) ), 
    days=365*3, 
    R0=14, pc=0.25,
    Pop_whole=N_tot/prop_child, 
)

# +
pl = plot(fmt=:png)
df = multiple_run(pars_pc25)
plot_cum!(pl, df; label="N_tot=100_000", color=:blue)

df = multiple_run(pars_10000)
plot_cum!(pl, df; label="N_tot=10_000", color=:orange)

df = multiple_run(pars_5000)
plot_cum!(pl, df; label="N_tot=5_000", color=:red)
# -

include("model.jl")

function obtain_detection_curve(pars)
    prev_sum = []
    for i in 1:5000
        rec, outcome = run_sim(pars; rec_flag=true)
        if isnan(outcome.t_ES) == false
            prev =  rec.I[outcome.t_ES]/(pars.Pop_whole*pars.pc)*100_000
            push!(prev_sum, prev)
        end
    end
    lognorm = fit(LogNormal, Float64.(prev_sum))
    return lognorm
end


lognorm_pc25 = obtain_detection_curve(pars_pc25)
lognorm_10000 = obtain_detection_curve(pars_10000)
lognorm_5000 = obtain_detection_curve(pars_5000)

x = 0:30
y_def = cdf(LogNormal(pars.ES_μ, pars.ES_σ), x)
y_pc25 = cdf(lognorm_pc25, x)
y_10000 = cdf(lognorm_10000, x)
y_5000 = cdf(lognorm_5000, x)
pl = plot(
    xlabel="Infects per 100,000 individuals", 
    ylabel="Probability of detection",
    title="Refit with Infects/pop within grid (not coverage)",
    fmt=:png,
)
plot!(x, y_def, label="default")
plot!(x, y_pc25, label="Pop with 100,000")
plot!(x, y_10000, label="Pop with 10,000")
plot!(x, y_5000, label="Pop with 5,000")

# ## Check the basic model behaviour

pl = plot(fmt=:png)
for i in 1:10
    rec, outcome = run_sim(pars_pc25; rec_flag=true)
    plot!(rec.S)
end
display(pl)

pl = plot(fmt=:png)
for i in 1:10
    rec, outcome = run_sim(pars_5000; rec_flag=true)
    plot!(rec.S)
end
display(pl)

# +
pl = plot(fmt=:png)
for i in 1:50
    rec, outcome = run_sim(pars_pc25; rec_flag=true)
    plot!(cumsum(rec.E))
end
display(pl)

pl = plot(fmt=:png)
for i in 1:50
    rec, outcome = run_sim(pars_5000; rec_flag=true)
    plot!(cumsum(rec.E))
end
display(pl)
# -






