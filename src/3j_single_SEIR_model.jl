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
include("model.jl")
include("visualise_fig.jl")

path = "../dt_tmp/spatial_params_agg230.jld2"
spatial_p = load(path)["data"]
@unpack pop, unvac, π_mat = spatial_p
nothing

pop_child = spatial_p.df[:,:value] |> sum
pop_whole = spatial_p.df[:,:value_whole] |> sum
prop_child = Float64(pop_child / pop_whole)
println("""
Proportion of children: $(prop_child)
with SA national populaion size for children and whole""")

# ## Weighted population coverage

w_pop_cov = mean((pop .- unvac)./pop, weights(pop))
println("Weighted vaccinated populations: ", w_pop_cov)

# +
Nc = 100_000
N_unvac = round(Nc * (1- w_pop_cov) ) |> Int64
# Whole population is defiened to account for 9.5% of children under 5 years old.
Pop_whole = Nc/prop_child

base_kwds = Dict(
    :Nc=>Nc, :N_unvac=>N_unvac,
    :pc=>0.25, :Pop_whole=>Pop_whole
)
nothing
# -

pars = SEIRModelParams(;base_kwds...)
dump(pars)

Re0 = pars.R0 * (1-w_pop_cov)
println("Effective reproduction number from the national level:\n", Re0)

# ## Check AFP incubation time

pars_inc = SEIRModelParams(I0_init=0 ;base_kwds...)
dump(pars_inc)

# +
pl = plot(xlim=[0,40], #ylim=[0, 0.1],
    xticks=[0, 7, 15, 21, 25, 33, 39], fmt=:png,
    xlabel="day", ylabel="Prob",
)
dist = fill(0, pars_inc.days)

n_sim = 100
n_init_inf = 1000
for i in 1:n_sim
    # when check_incubation pattern, the initial infectious individuals are given for 1000.
    rec, outcome = run_sim(pars_inc; rec_flag=true, pattern="check_incubation")
    dist += rec.Z_A5_6
    plot!(pl, 1:rec.days, rec.Z_A5_6/n_init_inf, label="New AFP", legend=:none, alpha=0.1)
    if i == n_sim
        plot!(pl, 1:rec.days, dist/n_init_inf/n_sim, label="New AFP", legend=:none)
    end
end
display(pl)
# -

# # Check single run

Random.seed!(43) # 32
rec, outcome = run_sim(pars; rec_flag=true)
println(outcome)
pl = plot(xlim=[0,365*3], ylim=[0, 200])
plot!(pl, 1:rec.days, rec.E, label="E")
plot!(pl, 1:rec.days, rec.Ic + rec.Inc, label="I")
plot!(pl, 1:rec.days, rec.Z_A5_6*100, label="New AFP")
plot!(pl, 1:rec.days, rec.R, label="R")

function multiple_run(pars_)
    Random.seed!(1234)
    # Simulation part
    outcomes = []
    N_sim = 20000
    for i in 1:N_sim
        rec, outcome = run_sim(pars_; rec_flag=false)
        push!(outcomes, outcome)
    end
    df = DataFrame(outcomes)
    return df
end

# # Effect of pc

cm = palette(:tab10)
colors3 = cm[1:3]
nothing

# Baseline
pars_base = SEIRModelParams(; base_kwds...)
days = pars_base.days

# +
# pc sensitivity
pars_pc100 = SEIRModelParams(; base_kwds...)
pars_pc100.pc = 1.0

pars_pc001 = SEIRModelParams(; base_kwds...)
pars_pc001.pc = 0.01

pars_pc_lis = [pars_pc100, pars_base, pars_pc001]
nothing

# +
# ES detection sensitivity
pars_sens10_high = SEIRModelParams(
    ES_μ=-0.994, ES_σ=1.493
    ;base_kwds...)

pars_sens10_low = SEIRModelParams(
    ES_μ=3.611, ES_σ=1.493
    ;base_kwds...)

pars_sens_lis = [pars_sens10_low, pars_base, pars_sens10_high]
nothing
# -

# Population size
Nc = 20_000
pars_20000 = SEIRModelParams(
    Nc=Nc, N_unvac=round(Nc * (1- w_pop_cov) ),
    Pop_whole=Nc/prop_child,
    days=365*3,
    R0=14, pc=0.25,
)
Nc = 5_000
pars_5000 = SEIRModelParams(
    Nc=Nc, N_unvac=round(Nc * (1- w_pop_cov) ),
    days=365*3,
    R0=14, pc=0.25,
    Pop_whole=Nc/prop_child,
)
pars_pop_lis = [pars_5000, pars_20000, pars_base]
nothing

# +
# Sampling frequency
pars_freq_1 = SEIRModelParams(ES_n_freq=1 ; base_kwds...)
pars_freq_60 = SEIRModelParams(ES_n_freq=60 ; base_kwds...)

pars_freq_lis = [pars_freq_1, pars_base, pars_freq_60]

# R0
pars_R0_10 = SEIRModelParams(R0=10; base_kwds...)
pars_R0_18 = SEIRModelParams(R0=18; base_kwds...)

pars_R0_lis = [pars_R0_10, pars_base, pars_R0_18]
nothing
# -

# ## Figure settings

# +
# Summarise the figures
pl_R0 = draw_cumulative_incidence(pars_R0_lis, colors3;
    legendtitle="R0", labels=["10", "14", "18"],
    xlabel="Day",
    legend=:topright,
)
plot!(pl_R0, left_margin=5Plots.mm)

pl_freq = draw_cumulative_incidence(pars_freq_lis, colors3;
    legendtitle="Sampling freq.", labels=["1 day", "30 day", "60 day"],
    xlabel="Day",
    legend=:topright,
)

pl_sens = draw_cumulative_incidence(pars_sens_lis, colors3;
    legendtitle="ES sensitivity", labels=["10x lower", "Baseline", "10x higher"],
    xlabel="Day",
    legend=:topright,
)

pl_pop  = draw_cumulative_incidence(pars_pop_lis, colors3;
    legendtitle="Nc", labels=["5000", "20,000", "100,000"],
    xlabel="Day",
    legend=:topright,
)

pl_pc = draw_cumulative_incidence(pars_pc_lis, colors3;
    xlabel="Day",
    legendtitle="pc", labels=["100%", "  25%", "    1%"],
)
nothing
# -

# # Detection pattern figure

function plot_prop_any_and_prop_pattern(
        par_lis, cate, xlabel
    )
    p_any_lis, freq_lis = obtain_p_any_and_freq_list(par_lis)
    pl_bar = bar(cate, p_any_lis.*100,
        ylabel="Simulated probability \nof detection (%)",
        xlabel=xlabel,
        label=:none
    )
    pl_grpbar = plot_groupedbar(cate, freq_lis; xlabel=xlabel,
        ylabel="Proportion of \ndetection pattern (%)",
    )
    return(pl_bar, pl_grpbar)
end

include("visualise_fig.jl")

# +
pl_bar_R0, pl_grpbar_R0 = plot_prop_any_and_prop_pattern(
    pars_R0_lis, ["10", "14", "18"], "R0")
pl_bar_freq, pl_grpbar_freq = plot_prop_any_and_prop_pattern(
    pars_freq_lis, ["1 day", "30 day", "60 day"], "Sampling freq.")
pl_bar_sens, pl_grpbar_sens = plot_prop_any_and_prop_pattern(
    pars_sens_lis, ["10x lower", "Baseline", "10x higher"], "ES sensitivity")
pl_bar_pop, pl_grpbar_pop = plot_prop_any_and_prop_pattern(
    pars_pop_lis, ["5,000", "20,000", "100,000"], "Nc")
pl_bar_pc, pl_grpbar_pc= plot_prop_any_and_prop_pattern(
    pars_pc_lis, ["100%", "  25%", "    1%"], "pc")

nothing

# +
plot!(pl_R0, top_margin=10Plots.mm)
annotate!(pl_R0, (0.05, 1.0), text("A", :left, :bottom, :black, 24))
annotate!(pl_bar_R0, (0.05, 1.0), text("B", :left, :bottom, :black, 24))
annotate!(pl_grpbar_R0, (0.05, 1.0), text("C", :left, :bottom, :black, 24))

prefixes = ["", "bar_", "grpbar_"]
pls = []
for var in ["R0", "freq", "sens", "pop", "pc"]
    for prefix in prefixes
        push!(pls, eval(Meta.parse("pl_$(prefix)$(var)")))
    end
end
layout = @layout [a{0.5w} b c; d e f; g h i; j k l; m n o]
plot(pls..., layout=layout,
    fmt=:png, dpi=300, size=(1200, 250*5),
)
# -

bin_labels = ["AFP only", "<-60 LT", "-60 ~ -1 LT", "0 ~ 59 LT", "≥60 LT", "ES only"]
colors = discretise_balance_color(bin_labels)
pl = groupedbar([""], [25 12.5 12.5 12.5 12.5 25],
    bar_position=:stack,
    labels=reshape(bin_labels, 1, 6),
    color=colors[:, end:-1:begin],
    legend=:none,
    ylim=[0,100],
    ylabel="Proportion of detection pattern (%)",
    ylabelfontsize=22,
    ytickfontsize=18,
    right_margin=30Plots.mm,
    left_margin=10Plots.mm,
    fmt=:png, dpi=200,
    bar_width=0.7,
    size=(300,600),
)



# # Fit a simulated data to re-estimate the parameters

include("model.jl")

function obtain_detection_curve(pars)
    prev_sum = []
    for i in 1:5000
        rec, outcome = run_sim(pars; rec_flag=true)
        if isnan(outcome.t_ES) == false
            prev =  rec.Ic[outcome.t_ES]/(pars.Pop_whole*pars.pc)*100_000
            push!(prev_sum, prev)
        end
    end
    lognorm = fit(LogNormal, Float64.(prev_sum))
    return lognorm
end


lognorm_pc25 = obtain_detection_curve(pars_base)
lognorm_20000 = obtain_detection_curve(pars_20000)
lognorm_5000 = obtain_detection_curve(pars_5000)

# +
x = 0:0.1:30
y_def = cdf(LogNormal(pars.ES_μ, pars.ES_σ), x)
y_pc25 = cdf(lognorm_pc25, x)
y_10000 = cdf(lognorm_20000, x)
y_5000 = cdf(lognorm_5000, x)

ytick_v = [0, 0.25, 0.5, 0.75, 1.0] 
ytick_l = ["0", "25", "50", "75", "100"]
pl = plot(
    xlabel="Infectious individuals per 100,000 population",
    ylabel="ES sensitivity (%)",
    yticks=(ytick_v, ytick_l),
    #title="Refit with Infects/pop within grid (not coverage)",
    fmt=:png, dpi=300
)
plot!(x, y_def, label="Assumed ES sensitivity", lw=2)
plot!(x, y_pc25, label="Nc = 100,000")
plot!(x, y_10000, label="Nc = 20,000")
plot!(x, y_5000, label="Nc = 5,000")
# -

# ## Check the basic model behaviour

pl = plot(fmt=:png)
for i in 1:10
    rec, outcome = run_sim(pars_base; rec_flag=true)
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
    rec, outcome = run_sim(pars_base; rec_flag=true)
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






