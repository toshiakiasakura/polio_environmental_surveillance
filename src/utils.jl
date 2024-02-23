#using Accessors
using CategoricalArrays
using CSV
using ColorSchemes
using DataFrames
using Dates
using Distributions
using FreqTables
using Glob
using GLM
using Interpolations
using JLD2
using LinearAlgebra
using Optim
using Parameters
using Pipe
using ProgressMeter
using PyFormattedStrings
using Random
#using Serialization
using ShiftedArrays
using SparseArrays
using SpecialFunctions
using StatsBase
using StatsPlots

function rand_binom(n, p)::Int64
    if (n == 0) | (p == 0)
        return 0
    else
        return rand(Binomial(n, p))
    end
end
get_today_time() = @pipe now() |> Dates.format(_, "yyyymmdd_HHMMSS")

function rand_beta_binom(n, p, ρ)::Int64
    s = 1 / ρ - 1
    if (n == 0) | (p == 0)
        return 0
    else
        return rand(BetaBinomial(n, p * s, (1 - p) * s))
    end
end

"""
# Arguments
- `days`: Day for each datapoint.
- `virus` : Amount of virus shedding.
"""
function relative_virus_shedding(days::Vector, virus::Vector)
    itp = extrapolate(interpolate((days,), virus, Gridded(Linear())), Line())
    xs = [i for i in 0:100]
    ys = [itp(x) for x in xs]
    ys = ys / sum(ys)
    return Dict(x => y for (x, y) in zip(xs, ys))
end

function hazard_function(I_new::Vector, gs::Dict,
    t::Int64; λ0=10)
    ind_end = minimum((100, t - 1))
    λ = (I_new[t-s] * gs[s] for s in 0:ind_end) |> sum
    return λ * λ0
end

function outcome_num_prop(df::DataFrame)::DataFrame
    cond = isnan.(df) .== false
    num = sum.(eachcol(cond))
    prop = num ./ size(df)[1] * 100
    DataFrame((num=num, prop=prop, outcome=names(df)))
end

function detect_pattern(df::DataFrame)::Vector
    cond = isnan.(df[:, [:t_ES, :t_AFP]]) .== false
    tp_sum = []
    for r in eachrow(cond)
        ES = r["t_ES"]
        AFP = r["t_AFP"]
        tp = "missing"
        if (ES == true) & (AFP == true)
            tp = "Both"
        elseif (ES == true) & (AFP == false)
            tp = "ES only"
        elseif (ES == false) & (AFP == true)
            tp = "AFP only"
        else
            tp = "Neither"
        end
        push!(tp_sum, tp)
    end
    return tp_sum
end

function create_indicater_vals(data::Vector, days)
    N_sim = length(data)
    Ys = fill(0, N_sim, days)
    for (i, d) in enumerate(data)
        if isnan(d) == false
            d = Int64(d)
            Ys[i, d:end] .= 1
        end
    end
    return Ys
end

function create_indicater_U(data::Vector, days)
    N_sim = length(data)
    U = fill(0, N_sim, days)
    for (i, d) in enumerate(data)
        if isnan(d) == true
            d = days + 1
        end
        d = Int64(d)
        U[i, begin:(d-1)] .= 1
    end
    return U
end

function cumulative_counts(data::Vector, days; prop::Bool=false)
    Ys = create_indicater_vals(data, days)
    cum = sum(Ys, dims=1)[1, :]
    if prop == true
        N_sim = length(data)
        cum = cum ./ N_sim
    end
    return cum
end

function conditional_cumulative_prob(ts::Vector, t_extinct::Vector, days)
    Ys = create_indicater_vals(ts, days)
    U = create_indicater_U(t_extinct, days)
    cum = sum(Ys .* U, dims=1) ./ sum(U, dims=1)
    return cum[1, :]
end

function vis_cumulative_prob(df::DataFrame, days::Int64; kwargs...)
    ts_ES = df[:, "t_ES"]
    ts_AFP = df[:, "t_AFP"]
    t_extinct = df[:, "t_extinct"]
    cum_ES = cumulative_counts(ts_ES, days; prop=true)
    cum_AFP = cumulative_counts(ts_AFP, days; prop=true)
    pl = plot(ylabel="Probability", xlabel="Days";
        kwargs...
    )
    plot!(pl, 1:pars.days, cum_ES, label="ES", fmt=:png)
    plot!(pl, 1:pars.days, cum_AFP, label="AFP")

    cum_ES = conditional_cumulative_prob(ts_ES, t_extinct, days)
    cum_AFP = conditional_cumulative_prob(ts_AFP, t_extinct, days)
    plot!(pl, 1:pars.days, cum_ES, label="Pc, ES")
    plot!(pl, 1:pars.days, cum_AFP, label="Pc, AFP")
end

function leadtime_diff(df::DataFrame)::Vector
    cond_df = isnan.(df[:, ["t_ES", "t_AFP"]]) .== false
    cond = cond_df[:, "t_ES"] .& cond_df[:, "t_AFP"]
    dfM = df[cond, :]
    diff = dfM[:, "t_AFP"] - dfM[:, "t_ES"]
    return diff
end

function leadtime_diff_statistics(diff::Vector)
    describe(diff)
    n_delay = (diff .> 0) |> sum
    print("\n% of early detect by ES: ", round(n_delay / length(diff) * 100, digits=2))
    boxplot(diff, size=(300, 500), fmt=:png) |> display
end


function crosstab(df::DataFrame, row, col)::DataFrame
    tab = freqtable(df, row, col)
    df_tab = DataFrame(tab |> Array, names(tab)[2] .|> string)
    df_tab[:, row] = names(tab)[1]
    return df_tab
end

quantile_tuple(x) = (
    q05=quantile(x, 0.05),
    q25=quantile(x, 0.25),
    q50=quantile(x, 0.50),
    q75=quantile(x, 0.75),
    q95=quantile(x, 0.95),
)

function convert_daily_to_weekly_obs(I_obs::Array{Int64,2}; index=1)
    n_point = floor(size(I_obs)[2] / 7) |> Int64
    n_sim = size(I_obs)[1]
    I_obs7 = fill(0, n_sim, n_point)
    for i in 1:n_point
        st = index + (i - 1) * 7
        fin = index - 1 + i * 7
        I_obs7[:, i] = sum(I_obs[:, st:fin], dims=2)
    end
    return I_obs7
end

function obtain_EVP()
    path = "../data/zaf_OPV_HEXA_vaccine_coverage_2020.csv"
    df_vac = CSV.read(path, DataFrame)
    VE = 0.63
    VE1 = 1 - (1 - VE)
    VE2 = 1 - (1 - VE)^2
    VE3 = 1 - (1 - VE)^3
    VE4 = 1 - (1 - VE)^4

    CV4 = df_vac.HEXA4
    CV3 = df_vac.HEXA3
    CV2 = df_vac.HEXA2
    CV1 = df_vac.HEXA1

    dif3 = [i < 0 ? 0.0 : i for i in CV3 .- CV4]
    dif2 = [i < 0 ? 0.0 : i for i in CV2 .- CV3]
    dif1 = [i < 0 ? 0.0 : i for i in CV1 .- CV2]
    EVP = CV4 .* VE4 .+ dif3 .* VE3 .+ dif2 .* VE2 .+ dif1 .* VE1
    df_vac[!, "EVP"] = EVP
    return df_vac
end

function borel_tanner_dist(x, R)
    log_p = (x - 2) * log(x) + (x - 1) * log(R) - x * R - loggamma(x)
    return exp(log_p)
end

function observe_more_than_y_cases_in_final_dist(y, R)
    p = 1
    for i in 1:y
        p -= borel_tanner_dist(i, R)
    end
    return p
end

_standard_normal = Normal()
function probit_func(x::Float64; a::Float64=1.0, b::Float64=1.0)
    return cdf(_standard_normal, a + b * x)
end

"""Read the ES related spatial parameter file.
This file choice is determined by the way to layout ES.

Args:
- `ES_pattern`: Takes `ES_population_size` or `ES_mozambique_imp_risk`

"""
function read_spatial_params_file(ES_pattern)
    spatial_p = :none
    if (ES_pattern == "ES_population_size")
        path = "../dt_tmp/spatial_params_agg230.jld2"
        spatial_p = load(path)["data"]
    elseif ES_pattern == "ES_mozambique_imp_risk"
        path = "../dt_tmp/spatial_params_agg230_moz_sorted.jld2"
        spatial_p = load(path)["data"]
    else
        error("Specify type from the follwoing: " *
              "ES_population_size, ES_mozambique_imp_risk"
        )
    end
    return spatial_p
end

function read_imp_ws_data(pattern, ES_pattern)
    if pattern == "population_size"
        imp_ws = [1]
    elseif pattern == "airport"
        if ES_pattern == "ES_population_size"
            imp_ws = CSV.read("../data/imp_ws_airport.csv", DataFrame, header=false)[:,1]
        elseif ES_pattern == "ES_mozambique_imp_risk"
            imp_ws = CSV.read("../data/imp_ws_airport_sorted.csv", DataFrame, header=false)[:,1]
        end
    elseif pattern == "mozambique"
        if ES_pattern == "ES_population_size"
            imp_ws = CSV.read("../data/imp_ws_moz.csv", DataFrame, header=false)[:,1]
        elseif ES_pattern == "ES_mozambique_imp_risk"
            imp_ws = CSV.read("../data/imp_ws_moz_sorted.csv", DataFrame, header=false)[:,1]
        end
    end
    return imp_ws
end