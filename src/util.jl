using Dates
using Distributions
using FreqTables
using GeoDataFrames
using Interpolations
using LinearAlgebra
using Pipe
using ProgressMeter
using Rasters
using Serialization
using Shapefile
using StatsPlots

function rand_binom(n, p)::Int64
    if (n == 0) | (p == 0)
        return 0
    else
        return rand(Binomial(n,p))
    end
end
get_today_time() = @pipe now() |> Dates.format(_, "yyyymmdd_HHMMSS")


"""
# Arguments
- `days`: Day for each datapoint.
- `virus` : Amount of virus shedding.
"""
function relative_virus_shedding(days::Vector, virus::Vector)
    itp = extrapolate(interpolate((days,), virus, Gridded(Linear())), Line())
    xs = [i for i in 0:100]
    ys = [itp(x) for x in xs]
    ys = ys/sum(ys)
    return Dict( x => y for (x,y) in zip(xs,ys) )
end

function hazard_function(I_new::Vector, gs::Dict, 
                         t::Int64; λ0=10)
    ind_end = minimum((100, t -1))
    λ = (I_new[t-s]*gs[s] for s in 0:ind_end) |> sum
    return λ*λ0
end

function outcome_num_prop(df::DataFrame)::DataFrame
    cond = isnan.(df) .== false
    num = sum.(eachcol(cond))
    prop = num ./size(df)[1]*100
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
    for (i,d) in enumerate(data)
        if isnan(d) == false
            d = Int64(d)
            Ys[i, d:end] .= 1
        end
    end
    return Ys
end

function create_indicater_U(data:: Vector, days)
    N_sim = length(data)
    U = fill(0, N_sim, days)
    for (i,d) in enumerate(data)
        if isnan(d) == true 
            d = days + 1
        end
        d = Int64(d)
        U[i, begin:(d-1)] .=1
    end
    return U
end

function cumulative_counts(data::Vector, days; prop::Bool=false)
    Ys = create_indicater_vals(data, days)
    cum = sum(Ys, dims=1)[1, :]
    if prop == true
        N_sim = length(data)
        cum = cum./N_sim
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
    print("\n% of early detect by ES: ", round(n_delay/length(diff)*100, digits=2))
    boxplot(diff, size=(300,500), fmt=:png) |> display
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

function convert_daily_to_weekly_obs(I_obs::Array{Int64, 2}; index=1)
    n_point = floor(size(I_obs)[2]/7) |> Int64
    n_sim = size(I_obs)[1]
    I_obs7 = fill(0, n_sim, n_point)
    for i in 1:n_point
        st = index + (i-1)*7 
        fin = index - 1 + i*7
        I_obs7[:, i] = sum(I_obs[:, st:fin], dims=2)
    end
    return I_obs7
end
