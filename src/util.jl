using Distributions
using Interpolations
using Rasters

rand_binom(n,p) = rand(Binomial(n,p))

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

function detect_pattern(df::DataFrame)::Dict
    cond = isnan.(df) .== false
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
    tp_sum |> countmap
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

function create_indicater_U(data:: Vector, day)
    N_sim = length(data)
    U = fill(0, N_sim, days)
    for (i,d) in enumerate(data)
        if isnan(d) == true 
            d = params.days + 1
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

function leadtime_diff(df::DataFrame)::Vector
    cond_df = isnan.(df) .== false
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

function raster_to_df(ras::Raster)::DataFrame
    coords = []
    lon = lookup(ras, X)
    lat = lookup(ras, Y)
    for lo in lon, la in lat
        push!(coords, (lon=lo, lat=la, value=ras[At(lo), At(la)]))
    end
    df = DataFrame(coords)
    return df
end

"""
    harvesine_dist(ϕ1, λ1, ϕ2, λ2)::Float64

Calculate distance between two geographical points using Vincenty Formula.

# Arguments
- `ϕ1::Float64`: latitude 1 [deg].
- `λ1::Float64`: longitude 1 [deg].
- `ϕ2::Float64`: latitude 2 [deg].
- `λ2::Float64`: longitude 2 [deg].
"""
function harversine_dist(ϕ1::Float64, λ1::Float64, ϕ2::Float64, λ2::Float64)::Float64
    r = 6371
    Δϕ = ϕ1 - ϕ2
    Δλ = λ1 - λ2
    a = sind(Δϕ/2.0)^2 + cosd(ϕ1)*cosd(ϕ2)*sind(Δλ/2.0)^2
    d = 2*r*atand(√a, √(1.0-a))*2π/360 # convert degree to radian
    return d
end