using Accessors
using Glob
using Parameters
using ProtoStructs
using Random
using Serialization
using StatsBase

include("util.jl")

mutable struct SEIRMetaModelParams
    R0::Float64 
    γ1::Float64 
    γ2::Float64
    γ3::Float64 
    P_AFP::Float64 
    N0::Vector{Int64}
    I0_init::Int64 
    I0_ind::Int64
    days::Int64
    n_site::Int64
    α::Float64
    π_mat::Matrix{Float64}
    β::Float64
    function SEIRMetaModelParams(;
            R0=10.0, γ1=1/4.0, γ2=1/15.02, γ3=1/7.0,
            P_AFP=1/200, 
            N0=[1], I0_init=1, I0_ind=1,
            days=365, n_site=1,
            α=0.05, π_mat=fill(0,1,1)
        )
        new(R0, γ1, γ2, γ3, P_AFP,
            N0, I0_init, I0_ind, 
            days, n_site,
            α, π_mat, R0*γ2
            )
    end
end

mutable struct ESParams
    g::Float64
    P_test::Float64
    area::Vector{Bool}
    n_freq::Int64
    function ESParams(;g=0.5, P_test=0.97, area=[1], n_freq=30)
        new(g, P_test, area, n_freq)
    end
end
Base.copy(x::ESParams) = ESParams(g=x.g, P_test=x.P_test, area=x.area, n_freq=x.n_freq)

mutable struct AFPSurParams
    P_test::Float64
    P_sample::Float64
    P_H::Float64
    function AFPSurParams(;P_test=0.97, P_sample=0.8, P_H=0.9)
        new(P_test, P_sample, P_H)
    end
end

@proto struct SEIRMetaModelOneStep
    S::Vector{Int64}
    E::Vector{Int64}
    Ia::Vector{Int64}
    I_AFP::Vector{Int64}
    R::Vector{Int64}
    H_AFP::Vector{Int64}
end

@proto struct SEIRMetaModelRecord
    days::Int64
    n_site::Int64
    #S::Array{Int64, 2} = fill(-999, n_site, days)
    #E::Array{Int64, 2} = fill(-999, n_site, days)
    Ia::Array{Int64, 2} = fill(-999, n_site, days)
    I_AFP::Array{Int64, 2} = fill(-999, n_site, days)
    #R::Array{Int64, 2} = fill(-999, n_site, days)
    H_AFP::Array{Int64, 2} = fill(-999, n_site, days)
end

function set_values!(
        rec::SEIRMetaModelRecord, 
        model::SEIRMetaModelOneStep, 
        ind::Int64
    )
    @unpack S, E, Ia, I_AFP, R, H_AFP = model
    #rec.S[ind] = S
    #rec.E[ind] = E
    rec.Ia[:, ind] = Ia
    rec.I_AFP[:, ind] = I_AFP
    #rec.R[ind] = R
    rec.H_AFP[:, ind] = H_AFP
end

function initialize_model(params::SEIRMetaModelParams, rec_flag::Bool)
    @unpack N0, I0_init, n_site, N0, days = params
    I0_ind = wsample(1:n_site, N0)

    S = copy(N0)
    S[I0_ind] -= I0_init
    Ia = fill(0, n_site)
    Ia[I0_ind] += I0_init

    model = SEIRMetaModelOneStep(
        S = S,
        E = fill(0, n_site),
        Ia = Ia,
        I_AFP = fill(0, n_site),
        R = fill(0, n_site),
        H_AFP = fill(0, n_site),  # Newly seeking
    )

    if rec_flag == true
        rec = SEIRMetaModelRecord(days=days, n_site=n_site)
        set_values!(rec, model, 1)
    else
        # Empty record
        rec = SEIRMetaModelRecord(days=1, n_site=1)
    end

    return rec, model
end

function update_model(
        model::SEIRMetaModelOneStep, 
        params::SEIRMetaModelParams
    )::SEIRMetaModelOneStep

    @unpack γ1, γ2, γ3, P_AFP, days, β, N0, π_mat, α = params
    @unpack S, E, Ia, I_AFP, R = model

    ext = sum(π_mat .* (Ia .+ I_AFP)', dims=2)
    λ = β./N0.*( (1-α) .* (Ia .+ I_AFP) .+ α.*ext)
    new_inf = rand_binom.(S, 1 .- exp.(-λ))
    rec_E = rand_binom.(E, 1 .- exp(-γ1))
    rec_E_Ia = rand_binom.(rec_E, 1 - P_AFP)
    rec_E_I_AFP = rec_E .- rec_E_Ia

    rec_Ia = rand_binom.(Ia, 1 .- exp(-γ2))
    rec_I_AFP = rand_binom.(I_AFP, 1 .- exp(-γ3))

    S .+= .- new_inf
    E .+= .+ new_inf .- rec_E
    Ia .+=           .+ rec_E_Ia    .- rec_Ia
    I_AFP .+=        .+ rec_E_I_AFP .- rec_I_AFP
    R .+=                           .+ rec_Ia .+ rec_I_AFP
    model = SEIRMetaModelOneStep(
        S=S, E=E, Ia=Ia, I_AFP=I_AFP,R=R, 
        H_AFP=rec_I_AFP
    )
    return model
end

function run_sim(pars::SEIRMetaModelParams; rec_flag::Bool=false)
    @unpack γ1, γ2, γ3, P_AFP, days = pars 

    rec, model = initialize_model(pars, rec_flag)

    t_extinct = NaN
    R_final_num = NaN
    R_final_site = NaN
    
    for t in 2:days
        # Stop condition
        if sum(model.E) + sum(model.Ia) + sum(model.I_AFP) == 0
            t_extinct = isnan(t_extinct) == true ? t : t_extinct
            if rec_flag == true
                set_values!(rec, model, t)
                continue
            else
                break
            end
        end
        model = update_model(model, pars)
        if rec_flag == true
            set_values!(rec, model, t)
        end
    end
    R_final_num = sum(model.R)
    R_final_site = sum(model.R .> 0)
    R_final_AFP = sum(rec.H_AFP)
    outcome = (
        t_extinct=Float64(t_extinct), 
        R_final_num=R_final_num, 
        R_final_site=R_final_site,
        R_final_AFP=R_final_AFP,
    )
    return (rec=rec, outcome=outcome, pars=pars)
end

function AFP_surveillance(rec::SEIRMetaModelRecord, pars::AFPSurParams)::Float64
    @unpack P_test, P_sample, P_H = pars
    I_tot = sum(rec.H_AFP, dims=1)[1,: ]
    days = length(I_tot)
    t_AFP = NaN
    for t in 1:days
        w = rand_binom(I_tot[t], P_test*P_sample*P_H)
        if w > 0
            t_AFP = t
            break
        end
    end
    return t_AFP
end

function enviro_surveillance(rec::SEIRMetaModelRecord, pars::ESParams)::Float64
    @unpack g, P_test, area, n_freq = pars
    n_site, days = size(rec.Ia)
    
    # TODO: synchronize sampling time?
    nt = fill(0, n_site, days)
    st_ind = rand(1:n_freq)
    sample_ind = st_ind:n_freq:days
    nt[area, sample_ind] .= 1
    
    t_ES = NaN
    for t in sample_ind
        if nt[1, t] == 0
            continue
        end
        ωt = 1 .- exp.( -g .* (rec.Ia[:, t] .+ rec.I_AFP[:, t]))
        wt = rand_binom.(nt[:,t], ωt .* P_test)
        if sum(wt) > 0
            t_ES = t
            break
        end
    end
    return t_ES
end

function heatmap_meta_pop(I::Matrix{Int64})
    n_site, days = size(I)
    pl = plot(xlabel="Days", ylabel="Location")
    heatmap!(pl, 1:days, 1:n_site, log10.(I .+1) ) 
    return pl
end

function run_and_save_sim(pars::SEIRMetaModelParams, ; n_sim=10)
    now_str = get_today_time()
    path = "../dt_tmp/$now_str"
    mkdir(path)
    println(path)
    @showprogress for i in 1:n_sim
        res = run_sim(pars; rec_flag=true)
        serialize("$(path)/$(i).ser", res)
    end
end

function fetch_sim_paths(path::String) 
    return glob("$(path)/*.ser")
end

function collect_summary_statistics(
        path::String,  par_AFP::AFPSurParams, par_ES::ESParams
    )
    path_objs = fetch_sim_paths(path)
    n = length(path_objs)
    sim_res = DataFrame()
    @showprogress for i in 1:n
        res = deserialize(path_objs[i])
        rec, outcome, pars = res
        t_AFP = AFP_surveillance(rec, par_AFP)
        t_ES = enviro_surveillance(rec, par_ES)
        res_t = (
            t_AFP=t_AFP,
            t_ES=t_ES,
            outcome...
        )
        push!(sim_res,res_t)
    end
    return sim_res
end

"""
    sensitivity_ana_ES

Adaptor for sensitivity analysis.

# Arguments
- `f`:: Inside function, the sensitivity analysis setting is determined.
    and results are added to sim_res
"""
function sensitivity_ana_ES(
        path::String, par_AFP::AFPSurParams, par_ES::ESParams,
        f::Function
    )
    path_objs = fetch_sim_paths(path)
    n_sim = length(path_objs)
    sim_res = DataFrame()
    par_ES = copy(par_ES)
    @showprogress for i in 1:n_sim
        res = deserialize(path_objs[i])
        sim_res = f(res, par_AFP, par_ES, sim_res)
    end
    return sim_res
end

function sensitivity_ana_all(
        path::String, par_AFP::AFPSurParams, par_ES::ESParams,
    )
    path_objs = fetch_sim_paths(path)
    n_sim = length(path_objs)
    sim_res1 = DataFrame()
    sim_res2 = DataFrame()
    sim_res3 = DataFrame()
    @showprogress for i in 1:n_sim
        res = deserialize(path_objs[i])
        par_ES_tmp = copy(par_ES)
        sim_res1 = sensitivity_hazard(
            res, par_AFP, par_ES_tmp, sim_res1)

        par_ES_tmp = copy(par_ES)
        sim_res2 = sensitivity_frequency_sampling(
            res, par_AFP, par_ES_tmp, sim_res2)

        par_ES_tmp = copy(par_ES)
        sim_res3 = sensitivity_ES_catchment_area(
            res, par_AFP, par_ES_tmp, sim_res3)
    end
    return (sim_res1, sim_res2, sim_res3)
end

"""Function for `sensitivity_ana_ES`
"""
function sensitivity_ES_catchment_area(
        res::NamedTuple, par_AFP::AFPSurParams,
        par_ES::ESParams, sim_res::DataFrame,
    )
    rec, outcome, pars = res
    t_AFP = AFP_surveillance(rec, par_AFP)

    for ind_site in 1:pars.n_site
        area = fill(0, n_site)
        area[1:ind_site] .= 1
        par_ES.area = area
        t_ES = enviro_surveillance(rec, par_ES)
        res_t = (
            t_AFP=t_AFP,
            t_ES=t_ES,
            ind_site=ind_site,
            outcome...
        )
        push!(sim_res, res_t)
    end
    return sim_res
end

"""Function for `sensitivity_ana_ES`
"""
function sensitivity_hazard(
        res::NamedTuple, par_AFP::AFPSurParams,
        par_ES::ESParams, sim_res::DataFrame,
    )
    rec, outcome, pars = res
    t_AFP = AFP_surveillance(rec, par_AFP)
    pop90_list = [
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
        20, 30, 40, 50, 60, 70, 80, 90, 100, 
        200, 300, 400, 500, 600, 700, 800, 900, 1000
    ]
    for pop90 in pop90_list
        par_ES.g = - log(1 - 0.9)/pop90
        t_ES = enviro_surveillance(rec, par_ES)
        res_t = (
            t_AFP=t_AFP,
            t_ES=t_ES,
            pop90=pop90,
            g=par_ES.g,
            outcome...
        )
        push!(sim_res, res_t)
    end
    return sim_res
end

"""Function for `sensitivity_ana_ES`
"""
function sensitivity_frequency_sampling(
        res::NamedTuple, par_AFP::AFPSurParams,
        par_ES::ESParams, sim_res::DataFrame,
    )
    rec, outcome, pars = res
    t_AFP = AFP_surveillance(rec, par_AFP)
    freq_list = [
        7, 14, 21, 28, 30, 35, 42, 49, 56, 63, 70, 77, 84, 91
    ]
    for n_freq in freq_list
        par_ES.n_freq = n_freq
        t_ES = enviro_surveillance(rec, par_ES)
        res_t = (
            t_AFP=t_AFP,
            t_ES=t_ES,
            n_freq=n_freq,
            outcome...
        )
        push!(sim_res, res_t)
    end
    return sim_res
end

function detection_pattern_sensitivity(df_res::DataFrame, col)::DataFrame
    df_res[:, :pattern] = detect_pattern(df_res)
    tab = crosstab(df_res, col, :pattern)
    tab[:, "Detect"] = tab[:, "ES only"] .+ tab[:, "AFP only"] .+ tab[:, "Both"]
    return tab
end

function vis_detection_pattern(
        tab::DataFrame, n_sim::Int64, col::Union{String, Symbol};
        kwargs...
    )
    pl = plot(
        xlabel=col, 
        ylabel="Probability (%)",
        title="Prop. of polio detection (# of sim = $n_sim)", 
        fmt=:png;
        kwargs...
    )
    x = tab[:, col]
    plot!(pl, x, tab[:, "Detect"]./n_sim*100, 
        label="Detect by AFP or ES", marker=:circle, markersize=1)
    plot!(pl, x, tab[:, "Both"]./n_sim*100, label="Both")
    plot!(pl, x, tab[:, "ES only"]./n_sim*100, label="ES only detect")
    plot!(pl, x, tab[:, "AFP only"]./n_sim*100, label="AFP only detect")
    display(pl)
end

function leadtime_diff_sensitivity(df_res::DataFrame, col)::DataFrame
    df_res[:, "diff"] = df_res[:, "t_AFP"] - df_res[:, "t_ES"]
    cond = isnan.(df_res[:, "diff"]) .== false
    df_fil = df_res[cond, :]
    df_diff = DataFrame()
    for uni in unique(df_res[:, col])
        dfM = filter(x -> x[col] == uni, df_fil)
        qs = quantile_tuple(dfM[:, :diff])
        qs_new = (qs..., uni=uni)
        push!(df_diff, qs_new)
    end
    rename!(df_diff, :uni => col)
    return df_diff
end

function vis_leadtime_diff_sensitivity(df_diff::DataFrame, col; kwargs...)
    x = df_diff[:, col]
    pl = plot(
        title="Leadtime of ES detection", 
        xlabel=col, ylabel="Lead time (day)", 
        fmt=:png; 
        kwargs...
    )
    plot!(pl, x, df_diff[:, :q05], label="5th")
    plot!(pl, x, df_diff[:, :q25], label="25th")
    plot!(pl, x, df_diff[:, :q50], label="50th", marker=:circle, markersize=1)
    plot!(pl, x, df_diff[:, :q75], label="75th")
    plot!(pl, x, df_diff[:, :q95], label="95th")
    display(pl)
end