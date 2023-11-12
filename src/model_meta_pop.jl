using Accessors
using CategoricalArrays
using Glob
using Parameters
using ProtoStructs
using Random
using Serialization
using SparseArrays
using StatsBase

include("util.jl")
mutable struct SEIRMetaModelParams
    R0::Float64 
    γ1::Float64 
    γ2::Float64
    σ::Float64 
    P_AFP::Float64 
    N_tot::Vector{Int64}
    N_unvac::Vector{Int64}
    I0_init::Int64 
    I0_ind::Int64
    days::Int64
    n_site::Int64
    α::Float64
    π_mat::Matrix{Float64}
    μ::Float64
    β::Float64
    pc::Float64
    imp_ws::Vector{Float64}
    function SEIRMetaModelParams(;
            R0=10.0, γ1=1/4.0, γ2=1/15.02, σ=0.3291,
            P_AFP=1/200, 
            N_tot=[1], N_unvac=[1], I0_init=1, I0_ind=1,
            days=365, n_site=1,
            α=0.05, π_mat=fill(0,1,1), μ=1/(5*365),
            pc=1.0, imp_ws=[1.0]
        )
        new(R0, γ1, γ2, σ, P_AFP,
            N_tot, N_unvac, I0_init, I0_ind, 
            days, n_site,
            α, π_mat, μ, R0*γ2, 
            pc, imp_ws
            )
    end
end

mutable struct ESParams
    g::Float64
    P_test::Float64
    area::Vector{Bool}
    n_freq::Int64
    function ESParams(;g=0.23, P_test=0.97, area=[1], n_freq=30)
        new(g, P_test, area, n_freq)
    end
end
Base.copy(x::ESParams) = ESParams(g=x.g, P_test=x.P_test, area=x.area, n_freq=x.n_freq)

mutable struct AFPSurParams
    P_test::Float64
    P_sample::Float64
    P_H::Float64
    function AFPSurParams(;P_test=0.97, P_sample=0.53, P_H=0.9)
        new(P_test, P_sample, P_H)
    end
end

@proto struct SEIRMetaModelOneStep
    S::Vector{Int64}
    E::Vector{Int64}
    Ic::Vector{Int64}
    Inc::Vector{Int64}
    R::Vector{Int64}
    A::Array{Int64, 2}
    Z_A5_6::Int64
    Z_S_E::Int64
end

@proto struct SEIRMetaModelRecord
    days::Int64
    n_site::Int64
    I::Array{Int64, 2} = fill(0, n_site, days)
    Z_A5_6::Array{Int64, 1} = fill(0, days)
    Z_S_E::Array{Int64, 1} = fill(0, days)
end

@proto struct SEIRMetaModelRecordSparse
    days::Int64
    n_site::Int64
    I::SparseMatrixCSC{Int64, Int64} 
    Z_A5_6::SparseVector{Int64, Int64} 
    Z_S_E::SparseVector{Int64, Int64} 
end

function set_values!(
        rec::SEIRMetaModelRecord, 
        model::SEIRMetaModelOneStep, 
        ind::Int64
    )
    @unpack S, E, Ic, Inc, Z_A5_6, Z_S_E = model
    rec.I[:, ind] = Ic
    rec.Z_A5_6[ind] = Z_A5_6
    rec.Z_S_E[ind] = Z_S_E
end

function initialize_model(params::SEIRMetaModelParams, rec_flag::Bool;
     pattern=""
     )
    @unpack I0_init, n_site, N_tot, N_unvac, days, pc, imp_ws = params
    if pattern == "population_size"
        I0_ind = wsample(1:n_site, N_tot)
    elseif pattern == "airport"
        # For international airport introduction version.
        I0_ind = wsample([11, 7, 62], [4_342_611, 1_156_996, 188_243])
    elseif pattern == "mozambique"
        I0_ind = wsample(1:n_site, imp_ws)
    else
        error("Specify pattern: 'population_size', 'airport' or 'mozambique'")
    end

    S = copy(N_unvac) 
    S[I0_ind] -= I0_init

    # initialise Ic and Inc.
    Ic_init = rand_binom(I0_init, pc)
    Inc_init = I0_init - Ic_init
    Ic = fill(0, n_site) 
    Ic[I0_ind] += Ic_init
    Inc = fill(0, n_site) 
    Inc[I0_ind] += Inc_init 

    model = SEIRMetaModelOneStep(
        S = S,
        E = fill(0, n_site),
        Ic = Ic,
        Inc = Inc,
        R = fill(0, n_site),
        A = fill(0, n_site, 6),
        Z_A5_6 = 0,
        Z_S_E = 0,
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

    @unpack γ1, γ2, σ, μ, P_AFP, days, β, N_tot, N_unvac, π_mat, α, n_site, pc = params
    @unpack S, E, Ic, Inc, R, A = model

    new_born = rand_binom.(N_unvac, μ)
    ext = sum(π_mat' .* (Ic .+ Inc)', dims=2)[:, 1]
    λ = β./N_tot.*( (1-α) .* (Ic + Inc) .+ α.*ext)

    rem_S = rand_binom.(S, 1 .- exp.(-λ .- μ))
    new_E = rand_binom.(rem_S, λ./(λ .+ μ))

    rem_E = rand_binom.(E, 1 .- exp(-γ1 - μ))
    new_I = rand_binom.(rem_E, γ1/(γ1 + μ))
    new_Ic = rand_binom.(new_I, pc)
    new_Inc = new_I .- new_Ic

    rem_Ic = rand_binom.(Ic, 1 .- exp(-γ2 - μ))
    rem_Inc = rand_binom.(Inc, 1 .- exp(-γ2 - μ))
    new_R = rand_binom.(rem_Ic .+ rem_Inc, γ2/(γ2 + μ))

    new_AFP = fill(0, n_site, 6)
    new_AFP[:, 1] = rand_binom.(new_E, P_AFP)
    rate = 1 - exp(-σ)
    for i in 1:5
        new_AFP[:, i+1] .= rand_binom.(A[:, i], rate)
    end

    S .+= new_born .- rem_S
    E .+=          .+ new_E .- rem_E
    Ic .+=                  .+ new_Ic  .- rem_Ic
    Inc .+=                   .+ new_Inc .- rem_Inc
    R .+=                              .+ new_R

    for i in 1:5
        A[:, i] .+= new_AFP[:, i] .- new_AFP[:, i+1]
    end
    A[:, 6] .+= new_AFP[:, 6]

    model = SEIRMetaModelOneStep(
        S=S, E=E, Ic=Ic, Inc=Inc, R=R, A=A,
        Z_A5_6 = sum(new_AFP[:, 6]), Z_S_E=sum(new_E) .|> Int64,
    )
    return model
end

function run_sim(pars::SEIRMetaModelParams; 
        rec_flag::Bool=false, pattern=""
        )
    @unpack γ1, γ2, σ, P_AFP, days = pars 

    rec, model = initialize_model(pars, rec_flag, pattern=pattern)

    t_extinct = NaN
    R_final_num = NaN
    R_final_site = NaN
    
    for t in 2:days
        # Stop condition
        if sum(model.E) + sum(model.Ic) + sum(model.Inc) + sum(model.A[:, 1:5]) == 0
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
    R_final_AFP = sum(rec.Z_A5_6)
    outcome = (
        t_extinct=Float64(t_extinct), 
        R_final_num=R_final_num, 
        R_final_site=R_final_site,
        R_final_AFP=R_final_AFP,
    )
    rec_sparse = SEIRMetaModelRecordSparse(
        days=rec.days,
        n_site=rec.n_site,
        I=sparse(rec.I),
        Z_A5_6=sparse(rec.Z_A5_6),
        Z_S_E=sparse(rec.Z_S_E),
    )
    return (rec=rec_sparse, outcome=outcome, pars=pars)
end

function AFP_surveillance(rec::SEIRMetaModelRecordSparse, pars::AFPSurParams)::Float64
    @unpack P_test, P_sample, P_H = pars
    I_tot = rec.Z_A5_6
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

function enviro_surveillance(rec::SEIRMetaModelRecordSparse, pars::ESParams)::Float64
    @unpack g, P_test, area, n_freq = pars
    n_site, days = size(rec.I)
    
    nt = fill(0, n_site, days)
    st_ind = rand(1:n_freq)
    sample_ind = st_ind:n_freq:days
    nt[area, sample_ind] .= 1
    
    t_ES = NaN
    for t in sample_ind
        if nt[1, t] == 0
            continue
        end
        ωt = 1 .- exp.( -g .* rec.I[:, t] )
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

function run_and_save_sim(pars::SEIRMetaModelParams; n_sim=10, pattern="")
    now_str = get_today_time()
    path = "../dt_tmp/$now_str"
    mkdir(path)
    println(path)
    @showprogress for i in 1:n_sim
        res = run_sim(pars; rec_flag=true, pattern=pattern)
        serialize("$(path)/$(i).ser", res)
    end
    return path
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
        R_inf_ES = isnan(t_ES) == false ? sum(rec.Z_S_E[begin:Int64(t_ES)]) : NaN
        res_t = (
            t_AFP=t_AFP,
            t_ES=t_ES,
            R_inf_ES = R_inf_ES,
            outcome...
        )
        push!(sim_res,res_t, promote=true)
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
        path::String, par_AFP::AFPSurParams, par_ES::ESParams; 
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

function obtain_ES_sensitivity_index(pop::Vector, inc_prop::Float64)::Vector
    n_site = length(pop)
    cum_prop = cumsum(pop/sum(pop))
    pre_prop = 0
    index = []
    for i in 1:n_site
        if i > 1
            dif = cum_prop[i] - pre_prop
            if dif < inc_prop
                continue
            end
        end
        push!(index, i)
        pre_prop = cum_prop[i]
    end
    return index
end

"""Function for `sensitivity_ana_ES`
"""
function sensitivity_ES_catchment_area(
        res::NamedTuple, par_AFP::AFPSurParams,
        par_ES::ESParams, sim_res::DataFrame
    )
    rec, outcome, pars = res
    t_AFP = AFP_surveillance(rec, par_AFP)

    sens_index = obtain_ES_sensitivity_index(res.pars.N_tot, 0.01) # 1% increment.
    for ind_site in sens_index
        area = fill(0, res.pars.n_site)
        area[1:ind_site] .= 1
        par_ES.area = area
        t_ES = enviro_surveillance(rec, par_ES)
        R_inf_ES = isnan(t_ES) == false ? sum(rec.Z_S_E[begin:Int64(t_ES)]) : NaN
        res_t = (
            t_AFP=t_AFP,
            t_ES=t_ES,
            R_inf_ES=R_inf_ES,
            ind_site=ind_site,
            outcome...
        )
        push!(sim_res, res_t, promote=true)
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
        1, 2, 3, 5, 7, 
        10, 20, 30, 50, 70, 
        100, 200, 300, 500, 700, 1000
    ]
    for pop90 in pop90_list
        par_ES.g = - log(1 - 0.9)/pop90
        t_ES = enviro_surveillance(rec, par_ES)
        R_inf_ES = isnan(t_ES) == false ? sum(rec.Z_S_E[begin:Int64(t_ES)]) : NaN
        res_t = (
            t_AFP=t_AFP,
            t_ES=t_ES,
            R_inf_ES=R_inf_ES,
            pop90=pop90,
            g=par_ES.g,
            outcome...
        )
        push!(sim_res, res_t, promote=true)
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
        1, 7, 14, 21, 28, 30, 35, 42, 49, 56,
    ]
    for n_freq in freq_list
        par_ES.n_freq = n_freq
        t_ES = enviro_surveillance(rec, par_ES)
        R_inf_ES = isnan(t_ES) == false ? sum(rec.Z_S_E[begin:Int64(t_ES)]) : NaN
        res_t = (
            t_AFP=t_AFP,
            t_ES=t_ES,
            R_inf_ES=R_inf_ES,
            n_freq=n_freq,
            outcome...
        )
        push!(sim_res, res_t, promote=true)
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
        ylabel="Detection pattern (%)",
        title="Prop. of polio detection (# of sim = $n_sim)", 
        ylim=[0,100],
        fmt=:png;
        kwargs...
    )
    x = tab[:, col]
    #plot!(pl, x, tab[:, "Detect"]./n_sim*100, 
    #    label="AFP surv. or ES", marker=:xcross, markersize=3)
    det = tab[:, "Detect"]
    plot!(pl, x, tab[:, "Both"]./det.*100, label="AFP surv. and ES", 
        marker=:circle, markersize=2, markerstrokewidth = 0,
    )
    plot!(pl, x, tab[:, "ES only"]./det.*100, label="ES only",
        marker=:circle, markersize=2, markerstrokewidth = 0,
    )
    plot!(pl, x, tab[:, "AFP only"]./det.*100, label="AFP surv. only",
        marker=:circle, markersize=2, markerstrokewidth = 0,
    )
    return pl
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
        xlabel=col, ylabel="Lead time of ES (day)", 
        fmt=:png; 
        kwargs...
    )
    plot!(pl, x, df_diff[:, :q05], label="5th", 
        color="blue", alpha=0.50, linestyle=:dot,
        )
    plot!(pl, x, df_diff[:, :q25], label="25th", 
        color="blue", alpha=0.75, linestyle=:dashdot,
    )
    plot!(pl, x, df_diff[:, :q50], label="50th", 
        color="blue", alpha=1.0,
        marker=:circle, markersize=2, markerstrokewidth = 0,
        )
    plot!(pl, x, df_diff[:, :q75], label="75th", 
        color="blue", alpha=0.75, linestyle=:dashdot,
    )
    plot!(pl, x, df_diff[:, :q95], label="95th", 
        color="blue", alpha=0.50, linestyle=:dot,
    )
    #hline!(pl, [0], color="black", linestyle=:dot, label=:none)
    return pl
end

function early_detect_prob_sensitivity(df_res::DataFrame, col)::DataFrame
    ES_nan = isnan.(df_res[:, "t_ES"])
    AFP_nan = isnan.(df_res[:, "t_AFP"])
    df_res[!, "ES_only"] .= (ES_nan .== false) .& AFP_nan
    df_res[!, "Either"] .= (ES_nan .== false) .| (AFP_nan .== false)
    df_res[!, "diff"] .= df_res[:, "t_AFP"] .- df_res[:, "t_ES"]
    df_res[!, "lead_0"] .= (df_res[:, "diff"] .> 0) .| (df_res[:, "ES_only"] == true)
    df_res[!, "lead_30"] .= (df_res[:, "diff"] .> 30) .| (df_res[:, "ES_only"] == true)
    df_res[!, "lead_60"] .= (df_res[:, "diff"] .> 60) .| (df_res[:, "ES_only"] == true)
    
    df_det = DataFrame()
    for  uni in unique(df_res[:, col])
        dfM = filter(x -> x[col] == uni, df_res)
        detect_prob = mean(dfM[:, "Either"])
        ES_only = mean(dfM[:, "ES_only"])
        lead_0 = mean(dfM[!, "lead_0"])
        lead_30 = mean(dfM[!, "lead_30"])
        lead_60 = mean(dfM[!, "lead_60"])
        prob_new = (
            lead_0=(lead_0 + ES_only)/detect_prob, 
            lead_30=(lead_30 + ES_only)/detect_prob, 
            lead_60=(lead_60 + ES_only)/detect_prob, 
            uni=uni
            )
        push!(df_det, prob_new)
    end
    rename!(df_det, :uni => col)
    return df_det
end

function vis_early_detect_prob(df_det::DataFrame, col; kwargs...)
    y_max = maximum(df_det[:, "lead_0"])*100
    xticks=[0, 20, 40, 60, 80, 100]
    pl = plot(
        xlabel="ES population coverage (%)", 
        ylabel="Probability of \nearly detection by ES (%)",
        ylim=[0, 100], xticks=xticks,
        fmt=:png;
        kwargs...
        
    )
    plot!(pl, df_det[:, col], df_det[:, "lead_0"]*100, label="LT >0 days")
    plot!(pl, df_det[:, col], df_det[:, "lead_30"]*100, label="LT >30 days",
        marker=:circle, markersize=2, markerstrokewidth = 0,
    )
    plot!(pl, df_det[:, col], df_det[:, "lead_60"]*100, label="LT >60 days")
    return pl
end

function vis_ES_sensitivity_res(res_all::Tuple, n_sim::Int64, sp_pars::NamedTuple)
    df_res1, df_res2, df_res3 = res_all

    # ES population coverage
    per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
    sens_index = obtain_ES_sensitivity_index(sp_pars.pop, 0.01)

    # ES cachment area
    df_res = df_res3
    tab = detection_pattern_sensitivity(df_res, :ind_site)
    tab[:, :per_pop] = per_pop[sens_index]
    pl1 = vis_detection_pattern(tab, n_sim, :per_pop; 
        xlabel="Population coverage (%)", title="",
    )

    df_diff = leadtime_diff_sensitivity(df_res, :ind_site)
    df_diff[:, :per_pop] = per_pop[sens_index]
    pl2 = vis_leadtime_diff_sensitivity(df_diff, :per_pop;  
        xlabel="Population coverage (%)", title="",
        )
        df_det = early_detect_prob_sensitivity(df_res, :ind_site) 

    df_det = early_detect_prob_sensitivity(df_res, :ind_site) 
    df_det[:, :per_pop] = per_pop[sens_index]
    pl3 = vis_early_detect_prob(df_det, :per_pop,
        xlabel="Population coverage (%)", title="",
    )
    plot!(pl3, [0, 100.0], [0,100.0], color=:black, label=:none, linestyle=:dot)
    # Sampling frequency
    df_res = df_res2
    tab = detection_pattern_sensitivity(df_res, :n_freq)
    pl4 = vis_detection_pattern(
        tab, n_sim, :n_freq; 
        xlabel="Sampling frequency (days)",
        title="",
        )
    df_diff = leadtime_diff_sensitivity(df_res, :n_freq)
    pl5 = vis_leadtime_diff_sensitivity(
        df_diff, :n_freq; 
        xlabel="Sampling frequency (days)",
        title="",
        )
    df_det = early_detect_prob_sensitivity(df_res, :n_freq) 
    pl6 = vis_early_detect_prob(df_det, :n_freq,
        xlabel="Sampling frequency (days)",
    )

    # ES sensitivity, g
    df_res = df_res1
    tab = detection_pattern_sensitivity(df_res, :pop90)
    xticks = (1, 3, 10, 30, 100, 300, 1000)
    xlabel = "ES sensitivity, Np90"
    pl7 = vis_detection_pattern(
        tab, n_sim, :pop90; 
        xscale=:log10, xticks=(xticks, xticks), 
        xlabel=xlabel, title="",
    )
    df_diff = leadtime_diff_sensitivity(df_res, :pop90)
    pl8 = vis_leadtime_diff_sensitivity(
        df_diff, :pop90; 
        xscale=:log10, xticks=(xticks, xticks),
        xlabel=xlabel, title="",
    )
    df_det = early_detect_prob_sensitivity(df_res, :pop90) 
    pl9 = vis_early_detect_prob(df_det, :pop90,
        xscale=:log10, xticks=(xticks, xticks),
        xlabel=xlabel,
    )

    pls = [pl1, pl2, pl3, pl4, pl5, pl6, pl7, pl8, pl9]
    plot(pls..., fmt=:png, 
        size=(1200, 300*3),  # 400 * 2
        layout=(3,3),
        left_margin=5Plots.mm, bottom_margin=5Plots.mm)
end

function part_plot_detection_pattern(pl, df_res::DataFrame, col; label="", kargs...)
    tab = nothing
    if col == :per_pop
        tab = detection_pattern_sensitivity(df_res, :ind_site)
        tab[:, :per_pop] = per_pop[sens_index]
    else
        tab = detection_pattern_sensitivity(df_res, col)
    end
    det = tab[:, :Detect]
    x = tab[:, col]
    y = tab[:, :Both]./det.*100
    plot!(pl, x, y, 
        label=label,
        ylim=[0, 100],
        ylabel=ylabel_prob; 
        kargs...
        )
    y = tab[:, "ES only"]./det.*100
    plot!(pl, x, y, label=:none, linestyle=:dot; 
        kargs...
    )
end

function part_plot_diff_sensitivity(pl, df_res::DataFrame, col; label="",kargs...)
    df_diff = nothing
    if col == :per_pop
        df_diff = leadtime_diff_sensitivity(df_res, :ind_site)
        df_diff[:, :per_pop] = per_pop[sens_index]
    else
        df_diff = leadtime_diff_sensitivity(df_res, col)
    end
    plot!(pl, df_diff[:, col], df_diff[:, :q50],
        label=label,
        ylabel=ylabel_lead; 
        kargs...
    )
end

function part_plot_early_det(pl, df_res::DataFrame, col; label="", kargs...)
    df_det = nothing
    if col == :per_pop
        df_det = early_detect_prob_sensitivity(df_res, :ind_site)
        df_det[:, :per_pop] = per_pop[sens_index]
    else
        df_det= early_detect_prob_sensitivity(df_res, col)
    end
    y_max = maximum(df_det[:, :lead_30]).*100
    xticks = [0, 20, 40, 60, 80, 100]
    plot!(pl, df_det[:, col], df_det[:, :lead_30].*100,
        label=label,
        ylim=[0, 100.0],
        xticks=xticks,
        ylabel=ylabel_early; 
        kargs...
    )
end

"""
    Following functions are to visualise the heatmap including "only" patterns.
    - create_lead_time_category
    - proportion_each_cate_by_group
    - discretise_balance_color
    - visualise_heatmap_include_only
    - add_reverse_order_legend!
"""
function create_lead_time_category(df_res::DataFrame)::Tuple{DataFrame, Vector}
    # Replace NaN with Inf for t_AFP and t_ES
    df_res.t_AFP = ifelse.(isnan.(df_res.t_AFP), Inf, df_res.t_AFP)
    df_res.t_ES = ifelse.(isnan.(df_res.t_ES), Inf, df_res.t_ES)

    # Calculate the lead_time
    df_res.lead_time = df_res.t_AFP .- df_res.t_ES
    df_fil = filter(x -> !isnan.(x.lead_time), df_res)

    # Use "extend" option for cut, and assin Inf as "ES only".
    #bin_edges = [-Inf, -10_000, -150, -100, -50, 0, 50, 100, 150, 10_000, Inf]
    #bin_labels = ["AFP surv. only", 
    #    "<-150 LT", "-150 ~ -101 LT", "-100 ~ -49 LT", "-50 ~ -1 LT",
    #    "0 ~ 49 LT", "50 ~ 99 LT",
    #    "100 ~ 149 LT", ">=150 LT", 
    #    "ES only"]
    bin_edges = [-Inf, -10_000, -50, 0, 50, 10_000, Inf]
    bin_labels = ["AFP only", "<-50", "-50- -1", "0-49", ">50", "ES only"]
    df_fil.lead_time_category= cut(df_fil.lead_time, bin_edges, extend=true, labels=bin_labels)
    return (df_fil, bin_labels)
end

function proportion_each_cate_by_group(df_fil::DataFrame, col::Symbol)::Tuple{Matrix, DataFrame}
    # Group by "ind_site" and calculate the proportion
    grp_prop = DataFrames.combine(groupby(df_fil, [col, :lead_time_category]), nrow => :count)
    uni = grp_prop[:, col] |> unique
    grp_p = DataFrame()
    grp_50 = DataFrame()
    for i in uni
        dfM = filter(x -> x[col] == i, grp_prop)
        dfM[!, :prop_cate] = dfM[:, :count]/ sum(dfM[:,:count])*100
        grp_p = vcat(grp_p, dfM)

        dfM = filter(x -> x[col] == i, df_fil)
        cond = dfM[:, :lead_time] .< 0
        prop_0 = nrow(dfM[cond, :])/nrow(dfM)*100
        grp_50 = vcat(grp_50, DataFrame(col=>i, :prop_0=>prop_0))
    end
    grp_p = sort(grp_p, col)
    grp_50 = sort(grp_50, col)
    df_ind_per = DataFrame(col => df_fil[:, col] |> unique )
    #grp_p = leftjoin(grp_p, df_ind_per, on=col)
    grp_p_unstack = unstack(grp_p, col, :lead_time_category, :prop_cate) 
    y = @pipe grp_p_unstack[:, bin_labels] |> Matrix .|> coalesce(_, 0.0)
    return (y, grp_50)
end

function discretise_balance_color(bin_labels::Vector)
    colors = cgrad(:balance, length(bin_labels), categorical = true, rev=true)
    colors = reshape([colors[i] for i in 1:length(colors)], 1,:)
    colors
end

function visualise_heatmap_include_only(x, y, grp_50, bin_labels)
    colors = discretise_balance_color(bin_labels)
    pl = areaplot(x, y, size=(800,600), label=nothing, fillalpha = [0.85 0.85], color=colors)
    plot!(x, grp_50[:, :prop_0], color="black", label=nothing, lw=1.5, 
        marker=:circle, markersize=3, markerstrokewidth = 0,
    )
    #ls=:dashdot)
    pl
end

function add_reverse_order_legend!(pl::Plots.Plot, bin_labels::Vector)
    # for color adjusting.
    colors = discretise_balance_color(bin_labels)
    areaplot!(pl, 1, reshape([0 for i in 1:length(bin_labels)], 1, :),
        label=reshape(reverse(bin_labels), 1, :), 
        color=reverse(colors)) 
end

function df_to_heatmap(df_res::DataFrame, x, col; 
        xticks=:none, add_zero::Bool=false, pc=1.0
        )
    df_fil, bin_labels = create_lead_time_category(df_res)
    colors = discretise_balance_color(bin_labels)
    y, grp_50 = proportion_each_cate_by_group(df_fil, col) 
    if add_zero == true
        n_cate = size(y)[2]
        y0 = hcat([100], fill(0, 1, n_cate-1))
        y = vcat(y0, y)
        grp_50 = vcat(DataFrame(Dict(col=>0, :prop_0 => 100)), grp_50)
        x = vcat([0], x)     
    end

    pl = visualise_heatmap_include_only(x, y, grp_50, bin_labels)
    ticks= [0, 20, 40 ,60, 80, 100]
    hline!(pl, ticks, color=:white, alpha=0.5, ls=:dash, label=:none)
    if (col == :ind_site)  & (maximum(x) <= 100)
        ticks_scaled = Int64.(pc.*ticks)
        plot!(pl, xticks=(ticks, ticks_scaled), yticks=(ticks, ticks), 
            xlim=[0,100], ylim=[0,100],
        )
        vline!(pl, ticks, color=:white, alpha=0.5, ls=:dash, label=:none) 
    end
    if xticks != :none
        plot!(pl, xticks=(xticks, xticks), yticks=(ticks, ticks),
            ylim=[0, 100],
        )
        vline!(pl, xticks, color=:white, alpha=0.5, ls=:dash, label=:none) 
    end
    return pl
end

@doc (@doc create_lead_time_category) proportion_each_cate_by_group, 
    discretise_balance_color, 
    visualise_heatmap_include_only, 
    create_reverse_order_legend!,
    df_to_heatmap

