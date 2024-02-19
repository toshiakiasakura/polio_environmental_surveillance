include("utils.jl")

Base.@kwdef mutable struct SEIRMetaModelParams
    R0::Float64 = 10.0
    γ1::Float64 = 1/4.0
    γ2::Float64 = 1/15.02
    σ::Float64 = 0.3291
    P_AFP::Float64 = 1/200
    N_tot::Vector{Int64}
    N_unvac::Vector{Int64}
    I0_init::Int64 = 1
    I0_ind::Int64 = 1
    days::Int64 = 365*3
    n_site::Int64
    α::Float64 = 0.05
    π_mat::Matrix{Float64}
    μ::Float64 = 1/(5*365)
    β::Float64 = R0*γ2
    pc::Float64 = 0.25
    imp_ws::Vector{Float64}
end

Base.@kwdef mutable struct ESParams
    ES_μ::Float64 = 1.308
    ES_σ::Float64 = 1.493
    Pop_whole::Vector{Float64}
    P_test::Float64 = 0.97
    area::Vector{Bool}
    n_freq::Int64 = 30
    pc::Float64
end

Base.copy(x::ESParams) = ESParams(
    ES_μ=x.ES_μ, ES_σ=x.ES_σ, Pop_whole=x.Pop_whole,
    P_test=x.P_test, area=x.area, n_freq=x.n_freq,
    pc=x.pc,
)

Base.@kwdef mutable struct AFPSurParams
    P_test::Float64 = 0.97
    P_sample::Float64 = 0.53
    P_H::Float64 = 0.9
end

Base.@kwdef mutable struct SEIRMetaModelOneStep
    S::Vector{Int64}
    E::Vector{Int64}
    Ic::Vector{Int64}
    Inc::Vector{Int64}
    R::Vector{Int64}
    A::Array{Int64, 2}
    Z_A5_6::Int64
    Z_S_E::Int64
end

Base.@kwdef mutable struct SEIRMetaModelRecord
    days::Int64
    n_site::Int64
    Ic::Array{Int64, 2} = fill(0, n_site, days)
    Z_A5_6::Array{Int64, 1} = fill(0, days)
    Z_S_E::Array{Int64, 1} = fill(0, days)
end

Base.@kwdef mutable struct SEIRMetaModelRecordSparse
    days::Int64
    n_site::Int64
    Ic::SparseMatrixCSC{Int64, Int64}
    Z_A5_6::SparseVector{Int64, Int64}
    Z_S_E::SparseVector{Int64, Int64}
end

function set_values!(
        rec::SEIRMetaModelRecord,
        model::SEIRMetaModelOneStep,
        ind::Int64
    )
    @unpack S, E, Ic, Inc, Z_A5_6, Z_S_E = model
    rec.Ic[:, ind] = Ic
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
        I0_ind = wsample(1:n_site, imp_ws)
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
    #new_E = rand_beta_binom.(S, 1 .- exp.(-λ), 1/5_000)
    #rem_S = rand_binom.(S .- new_E, 1 .- exp.(-μ)) .+ new_E

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
        Ic=sparse(rec.Ic),
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

function cdf_lognormal(lognorm::LogNormal, v)
    if v == 0
        return 0
    else
        return cdf(lognorm, v)
    end
end

function enviro_surveillance(rec::SEIRMetaModelRecordSparse, par_ES::ESParams)::Float64
    @unpack ES_μ, ES_σ, P_test, Pop_whole, area, n_freq, pc = par_ES
    _, days = size(rec.Ic)
    lognorm = LogNormal(ES_μ, ES_σ)
    n_area = sum(area)

    nt_area = fill(0, n_area, days)
    st_ind = rand(1:n_freq)
    sample_ind = st_ind:n_freq:days
    nt_area[:, sample_ind] .= 1
    Ic_area = rec.Ic[area, :]
    Pop_whole_area = Pop_whole[area] .* pc

    t_ES = NaN
    for t in sample_ind
        # ES is covered from the first, so if 1 is 0, this is not appropriate.
        if nt_area[1, t] == 0
            continue
        end
        ωt = cdf_lognormal.(lognorm, Ic_area[:, t] ./Pop_whole_area.*100_000)
        wt = rand_binom.(nt_area[:,t], ωt .* P_test)
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
        jldsave("$(path)/$(i).jld2"; data=res)
    end
    return path
end

function fetch_sim_paths(path::String)
    return glob("$(path)/*.jld2")
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

"""Run transmission model part given the parameters.
"""
function run_transmission_model(;
        R0=14.0, α=0.05, pc=0.25, imp_ws=[1.0],
        n_sim=10, pattern="", ES_pattern="ES_population_size"
    )
    sp_pars = read_spatial_params_file(ES_pattern)
    n_site = length(sp_pars.pop)
    println("# of sites: $(n_site)")
    pars = SEIRMetaModelParams(
        R0=R0, α=α, pc=pc,
        N_tot=sp_pars.pop,
        N_unvac=sp_pars.unvac,
        π_mat=sp_pars.π_mat,
        n_site=n_site,
        imp_ws=imp_ws,
    )

    pars |> dump
    Random.seed!(48)
    path_trans = run_and_save_sim(pars; n_sim=n_sim, pattern=pattern)
    return path_trans
end

"""Set AFP and ES parameters.
"""
function set_par_AFP_ES(;
        pattern="", pc=1.0, n_freq=30,
        ES_μ=1.31, ES_σ=1.49, ES_pattern="ES_population_size"
    )
    sp_pars = read_spatial_params_file(ES_pattern)
    n_site = length(sp_pars.pop)

    # Set AFP params.
    par_AFP = AFPSurParams()

    # Set ES params. The value comes from the empirical one.
    ES_obs_cov = 0.0859
    # ES coverage rate for the baseline.
    cov_rate = ES_obs_cov/pc

    # Catchment area
    pop = sp_pars.pop
    n_site = length(pop)
    cum_prop = cumsum(pop/sum(pop))
    cov = abs.(cum_prop .- cov_rate) |> argmin
    println("Cum index:", cov)
    println("Actual ES coverage: ", sum(pop[1:cov])/sum(pop)*100)
    area = fill(0, n_site)
    area[1:cov] .= 1

    # Set the whole population.
    Pop_whole = sp_pars.df[:, :value_whole]

    par_ES = ESParams(area=area, n_freq=n_freq,
        ES_μ=ES_μ, ES_σ=ES_σ, Pop_whole=Pop_whole,
        pc=pc
    )
    dump(par_AFP)
    dump(par_ES)
    return par_AFP, par_ES
end

"""Obtain a vector of index each of which detemines the cumulative
ES coverage. `inc_prop` determines the unit of increment in
the ES coverage.
"""
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

"""Sensitivity analysis for ES catchment area.

Args:
- `res`: a single meta-population model simulation result.
"""
function sensitivity_ES_catchment_area(
        res::NamedTuple, par_AFP::AFPSurParams,
        par_ES::ESParams, sim_res::DataFrame;
        inc_prop=0.01
    )
    rec, outcome, pars = res
    t_AFP = AFP_surveillance(rec, par_AFP)

    sens_index = obtain_ES_sensitivity_index(res.pars.N_tot, inc_prop) # 1% increment.
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

function save_sensitivity_ES_catchment_area(
        par_AFP, par_ES, path_trans;
        inc_prop=0.01, pattern="", ES_pattern="",
    )
    path_objs = fetch_sim_paths(path_trans)
    sim_res = DataFrame()
    @showprogress for i in 1:length(path_objs)
        res = load(path_objs[i])["data"]
        par_ES_tmp = copy(par_ES)
        sim_res = sensitivity_ES_catchment_area(
            res, par_AFP, par_ES_tmp, sim_res;
            inc_prop=inc_prop
        )
    end
    # Fetch parameters used in the transmission model..
    paths = fetch_sim_paths(path_trans)
    res = load(paths[1])["data"]

    now_str = get_today_time()
    path_save = "../dt_tmp_res/sens_ES_catchment_$(now_str).jld2"
    jldsave(path_save,
        sim_res=sim_res, inc_prop=inc_prop, path_trans=path_trans,
        pars=res.pars, par_ES=par_ES,
        pattern=pattern, ES_pattern=ES_pattern
    )
    return path_save
end

#"""Function for `sensitivity_ana_ES`
#"""
#function sensitivity_sampling_frequency(
#        res::NamedTuple, par_AFP::AFPSurParams,
#        par_ES::ESParams, sim_res::DataFrame,
#    )
#    rec, outcome, pars = res
#    t_AFP = AFP_surveillance(rec, par_AFP)
#
#    freq_list = [
#        1, 7, 14, 21, 28, 30, 35, 42, 49, 56,
#    ]
#    for n_freq in freq_list
#        par_ES.n_freq = n_freq
#        t_ES = enviro_surveillance(rec, par_ES)
#        R_inf_ES = isnan(t_ES) == false ? sum(rec.Z_S_E[begin:Int64(t_ES)]) : NaN
#        res_t = (
#            t_AFP=t_AFP,
#            t_ES=t_ES,
#            R_inf_ES=R_inf_ES,
#            n_freq=n_freq,
#            outcome...
#        )
#        push!(sim_res, res_t, promote=true)
#    end
#    return sim_res
#end
#
#function save_sensitivity_sampling_frequency(
#        par_AFP, par_ES, path_trans
#    )
#    path_objs = fetch_sim_paths(path_trans)
#    sim_res = DataFrame()
#    @showprogress for i in 1:length(path_objs)
#        res = load(path_objs[i])["data"]
#        par_ES_tmp = copy(par_ES)
#        sim_res = sensitivity_sampling_frequency(
#            res, par_AFP, par_ES_tmp, sim_res
#        )
#    end
#    now_str = get_today_time()
#    path_save = "../dt_tmp_res/sens_sampling_frequency_$(now_str).jld2"
#    jldsave(path_save, sim_res=sim_res, path_trans=path_trans)
#    return path_save
#end

function remove_all_transmission_results(path_trans)
    rm(path_trans; recursive=true)
end

function remove_all_leaving_one(path_trans)
    path_objs = fetch_sim_paths(path_trans)
    for p in path_objs
        if split(p, "/")[end] == "1.jld2"
            continue
        end
        rm(p)
    end
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


