include("utils.jl")

Base.@kwdef mutable struct SEIRMetaModelParams
    R0::Float64 = 14.0
    γ1::Float64 = 1 / 4.0
    γ2::Float64 = 1 / 15.02
    σ::Float64 = 0.3291
    P_AFP::Float64 = 1 / 200
    Nc::Vector{Int64}
    N_unvac::Vector{Int64}
    I0_init::Int64 = 1
    # I0_ind::Int64 = 1
    days::Int64 = 365 * 3
    n_site::Int64
    α::Float64 = 0.05
    π_mat::Matrix{Float64}
    μ::Float64 = 1 / (5 * 365)
    β::Float64 = R0 * γ2
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
    A::Array{Int64,2}
    Z_A5_6::Int64
    Z_S_E::Int64
end

Base.@kwdef mutable struct SEIRMetaModelRecord
    days::Int64
    n_site::Int64
    Ic::Array{Int64,2} = fill(0, n_site, days)
    Z_A5_6::Array{Int64,1} = fill(0, days)
    Z_S_E::Array{Int64,1} = fill(0, days)
end

Base.@kwdef mutable struct SEIRMetaModelRecordSparse
    days::Int64
    n_site::Int64
    Ic::SparseMatrixCSC{Int64,Int64}
    Z_A5_6::SparseVector{Int64,Int64}
    Z_S_E::SparseVector{Int64,Int64}
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
    @unpack I0_init, n_site, Nc, N_unvac, days, pc, imp_ws = params
    if pattern == "population_size"
        I0_ind = wsample(1:n_site, Nc)
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
        S=S,
        E=fill(0, n_site),
        Ic=Ic,
        Inc=Inc,
        R=fill(0, n_site),
        A=fill(0, n_site, 6),
        Z_A5_6=0,
        Z_S_E=0,
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

    @unpack γ1, γ2, σ, μ, P_AFP, days, β, Nc, N_unvac, π_mat, α, n_site, pc = params
    @unpack S, E, Ic, Inc, R, A = model

    new_born = rand_binom.(N_unvac, 1 - exp(-μ))
    ext = sum(π_mat' .* (Ic .+ Inc)', dims=2)[:, 1]
    λ = β ./ Nc .* ((1 - α) .* (Ic + Inc) .+ α .* ext)

    rem_S = rand_binom.(S, 1 .- exp.(-λ .- μ))
    new_E = rand_binom.(rem_S, λ ./ (λ .+ μ))

    rem_E = rand_binom.(E, 1 .- exp(-γ1 - μ))
    new_I = rand_binom.(rem_E, γ1 / (γ1 + μ))
    new_Ic = rand_binom.(new_I, pc)
    new_Inc = new_I .- new_Ic

    rem_Ic = rand_binom.(Ic, 1 .- exp(-γ2 - μ))
    rem_Inc = rand_binom.(Inc, 1 .- exp(-γ2 - μ))
    new_R = rand_binom.(rem_Ic .+ rem_Inc, γ2 / (γ2 + μ))

    new_AFP = fill(0, n_site, 6)
    new_AFP[:, 1] = rand_binom.(new_E, P_AFP)
    rate = 1 - exp(-σ)
    for i in 1:5
        new_AFP[:, i+1] .= rand_binom.(A[:, i], rate)
    end

    S .+= new_born .- rem_S
    E .+= .+new_E .- rem_E
    Ic .+= .+new_Ic .- rem_Ic
    Inc .+= .+new_Inc .- rem_Inc
    R .+= .+new_R

    for i in 1:5
        A[:, i] .+= new_AFP[:, i] .- new_AFP[:, i+1]
    end
    A[:, 6] .+= new_AFP[:, 6]

    model = SEIRMetaModelOneStep(
        S=S, E=E, Ic=Ic, Inc=Inc, R=R, A=A,
        Z_A5_6=sum(new_AFP[:, 6]), Z_S_E=sum(new_E) .|> Int64,
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
        w = rand_binom(I_tot[t], P_test * P_sample * P_H)
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
    # Remove rows without ES surveillance.
    Ic_area = rec.Ic[area, :]
    Pop_whole_area = Pop_whole[area] .* pc

    t_ES = NaN
    for t in sample_ind
        # ES is covered from the first, so if 1 is 0, this is not appropriate.
        if nt_area[1, t] == 0
            continue
        end
        ωt = cdf_lognormal.(lognorm, Ic_area[:, t] ./ Pop_whole_area .* 100_000)
        wt = rand_binom.(nt_area[:, t], ωt .* P_test)
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
    heatmap!(pl, 1:days, 1:n_site, log10.(I .+ 1))
    return pl
end

function run_and_save_sim(pars::SEIRMetaModelParams; n_sim=10, pattern="")
    now_str = get_today_time()
    path = "../dt_tmp/$now_str"
    mkdir(path)
    #println(path)
    @showprogress for i in 1:n_sim
        res = run_sim(pars; rec_flag=true, pattern=pattern)
        jldsave("$(path)/$(i).jld2"; data=res)
    end
    return path
end

function fetch_sim_paths(path::String)
    return glob("$(path)/*.jld2")
end

"""Run transmission model part given the parameters.
"""
function run_transmission_model(;
    R0=14.0, α=0.05, pc=0.25, imp_ws=[1.0],
    n_sim=10, pattern="", ES_pattern="ES_population_size"
)
    sp_pars = read_spatial_params_file(ES_pattern)
    n_site = length(sp_pars.pop)
    pars = SEIRMetaModelParams(
        R0=R0, α=α, pc=pc,
        Nc=sp_pars.pop,
        N_unvac=sp_pars.unvac,
        π_mat=sp_pars.π_mat,
        n_site=n_site,
        imp_ws=imp_ws,
    )

    #println("# of sites: $(n_site)")
    #pars |> dump
    Random.seed!(48)
    path_trans = run_and_save_sim(pars; n_sim=n_sim, pattern=pattern)
    return path_trans
end

"""Set AFP and ES parameters.
"""
function set_par_AFP_ES(;
    pattern, pc, ES_pattern,
    n_freq=30, ES_μ=1.31, ES_σ=1.49,
)
    sp_pars = read_spatial_params_file(ES_pattern)
    n_site = length(sp_pars.pop)

    # Set AFP params.
    par_AFP = AFPSurParams()

    # Set ES params. The value comes from the empirical one.
    ES_obs_cov = 0.0859
    # ES coverage rate for the baseline.
    cov_rate = ES_obs_cov / pc

    # Set the baseline catchment area (but not used anymore).
    pop = sp_pars.pop
    n_site = length(pop)
    cum_prop = cumsum(pop / sum(pop))
    cov = abs.(cum_prop .- cov_rate) |> argmin
    #println("Cum index:", cov)
    #println("Actual ES coverage: ", sum(pop[1:cov]) / sum(pop) * 100)
    area = fill(0, n_site)
    area[1:cov] .= 1

    # Set the whole population.
    Pop_whole = sp_pars.df[:, :value_whole]

    par_ES = ESParams(area=area, n_freq=n_freq,
        ES_μ=ES_μ, ES_σ=ES_σ, Pop_whole=Pop_whole,
        pc=pc
    )
    #dump(par_AFP)
    #dump(par_ES)
    return par_AFP, par_ES
end

"""Obtain a vector of index each of which detemines the cumulative
ES coverage. `inc_prop` determines the unit of increment in
the ES coverage.
"""
function obtain_ES_sensitivity_index(pop::Vector, inc_prop::Float64)::Vector
    n_site = length(pop)
    cum_prop = cumsum(pop / sum(pop))
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

    sens_index = obtain_ES_sensitivity_index(res.pars.Nc, inc_prop) # 1% increment.
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

"""Save a sensitivity analysis result.
"""
function save_sensitivity_ES_catchment_area(
    par_AFP, par_ES, path_trans;
    inc_prop=0.01, pattern, ES_pattern,
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

function remove_all_transmission_results(path_trans)
    rm(path_trans; recursive=true)
end

function post_processing(path_trans, pattern, ES_pattern, pc)
    par_AFP, par_ES = set_par_AFP_ES(;
        pc=pc, pattern=pattern, ES_pattern=ES_pattern
    )
    path_save = save_sensitivity_ES_catchment_area(
        par_AFP, par_ES, path_trans;
        inc_prop=0.01,
        pattern=pattern, ES_pattern=ES_pattern,
    )
    remove_all_transmission_results(path_trans)
    return path_save
end

"""
    run_trans_detection_process

Run the transmission model, AFP surveillance model,
and ES surveillance model.
"""
function run_trans_detection_process(;
    R0, α, pc, n_sim, pattern, ES_pattern,
)::NamedTuple
    # For pattern == "population_size", imp_ws is passed but not used.
    imp_ws = read_imp_ws_data(pattern, ES_pattern)

    path_trans = run_transmission_model(;
        R0=R0, α=α, pc=pc, imp_ws=imp_ws,
        n_sim=n_sim,
        pattern=pattern, ES_pattern=ES_pattern
    )
    path_save = post_processing(path_trans, pattern, ES_pattern, pc)
    println(path_save)
    rec = (R0=R0, α=α, pc=pc, n_sim=n_sim,
        pattern=pattern, ES_pattern=ES_pattern,
        path=path_save)
    return (rec)
end



