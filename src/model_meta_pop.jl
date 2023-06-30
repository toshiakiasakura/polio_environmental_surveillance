using Accessors
using Parameters
using ProtoStructs
using Random
using StatsBase

include("util.jl")

struct SEIRMetaModelParams
    R0::Float64 
    γ1::Float64 
    γ2::Float64
    γ3::Float64 
    P_AFP::Float64 
    P_H::Float64 
    P_AFP_sample::Float64 
    P_AFP_test::Float64 
    P_ES_test::Float64
    N0::Vector{Int64}
    I0_init::Int64 
    I0_ind::Int64
    λ0::Float64 
    days::Int64
    ES_n_freq::Int64
    n_site::Int64
    ES_area::Vector{Bool}
    α::Float64
    π_mat::Matrix{Float64}
    β::Float64
    function SEIRMetaModelParams(;
            R0=10.0, γ1=1/4.0, γ2=1/15.02, γ3=1/7.0,
            P_AFP=1/200, P_H=0.9, 
            P_AFP_sample=0.8, P_AFP_test=0.97, P_ES_test=0.97,
            N0=[1], I0_init=1, I0_ind=1, λ0=0.5, 
            days=365, ES_n_freq=30, n_site=1, ES_area=[1],
            α=0.05, π_mat=fill(0,1,1)
        )
        new(R0, γ1, γ2, γ3, P_AFP, P_H, 
            P_AFP_sample, P_AFP_test, P_ES_test,
            N0, I0_init, I0_ind, 
            λ0, days, ES_n_freq, n_site, ES_area, 
            α, π_mat, R0*γ2
            )
    end
end

@proto struct SEIRMetaModelRecord
    days::Int64
    n_site::Int64
    S::Array{Int64, 2} = fill(-999, n_site, days)
    E::Array{Int64, 2} = fill(-999, n_site, days)
    Ia::Array{Int64, 2} = fill(-999, n_site, days)
    I_AFP::Array{Int64, 2} = fill(-999, n_site, days)
    R::Array{Int64, 2} = fill(-999, n_site, days)
end

@proto struct SEIRMetaModelOneStep
    S::Vector{Int64}
    E::Vector{Int64}
    Ia::Vector{Int64}
    I_AFP::Vector{Int64}
    R::Vector{Int64}
    H_AFP::Vector{Int64}
end

function set_values!(
        rec::SEIRMetaModelRecord, 
        model::SEIRMetaModelOneStep, 
        ind::Int64
    )
    @unpack S, E, Ia, I_AFP, R = model
    rec.S[ind] = S
    rec.E[ind] = E
    rec.Ia[ind] = Ia
    rec.I_AFP[ind] = I_AFP
    rec.R[ind] = R
end

function initialize_model(params::SEIRMetaModelParams, rec_flag::Bool)
    @unpack N0, I0_init, I0_ind, n_site, days = params
    S = copy(N0)
    S[I0_ind] -= I0_init
    Ia = fill(0, n_site)
    Ia[I0_ind] -= I0_init

    model = SEIRMetaModelOneStep(
        S = N0 - I0_init,
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

    @unpack γ1, γ2, γ3, P_AFP, days, β, N0 = params
    @unpack S, E, Ia, I_AFP, R = model

    p_inf = β/N0*(Ia + I_AFP)
    new_inf = rand_binom(S, 1 .- exp(-p_inf))
    rec_E = rand_binom(E, 1 .- exp(-γ1))
    rec_E_Ia = rand_binom(rec_E, 1 - P_AFP)
    rec_E_I_AFP = rec_E - rec_E_Ia

    rec_Ia = rand_binom(Ia, 1 .- exp(-γ2))
    rec_I_AFP = rand_binom(I_AFP, 1 .- exp(-γ3))

    S += - new_inf
    E += + new_inf - rec_E
    Ia +=          + rec_E_Ia    - rec_Ia
    I_AFP +=       + rec_E_I_AFP - rec_I_AFP
    R +=                         + rec_Ia + rec_I_AFP
    model = SEIRMetaModelOneStep(
        S=S, E=E, Ia=Ia, I_AFP=I_AFP,R=R, 
        H_AFP=rec_I_AFP
    )
    return model
end

function run_sim(params::SEIRMetaModelParams; rec_flag::Bool=false)
    @unpack γ1, γ2, γ3, P_AFP, P_H, P_AFP_sample, 
            P_AFP_test, P_ES_test, days, ES_n_freq, λ0 = params

    rec, model = initialize_model(params, rec_flag)
    # set ES scheudle
    nt = [0 for _ in 1:days]
    st_ind = rand(1:ES_n_freq)
    nt[st_ind:30:days] .= 1

    t_extinct = NaN
    t_AFP = NaN
    t_ES = NaN
    R_final = NaN
    
    for t in 2:days
        # Stop condition
        if (model.E + model.Ia + model.I_AFP) == 0
            t_extinct = isnan(t_extinct) == true ? t : t_extinct
            R_final = sum(model.R)
            if rec_flag == true
                set_values!(rec, model, t)
                continue
            else
                break
            end
        end

        model = update_model(model, params)
        if rec_flag == true
            set_values!(rec, model, t)
        end

        # ES surveillance
        if isnan(t_ES) == true
            ωt = 1 .- exp.(-λ0 .* (model.Ia .+ model.I_AFP))
            wt = rand_binom.(nt[t], ωt.*P_ES_test)
            t_ES = sum(wt) == 1 ? t : t_ES
        end
        # AFP surveillance
        if isnan(t_AFP) == true
            Ct = rand_binom(model.H_AFP, P_H*P_AFP_sample*P_AFP_test)
            t_AFP = Ct >= 1 ? t : t_AFP
        end

    end
    outcome = (t_ES=t_ES, t_AFP=t_AFP, t_extinct=t_extinct, R_final=R_final)
    return (rec, outcome)
end

