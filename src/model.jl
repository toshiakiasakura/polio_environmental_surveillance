using Accessors
using Parameters
using ProtoStructs
using Random
using StatsBase

include("util.jl")

struct SEIRModelParams
    R0::Float64 
    γ1::Float64 
    γ2::Float64
    σ::Float64 
    P_AFP::Float64 
    P_H::Float64 
    P_AFP_sample::Float64 
    P_AFP_test::Float64 
    P_ES_test::Float64
    N_tot::Int64
    N_unvac::Int64
    I0_init::Int64 
    g::Float64 
    days::Int64
    ES_n_freq::Int64
    μ::Float64
    β::Float64
    function SEIRModelParams(;
            R0=10.0, γ1=1/4.0, γ2=1/15.02, σ=0.329,
            P_AFP=1/200, P_H=0.9, 
            P_AFP_sample=0.8, P_AFP_test=0.97, P_ES_test=0.97,
            N_tot=10000, N_unvac=1000, I0_init=1, g=0.23, 
            ES_n_freq=30, days=365, μ=1/(5*365)
        )
        new(R0, γ1, γ2, σ, P_AFP, P_H, 
            P_AFP_sample, P_AFP_test, P_ES_test,
            N_tot, N_unvac, I0_init, g, days, ES_n_freq, μ, R0*γ2
            )
    end
end

@proto struct SEIRModelRecord
    days::Int64
    S::Array{Int64, 1} = fill(-999, days)
    E::Array{Int64, 1} = fill(-999, days)
    I::Array{Int64, 1} = fill(-999, days)
    R::Array{Int64, 1} = fill(-999, days)
    A::Array{Int64, 2} = fill(-999, days, 6)
    Z_A5_6::Array{Int64, 1} = fill(-999, days)
end

@proto struct SEIRModelOneStep
    S::Int64
    E::Int64
    I::Int64
    R::Int64
    A::Array{Int64, 1} 
    Z_A5_6::Int64
end

function set_values!(
        rec::SEIRModelRecord, 
        model::SEIRModelOneStep, 
        ind::Int64
    )
    @unpack S, E, I, R, A, Z_A5_6 = model
    rec.S[ind] = S
    rec.E[ind] = E
    rec.I[ind] = I
    rec.R[ind] = R
    rec.A[ind, :] = A
    rec.Z_A5_6[ind] = Z_A5_6
end

function initialize_model(params::SEIRModelParams, rec_flag::Bool)
    @unpack N_tot, N_unvac, I0_init, days = params
    model = SEIRModelOneStep(
        S = N_unvac - I0_init,
        E = 0,
        I = I0_init,
        R = 0,
        A = [0, 0, 0, 0, 0, 0],
        Z_A5_6 = 0, # Newly seeking
    )

    if rec_flag == true
        rec = SEIRModelRecord(days=days)
        set_values!(rec, model, 1)
    else
        rec = SEIRModelRecord(days=1)
    end
    return rec, model
end

function update_model(
        model::SEIRModelOneStep, 
        params::SEIRModelParams
    )::SEIRModelOneStep

    @unpack γ1, γ2, σ, P_AFP, days, β, N_tot, μ = params
    @unpack S, E, I, R, A = model

    new_born = rand_binom(N_unvac, μ)
    p_inf = β/N_tot*I
    rem_S = rand_binom(S, 1 .- exp(- p_inf- μ))
    new_E = rand_binom(rem_S, p_inf/(p_inf + μ))

    rem_E = rand_binom(E, 1 .- exp(-γ1 - μ))
    new_I = rand_binom(rem_E, γ1/(γ1 + μ))

    rem_I = rand_binom(I, 1 .- exp(- γ2 - μ))
    new_R = rand_binom(rem_I, γ2/(γ2 + μ))

    new_AFP1 = rand_binom(new_E, P_AFP)
    new_AFP = [new_AFP1]
    for i in 1:5
        new_A_stage = rand_binom(A[i], 1 - exp(-σ)) # Ai to next, including H.
        push!(new_AFP, new_A_stage)
    end

    S += new_born - rem_S
    E +=          + new_E - rem_E
    I +=                  + new_I - rem_I
    R +=                          + new_R

    for i in 1:5
        A[i] += new_AFP[i] - new_AFP[i+1]
    end
    A[6] += new_AFP[6]

    model = SEIRModelOneStep(
        S=S, E=E, I=I, R=R, A=A,
        Z_A5_6=new_AFP[6]
    )
    return model
end

function run_sim(params::SEIRModelParams; rec_flag::Bool=false)
    @unpack γ1, γ2, σ, P_AFP, P_H, P_AFP_sample, 
            P_AFP_test, P_ES_test, days, ES_n_freq, g = params

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
        if (model.E + model.I + sum(model.A)) == 0
            t_extinct = isnan(t_extinct) == true ? t : t_extinct
            R_final = model.R
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
            ωt = 1 - exp(-g * model.I)
            wt = rand_binom(nt[t], ωt*P_ES_test)
            t_ES = wt == 1 ? t : t_ES
        end
        # AFP surveillance
        if isnan(t_AFP) == true
            Ct = rand_binom(model.Z_A5_6, P_H*P_AFP_sample*P_AFP_test)
            t_AFP = Ct >= 1 ? t : t_AFP
        end

    end
    R_final_AFP = model.A[6]
    outcome = (t_ES=t_ES, t_AFP=t_AFP, t_extinct=t_extinct, 
        R_final=R_final, R_final_AFP=R_final_AFP,
        )
    return (rec, outcome)
end

