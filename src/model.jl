include("utils.jl")

Base.@kwdef mutable struct SEIRModelParams
    R0::Float64 = 14.0
    γ1::Float64 = 1 / 4.0
    γ2::Float64 = 1 / 15.02
    σ::Float64 = 0.329
    P_AFP::Float64 = 1 / 200
    P_H::Float64 = 0.9
    P_AFP_sample::Float64 = 0.53
    P_AFP_test::Float64 = 0.97
    P_ES_test::Float64 = 0.97
    Nc::Int64 = 10_000
    N_unvac::Int64
    I0_init::Int64 = 1
    days::Int64 = 365 * 3
    ES_n_freq::Int64 = 30
    μ::Float64 = 1 / (5 * 365)
    β::Float64 = R0 * γ2
    pc::Float64
    Pop_whole::Float64
    ES_μ::Float64 = 1.308
    ES_σ::Float64 = 1.493
end

Base.@kwdef mutable struct SEIRModelRecord
    days::Int64
    S::Array{Int64,1} = fill(-999, days)
    E::Array{Int64,1} = fill(-999, days)
    Ic::Array{Int64,1} = fill(-999, days)
    Inc::Array{Int64,1} = fill(-999, days)
    R::Array{Int64,1} = fill(-999, days)
    A::Array{Int64,2} = fill(-999, days, 6)
    Z_A5_6::Array{Int64,1} = fill(-999, days)
end

Base.@kwdef mutable struct SEIRModelOneStep
    S::Int64
    E::Int64
    Inc::Int64
    Ic::Int64
    R::Int64
    A::Array{Int64,1}
    Z_A5_6::Int64
end

function set_values!(
    rec::SEIRModelRecord,
    model::SEIRModelOneStep,
    ind::Int64
)
    @unpack S, E, Inc, Ic, R, A, Z_A5_6 = model
    rec.S[ind] = S
    rec.E[ind] = E
    rec.Ic[ind] = Ic
    rec.Inc[ind] = Inc
    rec.R[ind] = R
    rec.A[ind, :] = A
    rec.Z_A5_6[ind] = Z_A5_6
end

function initialize_model(params::SEIRModelParams, rec_flag::Bool; pattern="")
    @unpack Nc, N_unvac, I0_init, days, pc = params
    A = pattern == "check_incubation" ? [1000, 0, 0, 0, 0, 0] : [0, 0, 0, 0, 0, 0]

    Ic_init = rand_binom(I0_init, pc)
    Inc_init = I0_init - Ic_init
    model = SEIRModelOneStep(
        S=N_unvac - I0_init,
        E=0,
        Ic=Ic_init,
        Inc=Inc_init,
        R=0,
        A=A,
        Z_A5_6=0, # Newly seeking
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

    @unpack γ1, γ2, σ, P_AFP, days, β, Nc, N_unvac, μ, pc = params
    @unpack S, E, Ic, Inc, R, A = model

    new_born = rand_binom(N_unvac, 1 - exp(-μ))
    p_inf = β / Nc * (Ic + Inc)
    rem_S = rand_binom(S, 1 .- exp(-p_inf - μ))
    new_E = rand_binom(rem_S, p_inf / (p_inf + μ))

    rem_E = rand_binom(E, 1 .- exp(-γ1 - μ))
    new_I = rand_binom(rem_E, γ1 / (γ1 + μ))
    new_Ic = rand_binom(new_I, pc)
    new_Inc = new_I - new_Ic

    rem_Ic = rand_binom(Ic, 1 .- exp(-γ2 - μ))
    rem_Inc = rand_binom(Inc, 1 .- exp(-γ2 - μ))
    new_R = rand_binom(rem_Ic + rem_Inc, γ2 / (γ2 + μ))

    new_AFP1 = rand_binom(new_E, P_AFP)
    new_AFP = [new_AFP1]
    for i in 1:5
        new_A_stage = rand_binom(A[i], 1 - exp(-σ)) # Ai to next, including H.
        push!(new_AFP, new_A_stage)
    end

    S += new_born - rem_S
    E += +new_E - rem_E
    Ic += +new_Ic - rem_Ic
    Inc += +new_Inc - rem_Inc
    R += +new_R

    for i in 1:5
        A[i] += new_AFP[i] - new_AFP[i+1]
    end
    A[6] += new_AFP[6]

    model = SEIRModelOneStep(
        S=S, E=E, Ic=Ic, Inc=Inc, R=R, A=A,
        Z_A5_6=new_AFP[6]
    )
    return model
end

function run_sim(params::SEIRModelParams;
    rec_flag::Bool=false, pattern=""
)
    @unpack γ1, γ2, σ, P_AFP, P_H, P_AFP_sample, pc,
    P_AFP_test, P_ES_test, days, ES_n_freq, Pop_whole, ES_μ, ES_σ = params
    lognorm = LogNormal(ES_μ, ES_σ)

    rec, model = initialize_model(params, rec_flag; pattern=pattern)
    # set ES scheudle
    nt = [0 for _ in 1:days]
    st_ind = rand(1:ES_n_freq)
    nt[st_ind:ES_n_freq:days] .= 1


    t_extinct = NaN
    t_AFP = NaN
    t_ES = NaN
    R_final = NaN

    for t in 2:days
        # Stop condition
        if (model.E + model.Inc + model.Ic + sum(model.A[1:5])) == 0
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
            ωt = cdf(lognorm, model.Ic / (Pop_whole * pc) * 100_000)
            wt = rand_binom(nt[t], ωt * P_ES_test)
            t_ES = wt == 1 ? t : t_ES
        end

        # AFP surveillance
        if isnan(t_AFP) == true
            Ct = rand_binom(model.Z_A5_6, P_H * P_AFP_sample * P_AFP_test)
            t_AFP = Ct >= 1 ? t : t_AFP
        end

    end
    R_final_AFP = model.A[6]
    outcome = (t_ES=t_ES, t_AFP=t_AFP, t_extinct=t_extinct,
        R_final=R_final, R_final_AFP=R_final_AFP,
    )
    return (rec, outcome)
end

