using Accessors
using Distributions
using Parameters
using Plots
using ProtoStructs
using ProgressMeter
using Random
using StatsBase

include("util.jl")

@proto struct BaseParams
    R0::Float64 = 10
    γ1::Float64 = 1/4
    γ2::Float64 = 1/24
    γ3::Float64 = 1/5
    p_AFP::Float64 = 1/200
    β::Float64 = R0*γ2
    days::Int64 = 365
    N0::Int64 = 10_000
    I0_init::Int64 = 1
    λ0::Float64 = 76.0
end

@proto struct BaseSIRModel
    days::Int64
    S::Array{Int64, 1} = fill(-1000, days)
    E::Array{Int64, 1} = fill(-1000, days)
    Ia::Array{Int64, 1} = fill(-1000, days)
    I_AFP::Array{Int64, 1} = fill(-1000, days)
    R::Array{Int64, 1} = fill(-1000, days)
    H::Array{Int64, 1} = fill(-1000, days)
    I_new::Array{Int64, 1} = fill(-1000, days)
    H_new::Array{Int64, 1} = fill(-1000, days)
end

function set_values(model::BaseSIRModel,
                    ind::Int64,
                    S::Int64,
                    E::Int64,
                    Ia::Int64,
                    I_AFP::Int64,
                    R::Int64,
                    H::Int64,
                    I_new::Int64,
                    H_new::Int64,
                    )
    model.S[ind] = S
    model.E[ind] = E
    model.Ia[ind] = Ia
    model.I_AFP[ind] = I_AFP
    model.R[ind] = R
    model.H[ind] = H
    model.I_new[ind] = I_new
    model.H_new[ind] = H_new
    return nothing
end

function get_values(model::BaseSIRModel, ind::Int64)
    return (
        model.S[ind], model.E[ind], model.Ia[ind], model.I_AFP[ind],
        model.R[ind], model.H[ind], 
    )
end

function initialize_model(p::BaseParams)
    @unpack N0, I0_init, days = p
    model = BaseSIRModel(days=days)
    S0 = N0
    S0 -= I0_init
    Ia = I0_init
    set_values(model, 1, S0, 0, Ia, 0, 0, 0, 0, 0)
    return model
end
    
function run_sim(p::BaseParams; verbose=false)
    @unpack γ1, γ2, γ3, p_AFP, β, days, N0, I0_init = p
    model = initialize_model(p)

    for t in 2:days
        S, E, Ia, I_AFP, R, H = get_values(model, t-1)

        p_inf = β/N0*(Ia + I_AFP)
        new_inf = rand_binom(S, 1 .- exp(-p_inf))
        rec_E = rand_binom(E, 1 .- exp(-γ1))
        rec_E_Ia = rand_binom(rec_E, 1 - p_AFP)
        rec_E_I_AFP = rec_E - rec_E_Ia
        rec_Ia = rand_binom(Ia, 1 .- exp(-γ2))
        rec_I_AFP = rand_binom(I_AFP, 1 .- exp(-γ3))

        S = S - new_inf
        E = E + new_inf - rec_E
        Ia = Ia         + rec_E_Ia    - rec_Ia
        I_AFP = I_AFP   + rec_E_I_AFP - rec_I_AFP
        R = R + rec_Ia
        H = H + rec_I_AFP
        set_values(model, t, S, E, Ia, I_AFP, R, H, rec_E, rec_I_AFP)
    end
    return model
end

function environmental_sample_process(days, I_new::Vector; λ0=76)
    st_ind = rand(1:30)

    ts = 1:days
    λ = [hazard_function(res.I_new, gs, t; λ0=76)  for t in ts]
    ω = 1 .- exp.(-λ)
    ns = [0 for _ in 1:p.days]
    inds = st_ind:30:p.days
    ns[inds] .= 1
    w_det = rand_binom.(ns, ω)
    return w_det
end

function fetch_first_date(p::BaseParams, res::BaseSIRModel)
    d_afp1 = findfirst(x -> x > 0, res.H_new)
    wt = environmental_sample_process(p.days, res.I_new; λ0=76)
    d_wt1 = findfirst(x -> x > 0, wt)
    return d_afp1, d_wt1
end



