using Distributions
using Interpolations

rand_binom(n,p) = rand(Binomial(n,p))

"""
...
# Arguments
- `days`: Day for each datapoint.
- `virus` : Amount of virus shedding.
...
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