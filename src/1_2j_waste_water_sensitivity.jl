# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.7
#   kernelspec:
#     display_name: Julia 1.8.3
#     language: julia
#     name: julia-1.8
# ---

include("utils.jl")

# ## Lognormal for Fuqing Wu 2021

function cum_lognormal(x, data_x, data_y)
    if x[2] < 0; return Inf; end 
    lognorm = LogNormal(x[1], x[2])
    
    prob = cdf.(lognorm, data_x)
    prob = [ maximum([x, eps(0.0)]) for x in prob]
    lklh = [
        if data_y[i] == 1; prob[i]; else 1 - prob[i]; end
        for i in 1:length(data_y)
    ]
    ll = log.(lklh) |> sum
    return -ll
end

cov_RNA = "sars.cov.2_cocentration_copy_per_ml*"
col_det = :detection_status_SARSCoV2
col_new_cases = "new_cases***"

df = CSV.read("../data/fuqing_wu_2021.csv", DataFrame)
println("Original data: ", size(df))
n_sample_loc = df[:, :Sampling_location_ID] |> unique |> length
println("# of sampling locations (we substracted one for missing): ", n_sample_loc - 1)
df = filter(:Sampling_location_ID => x -> ismissing(x) == false, df)
df_fil = filter(col_det => x -> ismissing(x) == false, df)
df_fil = filter(:Estimated_population_served => x -> ismissing(x) == false, df)
println(size(df_fil))
df_fil[:, :flag_det] .=  [ 
    if x == "quantifiable"; 1; else 0; end 
    for x in df_fil[:, col_det] 
]
println("# of rows with no missing for flag and estimated pop: ", nrow(df_fil))
nothing

# Obtain location with positives and negatives
cnts = df_fil[:, :Sampling_location_ID] |> unique
cnts_pos_neg = []
cnts_virus = []
ind = 0
for c in cnts
    dfM = filter(:Sampling_location_ID => x -> x == c, df_fil)
    df_pos = filter(:flag_det => x -> x == 1, dfM)
    df_neg = filter(:flag_det => x -> x == 0, dfM)
    n_pos = nrow(df_pos)
    n_neg = nrow(df_neg)

    if (n_pos > 0) & (n_neg > 0) 
        push!(cnts_pos_neg, c)
    else
        continue
    end
    # further filter the virus concentration.
    pos_min = minimum(df_pos[:, col_new_cases])
    neg_max = maximum(df_neg[:, col_new_cases])
    if neg_max < pos_min
        push!(cnts_virus, c)
    end
end
filter!(:Sampling_location_ID => x -> x in cnts_virus, df_fil)
println("# of sites with pos neg: ", length(cnts_pos_neg))
println("# of rows with pos neg sites: ", length(cnts_pos_neg))
println("# of rows with consistent detection: ", length(cnts_virus))
println("# of rows : ", nrow(df_fil))
nothing

#data_x = df_fil[:, "new_cases***"] ./ df_fil[:, :Estimated_population_served] .* 100_000
data_x = df_fil[:, "new_cases***"] ./ df_fil[:, :county_population] .* 100_000
data_y = df_fil[:, :flag_det]
nothing

scatter(
    log10.(df_fil[:, :Estimated_population_served]), 
    log10.(df_fil[:, :county_population])
)

histogram(log10.(data_x), xlabel="new reported cases per 100,000", ylabel="Count")

inits = [0.5, 1.1]
fit_lognorm(data_x, data_y, inits) = optimize(
    x-> cum_lognormal(x, data_x, data_y), inits
) |> Optim.minimizer

res = fit_lognorm(data_x, data_y, [1.0, 5.0])

lognorm = LogNormal(res[1], res[2])
q50 = quantile(lognorm, 0.5)
q80 = quantile(lognorm, 0.8)
println("Quantil 50%: $(q50)")
println("Quantil 80%: $(q80)")

pl = plot()
scatter!(pl, data_x, data_y)
y = cdf.(lognorm, data_x)
scatter!(pl, data_x, y)

pl = plot(fmt=:png, dpi=150)
labels = ["10x lower", "3x lower", "Baseline", "3x higher", "10x higher"]
for (i, k) in enumerate([10, 3, 1, 1/3, 1/10])
    res = fit_lognorm(data_x .* k, data_y, [1.0, 5.0])
    println(f"Estimated params: {res[1]:.3f}, {res[2]:.3f}")

    lognorm = LogNormal(res[1], res[2])
    x = 0:0.1:50
    y = cdf(lognorm, x)
    plot!(pl, x, y, label=labels[i], lw=2,
        xlabel="Infectious individuals per 100,000 population",
        ylabel="ES sensitivity",
    )
end
display(pl)







