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
include("geo_ana.jl")
include("model_meta_pop.jl")

path = "../dt_tmp/spatial_params_agg230.jld2"
spatial_p = load(path)["data"]
spatial_p |> keys |> println
df = spatial_p.df
π_mat = spatial_p.π_mat
nothing

path_zaf_0_4 = "../data_pop/zaf_merged_0_4.tif"
zaf_0_4 = read(Raster(path_zaf_0_4))
nothing

# # Probability of importation for each scenario

function visualise_importation_prob(df_mer; 
        title::String = "", colorbar = true,
        colorbar_title = "",
    )
    zaf_map = copy(zaf_0_4)
    zaf_map[:] .= 0
    n = size(df_mer)[1]
    for i in 1:n
        x = df_mer[i, :lon]
        y = df_mer[i, :lat]
        zaf_map[ At(x), At(y)] = log10.(df_mer[i, :imp_prob])
    end
    pl = plot()
    plot!(zaf_map,
        xlim=[15,35], ylim=[-35, -21],
        axis=nothing, border=:none,
        logscale=:log10,
        #colorbar_title="log10(πij)",
        title=title,
        right_margin=5Plots.mm,
        top_margin = 0Plots.mm,
        colorbar = colorbar,
        colorbar_title = colorbar_title,
        color = :matter,
        clims = (-6.5, -0.5),
        fmt=:png, dpi=200,
    )
    add_zaf_borders!(pl)
    return pl
end

sp_pars = read_spatial_params_file("ES_population_size")
imp_pop = normalize(sp_pars.pop, 1)
imp_ws_moz = CSV.read("../data/imp_ws_moz.csv", DataFrame, header=false)[:,1]
imp_ws_airport = CSV.read("../data/imp_ws_airport.csv", DataFrame, header=false)[:,1]
nothing

log10.(imp_pop) |> describe
log10.(imp_ws_airport) |> describe
log10.(imp_ws_moz) |> describe

# +
add_annotation!(pl, text_) =  annotate!(pl, 16, -22, text(text_, :black, :left, 18))

df[:, :imp_prob] = imp_pop
pl_pop = visualise_importation_prob(df; colorbar=false)
add_annotation!(pl_pop, "A")

df[:, :imp_prob] = imp_ws_airport
pl_airport = visualise_importation_prob(df; colorbar=false)
add_annotation!(pl_airport, "B")

df[:, :imp_prob] = imp_ws_moz
pl_moz = visualise_importation_prob(df; colorbar=false)
add_annotation!(pl_moz, "C")
layout = @layout [a b; c c]
plot(pl_pop, pl_airport, pl_moz, layout=layout, 
    fmt=:png, dpi=300, size=(1200, 800))
# -

df[:, :imp_prob] = imp_pop
pl_pop = visualise_importation_prob(df; colorbar=true, colorbar_title="log10(Probability of importation)")
plot!(pl_pop, size=(800,600))

# ## For figure 1 

df[:, :imp_prob] = imp_pop
pl_pop = visualise_importation_prob(df; colorbar=false)
plot!(pl_pop, right_margin=-25Plots.mm)
nothing

# +
# International airport information
ind = [10, 7, 61]
travel_volume = [4_342_611, 1_156_996, 188_243]
#weight = normalize((travel_volume).^(1/2) , 1) # Travel volume
weight = [0.7, 0.5, 0.3]
cm = palette(:tab10)

zaf_map = copy(zaf_0_4)
zaf_map[:] .= 0
pl_airport = plot()
plot!(pl_airport, zaf_map,
    xlim=[15,35], ylim=[-35, -21],
    axis=nothing, border=:none,
    logscale=:log10,
    #colorbar_title="log10(πij)",
    title="",
    right_margin=-10Plots.mm,
    top_margin = 0Plots.mm,
    colorbar = false,
    color = :matter,
    clims = (-6.5, -0.5),
    fmt=:png, dpi=200,
)
add_zaf_borders!(pl_airport)
# Plot the internatinoal airport positions
scatter!(pl_airport, df[ind, :lon], df[ind, :lat], 
    label=false, color=cm[1:3], 
    markershape=:circle,
    markersize=weight*50, alpha=1,
)

nothing
# -

include("geo_ana.jl")

pl_moz = plot()
plot!(pl_moz, zaf_map,
    axis=nothing, border=:none,
    title="",
    xlim=[15,45], ylim=[-35, -10],
)
add_zaf_ADM0_borders!(pl_moz)
add_moz_borders!(pl_moz)
scatter!(pl_moz, 
    [24, 33, 40, 33], 
    [-30, -15, -13, -22], 
    label=false, color=cm[4],
    markershape=:rect,
    markersize=[12, 8, 8, 8],
)
nothing

plot(pl_pop, pl_airport, pl_moz, 
    fmt=:png, size=(1800, 600), dpi=150,
    layout=@layout[a b c],
)

# ## Prepare ES layout for figure 1

function plot_ES_covered_sites(
    ES_pattern, r_ind
)
    spatial_p = read_spatial_params_file(ES_pattern)
    df = spatial_p.df
    
    # Prepare data
    zaf_map = copy(zaf_0_4)
    zaf_map[:] .= 0
    
    ind = 1:r_ind
    for r in eachrow(df[ind,:])
        zaf_map[ At(r.lon), At(r.lat)] = 100
    end
    pl = plot(zaf_map,
        xlim=[15,35], ylim=[-35, -21],
        color=:blue,
        colorbar=false,
        title="",
        axis=nothing, border=:none,
        dpi=100, fmt=:png, size=(1200,1200),
    )
    add_zaf_borders!(pl)
    pl
end

# This setting is for Figure 1. 
pl_pop = plot_ES_covered_sites("ES_population_size", 50)
plot!(pl_pop, right_margin=-20Plots.mm)
pl_moz = plot_ES_covered_sites("ES_mozambique_imp_risk", 50)
plot(pl_pop, pl_moz, fmt=:png, dpi=150, size=(1200,600))

# +
# South Africa's empirical data.
path = "../data/ES_surveillance_information_20240318.xlsx"
df_ES = @pipe XLSX.readtable(path, "new_calculation","A:I"; first_row=1) |> DataFrame

pop_size = "Population size served by the facility"
cond = (df_ES[:,  "Plant name"] .== "Rooiwal Eastern") .| (df_ES[:,  "Plant name"] .== "Daspoort")
df_ES[:, :color] .= :blue
df_ES[cond, :color] .= :grey
df_ES[cond, pop_size] .= 1000
df_ES[:, "GPS longitude"] = @pipe df_ES[:, "GPS longitude"] .|> parse(Float64, _)
df_ES[:, "GPS latitude"] = @pipe df_ES[:, "GPS latitude"] .|> parse(Float64, _)

nothing
# -

# # Comparison between observed and simulated wastewater sites. 

col_ES_pop = "Population size served by the facility"
df_ES[:, ["Plant name", col_ES_pop]]
cates = ["0 - 200,000", "200,001 - 500,000", "500,001 - 1,500,000", "NA"]
markershapes = [:x, :circle, :rect, :utriangle]
colors = [:purple, :cyan3, :tan3, :grey]
df_ES[!, :ES_category] = map(x -> 
    if x <= 1000;  cates[4]
    elseif x <= 200_000; cates[1]
    elseif x <= 500_000; cates[2]
    else cates[3]
    end, df_ES[:, col_ES_pop])
df_ES[:, :ES_color] = replace.(df_ES[:, :ES_category],  Dict(cates[i]=>colors[i] for i in 1:4)...)
df_ES[:, :markershape] = replace.(df_ES[:, :ES_category], Dict(cates[i]=>markershapes[i] for i in 1:4)...)
sort!(df_ES, col_ES_pop, rev=true)
nothing

function add_observed_ES_sites()
    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

    # Main figure
    pl = plot(zaf_map,
        xlim=[15,35], ylim=[-35, -21],
        colorbar=false,
        title="",
        axis=nothing, border=:none,
        dpi=100, fmt=:png, size=(1200,1200),
    )
    x_rec = 26.5; y_rec = -26.8; x_Δ = 2.5; y_Δ = 1.7
    plot!(rectangle(x_Δ,y_Δ, x_rec,y_rec), color=:white, 
        strokecolor=:black, alpha=0.7, linewidth=1.5, label="")
    scatter!(pl, df_ES[:, "GPS longitude"], df_ES[:, "GPS latitude"],  
        markershape=df_ES[:, :markershape] .|> Symbol,
        color=df_ES[:, :ES_color], markersize=4, markerstrokewidth=0.8,
        label=""
    )    
    # Custom label
    for (label, ms, c) in zip(cates, markershapes, colors)
        scatter!(pl, [1],[1], markershape=ms, label=label, 
            markersize=12, markerstrokewidth=1.2, color=c,
            foreground_color_legend=nothing,
            legend=(1.2, 0.8), 
            legendtitle="ES-covered Population",
            legendtitlefontsize=8
        )
    end
    add_zaf_borders!(pl)
    plot!(pl, right_margin=-10Plots.mm)
    # Inset
    plot!(pl, inset=bbox(0.65, 0.3, 0.15, 0.8), subplot=2)
    # Draw outline of the inset.
    plot!(pl[2], rectangle(x_Δ,y_Δ, x_rec,y_rec), color=:white, 
        strokecolor=:black, alpha=1.0, linewidth=1.5, label="")
    plot!(pl[2], zaf_map,
        xlim=[x_rec,x_rec+x_Δ], ylim=[y_rec, y_rec + y_Δ],
        colorbar=false,
        title="", axis=nothing, border=:none,
    )
    add_zaf_borders!(pl[2])
    # Plot ES sites.
    scatter!(pl[2], df_ES[:, "GPS longitude"], df_ES[:, "GPS latitude"],  
        markershape=df_ES[:, :markershape] .|> Symbol,
        color=df_ES[:, :ES_color], markersize=6, markerstrokewidth=0.8,
        label="", alpha=1.0,
    )    

    return pl
end



# +
pl = add_observed_ES_sites()
annotate!(pl[1], (0.05, 0.8), text("A", :left, :bottom, :black,24))
plot!(pl, right_margin=-10Plots.mm)

# Corerspond to 11.3% national ES coverage for each scenario.
ES_POP_num = 58
ES_LBC_num = 154
pl_pop = plot_ES_covered_sites("ES_population_size", ES_POP_num)
plot!(pl_pop, right_margin=-30Plots.mm)
annotate!(pl_pop, (0.05, 0.8), text("B", :left, :bottom, :black,24))
pl_moz = plot_ES_covered_sites("ES_mozambique_imp_risk", ES_LBC_num)
plot!(pl_moz, left_margin=-30Plots.mm)
annotate!(pl_moz, (0.05, 0.8), text("C", :left, :bottom, :black,24))

layout = @layout[a ; b c]
plot(pl, pl_pop, pl_moz, layout=layout, fmt=:png, dpi=150, size=(1200,600))
# -

# ## Prepare multiple figures with ES coverage sites for gif

function create_gif_resources(ES_pattern, dir_, title; pc=0.25)
    spatial_p = read_spatial_params_file(ES_pattern)
    df = spatial_p.df
    # Prepare data
    pop = df[:, :value]
    cum_per = normalize(pop, 1) .* 100 |> cumsum
    sens_index = obtain_ES_sensitivity_index(pop, 0.01)
    zaf_map = copy(zaf_0_4)
    zaf_map[:] .= 0
    
    # Create multiple files.
    l_ind =  1
    for (i, r_ind) in enumerate(sens_index)
        ind = l_ind:r_ind
        for r in eachrow(df[ind,:])
            zaf_map[ At(r.lon), At(r.lat)] = 100
        end
        l_ind = r_ind + 1
        #annotate!(pl, r.lon, r.lat, text(r.cum_per, :left, 4))
        pl = plot(zaf_map,
            xlim=[15,35], ylim=[-35, -21],
            color=:blue,
            colorbar=false,
            axis=nothing, border=:none,
            dpi=100, fmt=:png, size=(1200,1200),
            title="",
        )
        per = round(cum_per[r_ind] * pc,digits=1)
        annotate!(pl, (0.5, 1.15), text(title, :center, :center, 30))
        annotate!(pl, (0.2, 1.0), text("Natioanl ES population coverage: $(per)%\nNumber of ES-covered patches: $(r_ind)", :left, :bottom, 24))
        add_zaf_borders!(pl)
        file = lpad(i, 4, "0")
        savefig(pl, "$(dir_)/$(file).png")
    end
    # Make gif in python
end

ES_pattern = "ES_population_size"
dir_ = "../dt_tmp/gif_map"
create_gif_resources(ES_pattern, dir_, "ES-POP")

ES_pattern = "ES_mozambique_imp_risk"
dir_ = "../dt_tmp/gif_map_ES_LBC"
create_gif_resources(ES_pattern, dir_, "ES-LBC")

# ## Calculate number of ES-covered patches corresponding to the observed national ES coverage

function correspond_ES_covered_sites(spatial_p)
    df = spatial_p.df
    pop = df[:, :value]
    cum_per = normalize(pop, 1) |> cumsum |> x -> x .* 100
    sens_index = obtain_ES_sensitivity_index(pop, 0.01)
    
    # These index are for gif creation process in python.
    # Use the maximum index for the inclusion criteria.
    
    # sum(sens_index .<= n) is used to specify r_ind. 
    national_ES_cov = 11.3 # for old value, 8.59
    n = (cum_per .<= national_ES_cov/0.25) |> sum
    println("pc 25, maximum index: ", sum(sens_index .<= n), ", Number of sites: ", n)
    n = (cum_per .<= national_ES_cov/0.30) |> sum
    println("pc 30, maximum index: ", sum(sens_index .<= n), ", Number of sites: ", n)
    n = (cum_per .<= national_ES_cov/0.5) |> sum
    println("pc 50, maximum index: ", sum(sens_index .<= n), ", Number of sites: ", n)
end

spatial_p = read_spatial_params_file("ES_population_size")
correspond_ES_covered_sites(spatial_p)

spatial_p = read_spatial_params_file("ES_mozambique_imp_risk")
correspond_ES_covered_sites(spatial_p)

df = spatial_p.df
pop = df[:, :value]
normalize(pop, 1) |> sum 

# ## Visualise the propability of mobilisation from 3 top populous areas

include("utils.jl")

function visualise_probs(df_mer, ind::Int64, title::String; colorbar_title="")
    zaf_map = copy(zaf_0_4)
    zaf_map[:] .= 0
    n = size(df_mer)[1]
    for i in 1:n
        x = df_mer[i, :lon]
        y = df_mer[i, :lat]
        zaf_map[ At(x), At(y)] = log10(π_mat[i, ind])
    end
    pl = plot()
    plot!(zaf_map,
        xlim=[15,35], ylim=[-35, -21],
        axis=nothing, border=:none,
        logscale=:log10,
        colorbar_title=colorbar_title, #"log10(πij)",
        color=:matter,
        title=title,
    )
    add_zaf_borders!(pl)
    return pl
end

# +
ind = 1
pl1 = histogram(log10.(π_mat[ind, :]),
    legend=:none,
    xlabel = "log10(π1i)", ylabel="Number of patches",
)
annotate!(pl1, (0.1, 0.95), "A")
dist = df[ind, "shapeName"]
pl2 = visualise_probs(df, ind, "1st populous patch\n$(dist)";
    colorbar_title="log10(π1i)",
)
annotate!(pl2, (0.1, 0.95), "B")

ind = 2
pl3 = histogram(log10.(π_mat[ind, :]),
    legend=:none,
    xlabel = "log10(π2i)", ylabel="Number of patchs",
)
annotate!(pl3, (0.1, 0.95), "C")
dist = df[ind, "shapeName"]
pl4 = visualise_probs(df, ind, "2nd populous patch\n$(dist)",
    colorbar_title="log10(π2i)",
)
annotate!(pl4, (0.1, 0.95), "D")

ind = 3
pl5 = histogram(log10.(π_mat[ind, :]),
    legend=:none,
    xlabel = "log10(π3i)", ylabel="Number of patchs",
)
annotate!(pl5, (0.1, 0.95), "E")
dist = df[ind, "shapeName"]
pl6 = visualise_probs(df, ind, "3rd populous patch\n$(dist)",
    colorbar_title="log10(π3i)",
)
annotate!(pl6, (0.1, 0.95), "F")

l = @layout [
    a{0.3w} b
    c{0.3w} d
    e{0.3w} f
]
pl = plot(pl1, pl2, pl3, pl4, pl5, pl6,
    dpi=300, size=(800, 350 * 3), fmt=:png,
    layout=l, left_margin=5Plots.mm,
)
display(pl)
savefig(pl, "../res/fig_pi_map.png")
# -
# ## Unvaccinated base

path = "../dt_tmp/spatial_params_agg230_unvac.jld2"
spatial_p |> keys |> println
df = spatial_p.df
π_mat = spatial_p.π_mat
nothing

# +
ind = 1
pl1 = histogram(log10.(π_mat[:, ind]),
    legend=:none,
    xlabel = "log10(πi1)", ylabel="Number of locations",
)
annotate!(pl1, (0.1, 0.95), "(A)")
dist = df[ind, "shapeName"]
pl2 = visualise_probs(df, ind, "1st populous location\n$(dist)")
annotate!(pl2, (0.1, 0.95), "(B)")

ind = 2
pl3 = histogram(log10.(π_mat[:, ind]),
    legend=:none,
    xlabel = "log10(πi2)", ylabel="Number of locations",
)
annotate!(pl3, (0.1, 0.95), "(C)")
dist = df[ind, "shapeName"]
pl4 = visualise_probs(df, ind, "2nd populous location\n$(dist)")
annotate!(pl4, (0.1, 0.95), "(D)")

ind = 3
pl5 = histogram(log10.(π_mat[:, ind]),
    legend=:none,
    xlabel = "log10(πi3)", ylabel="Number of locations",
)
annotate!(pl5, (0.1, 0.95), "(E)")
dist = df[ind, "shapeName"]
pl6 = visualise_probs(df, ind, "3rd populous location\n$(dist)")
annotate!(pl6, (0.1, 0.95), "(F)")

l = @layout [
    a{0.3w} b
    c{0.3w} d
    e{0.3w} f
]
pl = plot(pl1, pl2, pl3, pl4, pl5, pl6,
    dpi=300, size=(800, 350 * 3),
    layout=l, left_margin=5Plots.mm,
)
display(pl)
savefig(pl, "../res/fig_pi_map_unvac.png")
# -






