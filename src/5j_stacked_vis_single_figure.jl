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
include("model_meta_pop.jl")
include("visualise_fig.jl")

# # Data preparation

path_res1 = "../dt_tmp_hpc/sens_ES_catchment_20240302_050027059.jld2"
path_res2 = "../dt_tmp_hpc/sens_ES_catchment_20240302_091756993.jld2"
path_res3 = "../dt_tmp_hpc/sens_ES_catchment_20240302_190345473.jld2"
path_res1_moz = "../dt_tmp_hpc/sens_ES_catchment_20240303_020804257.jld2"
path_res2_moz = "../dt_tmp_hpc/sens_ES_catchment_20240303_064157726.jld2"
path_res3_moz = "../dt_tmp_hpc/sens_ES_catchment_20240303_160629569.jld2"

fontsize = 8 
xlabel = "Number of ES-covered patchs"
ylabel = "Proportion of detection pattern (%)"
vis_kwds = (
    left_margin=0Plots.mm, right_margin=0Plots.mm,
    xlabelfontsize=fontsize, ylabelfontsize=fontsize, tickfontsize=fontsize -2, 
    foreground_color_legend=nothing, background_color_legend=nothing,
)
kwds = (x_var="site", xlim=[0,160], legend=false)
nothing

# ## For Figure 1

fontsize1 = 20
vis_kwds1 = (
    xlabelfontsize=fontsize1, ylabelfontsize=fontsize1, 
    tickfontsize=fontsize1-4,
)
pl1 = single_stacked_heatmap(path_res1; 
    vis_kwds=vis_kwds1,
    xlabel=xlabel, ylabel=ylabel, kwds...)
plot!(pl1, dpi=200, fmt=:png, 
    right_margin=5Plots.mm, left_margin=5Plots.mm,
    top_margin=5Plots.mm, size=(700, 600),
)

# ## Main Figure

# +
pl1 = single_stacked_heatmap(path_res1; vis_kwds=vis_kwds, 
    ylabel=ylabel, kwds...)
plot!(pl1, left_margin=5Plots.mm)
#add_annotation!(pl1, "(A) Population size")
add_annotation1!(pl1, "A", "IMP-POP/ES-POP")

pl2 = single_stacked_heatmap(path_res2; vis_kwds=vis_kwds, kwds...)
#add_annotation!(pl2, "(B) Airport")
add_annotation1!(pl2, "B", "IMP-AIR/ES-POP")

pl3 = single_stacked_heatmap(path_res3; vis_kwds=vis_kwds, 
    x_var="site", xlim=[0,160], legend=true)
add_annotation1!(pl3, "C", "IMP-LBC/ES-POP")
plot!(pl3, right_margin=25Plots.mm)

pl4 = single_stacked_heatmap(path_res1_moz; 
    vis_kwds=vis_kwds, ylabel=ylabel, xlabel=xlabel, kwds...)
plot!(pl4, bottom_margin=3Plots.mm, )
#add_annotation!(pl4, "(D) Population size")
add_annotation1!(pl4, "D", "IMP-POP/ES-LBC")

pl5 = single_stacked_heatmap(path_res2_moz; 
    vis_kwds=vis_kwds, xlabel=xlabel, kwds...)
add_annotation1!(pl5, "E", "IMP-AIR/ES-LBC")

pl6 = single_stacked_heatmap(path_res3_moz; 
    vis_kwds=vis_kwds, xlabel=xlabel, kwds...)
add_annotation1!(pl6, "F", "IMP-LBC/ES-LBC")

pls = [pl1, pl2, pl3, pl4, pl5, pl6]
pl = plot(pls..., layout = @layout[a b c; d e f],  
    fmt=:png, dpi=300, size=(1000,500),)
display(pl)
savefig(pl, "../res/Fig2_6scenarios.png")
# -
# ## Main figure against the national ES population coverage (%)

xlabel = "National ES population coverage (%)"
kwds_cov = (x_var="coverage", xlim=[0,100], legend=false)
nothing

include("visualise_fig.jl")

# +
pl1 = single_stacked_heatmap(path_res1; vis_kwds=vis_kwds, 
    ylabel=ylabel, kwds_cov...)
plot!(pl1, left_margin=5Plots.mm)

pl2 = single_stacked_heatmap(path_res2; vis_kwds=vis_kwds, kwds_cov...)

pl3 = single_stacked_heatmap(path_res3; vis_kwds=vis_kwds, 
    x_var="coverage", xlim=[0,100], legend=true)
plot!(pl3, right_margin=25Plots.mm)

pl4 = single_stacked_heatmap(path_res1_moz; 
    vis_kwds=vis_kwds, ylabel=ylabel, xlabel=xlabel, kwds_cov...)
plot!(pl4, bottom_margin=3Plots.mm, )

pl5 = single_stacked_heatmap(path_res2_moz; 
    vis_kwds=vis_kwds, xlabel=xlabel, kwds_cov...)

pl6 = single_stacked_heatmap(path_res3_moz; 
    vis_kwds=vis_kwds, xlabel=xlabel, kwds_cov...)

add_annotation1!(pl1, "A", "IMP-POP/ES-POP")
add_annotation1!(pl2, "B", "IMP-AIR/ES-POP")
annotate!(pl3, (0.05, 0.85), text("C", :white, :left, :bottom, 13))
annotate!(pl3, (0.15, 0.85), text("IMP-ES/ES-POP", :white, :left, :bottom, 11))
add_annotation1!(pl4, "D", "IMP-POP/ES-LBC")
add_annotation1!(pl5, "E", "IMP-AIR/ES-LBC")
add_annotation1!(pl6, "F", "IMP-LBC/ES-LBC")

pls = [pl1, pl2, pl3, pl4, pl5, pl6]
plot(pls..., layout = @layout[a b c; d e f],  
    fmt=:png, dpi=300, size=(1000,500),)
# -
# ## Main figure against the number of ES-covered patches including no detection pattern

include("visualise_fig.jl")

kwds1 = (x_var="site", xlim=[0,160], ylim=[0,35], )
n_sim = 10_000
xlabel= "Number of ES-covered patches"
ylabel= "Proportion of detection pattern (%)"

# +
pl1 = single_stacked_heatmap_with_no_detection(path_res1, n_sim; 
    vis_kwds=vis_kwds,  ylabel=ylabel, legend=false,
    kwds1...,
)
plot!(pl1, left_margin=5Plots.mm)

pl2 = single_stacked_heatmap_with_no_detection(path_res2, n_sim; 
    vis_kwds=vis_kwds, legend=false,
    kwds1...,
)
pl3 = single_stacked_heatmap_with_no_detection(path_res3, n_sim; 
    vis_kwds=vis_kwds, legend=true,
    kwds1...,
)
plot!(pl3, right_margin=25Plots.mm)

pl4 = single_stacked_heatmap_with_no_detection(path_res1_moz, n_sim; 
    vis_kwds=vis_kwds, xlabel=xlabel, ylabel=ylabel, legend=false,
    kwds1...,
)

pl5 = single_stacked_heatmap_with_no_detection(path_res2_moz, n_sim; 
    vis_kwds=vis_kwds, xlabel=xlabel, legend=false,
    kwds1...,
)
pl6 = single_stacked_heatmap_with_no_detection(path_res3_moz, n_sim; 
    vis_kwds=vis_kwds, xlabel=xlabel, legend=false,
    kwds1...,
)

add_annotation1_no_detect!(pl1, "A", "IMP-POP/ES-POP");  add_annotation1_no_detect!(pl2, "B", "IMP-AIR/ES-POP")
add_annotation1_no_detect!(pl3, "C", "IMP-LBC/ES-POP");  add_annotation1_no_detect!(pl4, "D", "IMP-POP/ES-LBC")
add_annotation1_no_detect!(pl5, "E", "IMP-AIR/ES-LBC");  add_annotation1!(pl6, "F", "IMP-LBC/ES-LBC")

pls = [pl1, pl2, pl3, pl4, pl5, pl6]
plot(pls..., layout = @layout[a b c; d e f],  
    fmt=:png, dpi=300, size=(1000,500), bottom_margin=5Plots.mm)
# -

# # Sensitivity results

include("visualise_fig.jl")

# ### R0

tuples = [(R0 = 10, α = 0.05, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_res/sens_ES_catchment_20240302_025116224.jld2"), (R0 = 12, α = 0.05, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_res/sens_ES_catchment_20240302_082109968.jld2"), (R0 = 16, α = 0.05, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_res/sens_ES_catchment_20240302_162056133.jld2"), (R0 = 18, α = 0.05, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_res/sens_ES_catchment_20240303_012857453.jld2")]
for t in tuples
    println("path_R0_$(t.R0) = \"$(t.path)\"")
end
path_R0_10 = "../dt_tmp_hpc/sens_ES_catchment_20240302_025116224.jld2"
path_R0_12 = "../dt_tmp_hpc/sens_ES_catchment_20240302_082109968.jld2"
path_R0_16 = "../dt_tmp_hpc/sens_ES_catchment_20240302_162056133.jld2"
path_R0_18 = "../dt_tmp_hpc/sens_ES_catchment_20240303_012857453.jld2"
path_Rs = [path_R0_10, path_R0_12, path_res1, path_R0_16, path_R0_18]

# ## α

tuples = [(R0 = 14, α = 0.005, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_064121817.jld2"), (R0 = 14, α = 0.01, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_121639592.jld2"), (R0 = 14, α = 0.1, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_184251996.jld2"), (R0 = 14, α = 0.5, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240304_011519518.jld2")]
for t in tuples
    println("path_α_$(t.α) = \"$(t.path)\"")
end
path_α_0005 = "../dt_tmp_hpc/sens_ES_catchment_20240303_064121817.jld2"
path_α_001 = "../dt_tmp_hpc/sens_ES_catchment_20240303_121639592.jld2"
path_α_01 = "../dt_tmp_hpc/sens_ES_catchment_20240303_184251996.jld2"
path_α_05 = "../dt_tmp_hpc/sens_ES_catchment_20240304_011519518.jld2"
path_αs = [path_α_0005, path_α_001, path_res1, path_α_01, path_α_05]

# ## freq

tuples = [(R0 = 14, α = 0.05, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_144128963.jld2", n_freq = 1), (R0 = 14, α = 0.05, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_185355250.jld2", n_freq = 7), (R0 = 14, α = 0.05, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_211900139.jld2", n_freq = 14), (R0 = 14, α = 0.05, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_220740874.jld2", n_freq = 60)]
for t in tuples
    println("path_freq_$(t.n_freq) = \"$(t.path)\"")
end
path_freq_1 = "../dt_tmp_hpc/sens_ES_catchment_20240303_144128963.jld2"
path_freq_7 = "../dt_tmp_hpc/sens_ES_catchment_20240303_185355250.jld2"
path_freq_14 = "../dt_tmp_hpc/sens_ES_catchment_20240303_211900139.jld2"
path_freq_60 = "../dt_tmp_hpc/sens_ES_catchment_20240303_220740874.jld2"
path_freqs = [path_freq_1, path_freq_7, path_freq_14, path_res1,  
    path_freq_60]

tuples = [(R0 = 14, α = 0.05, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_052731428.jld2", ES_μ = 3.121, ES_σ = 1.45), (R0 = 14, α = 0.05, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_070141277.jld2", ES_μ = 1.917, ES_σ = 1.45), (R0 = 14, α = 0.05, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_084258207.jld2", ES_μ = -0.281, ES_σ = 1.45), (R0 = 14, α = 0.05, pc = 0.25, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_100427100.jld2", ES_μ = -1.485, ES_σ = 1.45)]
for t in tuples
    println("path_det_$(t.ES_μ) = \"$(t.path)\"")
end
path_det_10low = "../dt_tmp_hpc/sens_ES_catchment_20240302_052731428.jld2"
path_det_3low = "../dt_tmp_hpc/sens_ES_catchment_20240302_070141277.jld2"
path_det_3high = "../dt_tmp_hpc/sens_ES_catchment_20240302_084258207.jld2"
path_det_10high = "../dt_tmp_hpc/sens_ES_catchment_20240302_100427100.jld2"
path_dets = [path_det_10low, path_det_3low, path_res1, 
    path_det_3high, path_det_10high]

include("visualise_fig.jl")

fontsize = 8
line5_set = (
    lw = [1, 1, 2.3, 2, 2],
    ls = [:dot, :solid, :solid, :solid, :dot],
    color = [:black, :black, :black, :grey, :grey],
    xlabelfontsize=fontsize, ylabelfontsize=fontsize, tickfontsize=fontsize - 2,
)
line5_set_freq = (
    lw = [1, 1, 2, 2.3, 2],
    ls = [:dot, :solid, :solid, :solid, :dot],
    color = [:black, :black, :grey, :black, :grey],
    xlabelfontsize=fontsize, ylabelfontsize=fontsize, tickfontsize=fontsize - 2,
)

# +
### Define summarised plots.
pl_R0 = plot_sens_adaptor(
    path_Rs, single_vis_R0!; 
    xlabel="",
    legendtitle="R0",
    line5_set...
)
plot!(pl_R0, left_margin=5Plots.mm)

pl_α = plot_sens_adaptor(
    path_αs, single_vis_α!; 
    legendtitle="α",
    line5_set...
)

pl_freq = plot_sens_adaptor(
    path_freqs, single_vis_sampling_freq!; 
    xlabel="",
    legendtitle="Sampling freq.",
    line5_set_freq...
)

pl_det = plot_sens_adaptor(
    path_dets, single_vis_ES_det!; 
    legend = (0.7, 0.4),
    legendtitle="ES sensitivity",
    labels = [
        "10x lower", "3x lower", "Baseline",
        "3x higher", "10x higher"
    ],
    line5_set...
)
nothing

# +
#####  Individual stacked plot.
xlabel = "Number of ES-covered patches" 
ylabel = "Proportion of detection pattern (%)"
### R0 
pl_R0_10= single_stacked_heatmap(path_R0_10; 
    ylabel=ylabel,
    vis_kwds=vis_kwds, kwds...)
add_annotation2!(pl_R0_10, "R0 = 10")
plot!(pl_R0_10, left_margin=3Plots.mm)

pl_R0_18= single_stacked_heatmap(path_R0_18; vis_kwds=vis_kwds, 
    x_var="site", xlim=[0, 160], legend=true)
add_annotation2!(pl_R0_18, "R0 = 18")
plot!(pl_R0_18, right_margin=25Plots.mm)

### α
pl_α_0005 = single_stacked_heatmap(path_α_0005; 
    ylabel=ylabel,
    vis_kwds=vis_kwds, xlabel=xlabel, kwds...)
add_annotation2!(pl_α_0005, "α = 0.005")

pl_α_05 = single_stacked_heatmap(path_α_05; 
    vis_kwds=vis_kwds, xlabel=xlabel, kwds...)
add_annotation2!(pl_α_05, "α = 0.500")
plot!(pl_α, left_margin=5Plots.mm)

### Sampling freq
pl_freq_1 = single_stacked_heatmap(path_freq_1; 
    ylabel=ylabel,
    vis_kwds=vis_kwds, kwds...)
add_annotation2!(pl_freq_1, "Sampling freq. = 1 day")

pl_freq_60= single_stacked_heatmap(path_α_05; vis_kwds=vis_kwds,
    x_var="site", xlim=[0, 160], legend=true)
add_annotation2!(pl_freq_60, "Sampling freq. = 60 day")

### ES detection sensitivity
pl_det_10low = single_stacked_heatmap(path_det_10low; 
    ylabel=ylabel, 
    vis_kwds=vis_kwds, xlabel=xlabel, kwds...)
#add_annotation2!(pl_det_10low, "10x lower ES sensitivity")
annotate!(pl_det_10low, (0.05, 0.85), text("10x lower ES sensitivity", :white, :left, :bottom, 11))

pl_det_10high = single_stacked_heatmap(path_det_10high; 
    vis_kwds=vis_kwds, xlabel=xlabel, kwds...)
add_annotation2!(pl_det_10high, "10x higehr ES sensitivity")
nothing
# -

plot!(pl_R0, top_margin=10Plots.mm)
annotate!(pl_R0, (0.05, 1.00), text("A", :black, :left, :bottom, 24))
annotate!(pl_R0_10, (0.05, 1.00), text("B", :black, :left, :bottom, 24))
annotate!(pl_R0_18, (0.05, 1.00), text("C", :black, :left, :bottom, 24))
nothing

### Synthesize above all.
pls = [
    pl_R0, pl_R0_10, pl_R0_18,
    pl_α, pl_α_0005, pl_α_05,
    pl_freq, pl_freq_1, pl_freq_60,
    pl_det, pl_det_10low, pl_det_10high,
]
layout = @layout [a b c; d e f; g h i; j k l]
pl = plot(pls..., layout=layout, size=(1100,250*4),
    fmt=:png, dpi=300,
)
display(pl)
savefig(pl, "../res/Fig3_sens.png")

# # Results for pc. 

include("visualise_fig.jl")

tuples = [(R0 = 14.0, α = 0.05, pc = 0.01, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_050431173.jld2"), (R0 = 14.0, α = 0.05, pc = 0.05, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_114232977.jld2"), (R0 = 14.0, α = 0.05, pc = 0.5, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_182700909.jld2"), (R0 = 14.0, α = 0.05, pc = 1.0, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_005806947.jld2")]
path_pcs = prepare_pc_paths_from_tuples(tuples, path_res1)

pls = three_panels_for_pc(path_pcs)
annotate!(pls[1], (0.05, 0.87), text("A", :left, :bottom, :black, 20))
annotate!(pls[2], (0.05, 0.87), text("B", :left, :bottom, :black, 20))
annotate!(pls[3], (0.05, 0.87), text("C", :left, :bottom, :black, 20))
add_three_scatters!(pls, path_pcs)
plot(pls...,
    fmt=:png, dpi=300, size=(960, 300), layout=@layout[a b c]
)

tuples = [(R0 = 14.0, α = 0.05, pc = 0.01, n_sim = 10000, pattern = "airport", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_025324354.jld2"), (R0 = 14.0, α = 0.05, pc = 0.05, n_sim = 10000, pattern = "airport", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_070622591.jld2"), (R0 = 14.0, α = 0.05, pc = 0.5, n_sim = 10000, pattern = "airport", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_112053541.jld2"), (R0 = 14.0, α = 0.05, pc = 1.0, n_sim = 10000, pattern = "airport", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_152900513.jld2")]
path_pcs_airport = prepare_pc_paths_from_tuples(tuples, path_res2)

tuples = [(R0 = 14.0, α = 0.05, pc = 0.01, n_sim = 10000, pattern = "mozambique", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_084141929.jld2"), (R0 = 14.0, α = 0.05, pc = 0.05, n_sim = 10000, pattern = "mozambique", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_184104084.jld2"), (R0 = 14.0, α = 0.05, pc = 0.5, n_sim = 10000, pattern = "mozambique", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_035611430.jld2"), (R0 = 14.0, α = 0.05, pc = 1.0, n_sim = 10000, pattern = "mozambique", ES_pattern = "ES_population_size", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_124509645.jld2")]
path_pcs_moz = prepare_pc_paths_from_tuples(tuples, path_res3)

pls_airport = three_panels_for_pc(path_pcs_airport)
pls_moz = three_panels_for_pc(path_pcs_moz)
add_three_scatters!(pls_airport, path_pcs_airport)
add_three_scatters!(pls_moz, path_pcs_moz)

for (i, letter) in enumerate(["A","B","C"])
    add_annotation3_scenario!(pls_airport[i], letter, "IMP-AIR/ES-POP")
end
for (i, (letter, pos)) in enumerate(zip(["D","E","F"], [:upper, :upper, :middle]))
    add_annotation3_scenario!(pls_moz[i], letter, "IMP-LBC/ES-POP"; position=pos)
end
plot(pls_airport..., pls_moz...,
    layout= @layout[a b c; d e f],
    fmt=:png, dpi=300, size=(960, 300*2),
    #right_margin=5Plots.mm,
)

tuples = [(R0 = 14.0, α = 0.05, pc = 0.01, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_mozambique_imp_risk", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_055244568.jld2"), (R0 = 14.0, α = 0.05, pc = 0.05, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_mozambique_imp_risk", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_130717570.jld2"), (R0 = 14.0, α = 0.05, pc = 0.5, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_mozambique_imp_risk", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_200904149.jld2"), (R0 = 14.0, α = 0.05, pc = 1.0, n_sim = 10000, pattern = "population_size", ES_pattern = "ES_mozambique_imp_risk", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_030935200.jld2")]
path_pcs_pop = prepare_pc_paths_from_tuples(tuples, path_res1_moz)

tuples = [(R0 = 14.0, α = 0.05, pc = 0.01, n_sim = 10000, pattern = "airport", ES_pattern = "ES_mozambique_imp_risk", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_032641945.jld2"), (R0 = 14.0, α = 0.05, pc = 0.05, n_sim = 10000, pattern = "airport", ES_pattern = "ES_mozambique_imp_risk", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_081557465.jld2"), (R0 = 14.0, α = 0.05, pc = 0.5, n_sim = 10000, pattern = "airport", ES_pattern = "ES_mozambique_imp_risk", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_125615040.jld2"), (R0 = 14.0, α = 0.05, pc = 1.0, n_sim = 10000, pattern = "airport", ES_pattern = "ES_mozambique_imp_risk", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_183123788.jld2")]
path_pcs_airport = prepare_pc_paths_from_tuples(tuples, path_res2_moz)

tuples = [(R0 = 14.0, α = 0.05, pc = 0.01, n_sim = 10000, pattern = "mozambique", ES_pattern = "ES_mozambique_imp_risk", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_091032122.jld2"), (R0 = 14.0, α = 0.05, pc = 0.05, n_sim = 10000, pattern = "mozambique", ES_pattern = "ES_mozambique_imp_risk", path = "../dt_tmp_hpc/sens_ES_catchment_20240302_193055492.jld2"), (R0 = 14.0, α = 0.05, pc = 0.5, n_sim = 10000, pattern = "mozambique", ES_pattern = "ES_mozambique_imp_risk", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_051152393.jld2"), (R0 = 14.0, α = 0.05, pc = 1.0, n_sim = 10000, pattern = "mozambique", ES_pattern = "ES_mozambique_imp_risk", path = "../dt_tmp_hpc/sens_ES_catchment_20240303_140256713.jld2")]
path_pcs_moz = prepare_pc_paths_from_tuples(tuples, path_res3_moz)

pls_pop = three_panels_for_pc(path_pcs_pop)
pls_airport = three_panels_for_pc(path_pcs_airport)
pls_moz = three_panels_for_pc(path_pcs_moz)

# +
x_rel = 0.71
y_rels = [0.20, 0.15, 0.09]
pos = (x_rel=x_rel, y_rels=y_rels)
add_three_scatters!(pls_pop, ["A", "B","C"], path_pcs_pop; pos...)
add_three_scatters!(pls_airport, ["D", "E","F"], path_pcs_airport; pos...)
add_three_scatters!(pls_moz, ["G", "H","I"], path_pcs_moz; pos...)
for (i, (letter, pos)) in enumerate(zip(["A","B","C"], [:upper, :upper, :lower]))
    add_annotation3_scenario!(pls_pop[i], letter, "IMP-POP/ES-LBC"; position=pos)
end
for (i, (letter, pos)) in enumerate(zip(["D","E","F"], [:upper, :upper, :upper]))
    add_annotation3_scenario!(pls_airport[i], letter, "IMP-AIR/ES-LBC"; position=pos)
end
for (i, (letter, pos)) in enumerate(zip(["G","H","I"], [:lower, :lower, :lower]))
    add_annotation3_scenario!(pls_moz[i], letter, "IMP-LBC/ES-LBC"; position=pos)
end

plot(pls_pop..., pls_airport..., pls_moz...,
    layout= @layout[a b c; d e f; g h i ],
    fmt=:png, dpi=300, size=(960, 300*3),
)
# -











