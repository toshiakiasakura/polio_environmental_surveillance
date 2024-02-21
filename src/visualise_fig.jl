##### Functions for a single population model #####

function plot_cum!(pl, df; label="", color=:blue)
    cum_ES = cumulative_counts(df[:, "t_ES"], days; prop=true)
    cum_AFP = cumulative_counts(df[:, "t_AFP"], days; prop=true)
    ES_label = label == :none ? :none : "ES $(label)"
    AFP_label = label == :none ? :none : "AFP $(label)"
    plot!(pl, 1:days, cum_ES, label=ES_label, color=color)
    plot!(pl, 1:days, cum_AFP, label=AFP_label, color=color, ls=:dash)
end

function draw_cumulative_incidence(
        par_lis, colors;
        legendtitle="", labels="",
        xlabel="Day", ylabel="Cumulative\ndetection probability",
        legend=:bottomright
    )
    Random.seed!(123)

    pl = plot(
        xlabel=xlabel, ylabel=ylabel,
        legend=legend,
    )
    for (par, c) in zip(par_lis, colors)
        df = multiple_run(par)
        plot_cum!(pl, df; color=c, label=:none)
    end

    # Label setting
    scatter!([1], [0], ms=0, mc=:white, label=legendtitle)
    for (label, c) in zip(labels, colors)
        plot!([1], [0], color=c, label=" "^3 * label)
    end

    #scatter!([1], [0], ms=0, mc=:white, label=" ")
    #scatter!([1], [0], ms=0, mc=:white, label="Detect. type")
    #plot!([1], [0], color=:black, label="   ES")
    #plot!([1], [0], color=:black, ls=:dot, label="   AFP surv.")

    return pl
end

"""Calculate any probability given each simulation
and the probability of each detection pattern.

Args:
- `dfM`: Simulation result including no detection.
"""
function calculate_any_prob_and_each_percentage(
        dfM::DataFrame
    )::Tuple{Float64, DataFrame}
    dfM_fil, bin_labels = create_lead_time_category(dfM)
    p_any = nrow(dfM_fil)/nrow(dfM)
    freq = combine(groupby(dfM_fil, :lead_time_category), nrow => :freq)
    freq[:, :prop] = freq[:, :freq]/nrow(dfM_fil)
    return(p_any, freq)
end

function obtain_p_any_and_freq_list(par_lis)
    p_any_lis = []
    freq_lis = []
    for par in par_lis
        dfM = multiple_run(par)
        p_any, freq = calculate_any_prob_and_each_percentage(dfM)
        push!(p_any_lis, p_any)
        push!(freq_lis, freq)
    end
    return(p_any_lis, freq_lis)
end

function plot_groupedbar(
        cate, freq_lis::Vector;
        xlabel="",
        ylabel="Probability of \ndetection pattern (%)"
    )
    bin_labels = ["AFP only", "<-60 LT", "-60 ~ -1 LT", "0 ~ 59 LT", "≥60 LT", "ES only"]
    colors = discretise_balance_color(bin_labels)
    x = [x[:, :prop] for x in freq_lis]
    x_grouped = mapreduce(permutedims, vcat, x) # reduce vector of vector to matrix
    pl = groupedbar(cate, x_grouped,
        xlabel=xlabel, ylabel=ylabel,
        bar_position=:stack,
        labels=reshape(bin_labels, 1, 6),
        color=colors[:, end:-1:begin],
        legend=(1.1,0.5),
        right_margin=30Plots.mm,
    )
    return pl
end

###### Functions for meta-population model #####

"""
    Following functions are to visualise the heatmap including "only" patterns.
    - create_lead_time_category
    - proportion_each_cate_by_group
    - discretise_balance_color
    - visualise_heatmap_include_only
    - add_reverse_order_legend!
"""
function create_lead_time_category(df_res::DataFrame)::Tuple{DataFrame, Vector}
    # Replace NaN with Inf for t_AFP and t_ES
    df_res.t_AFP = ifelse.(isnan.(df_res.t_AFP), Inf, df_res.t_AFP)
    df_res.t_ES = ifelse.(isnan.(df_res.t_ES), Inf, df_res.t_ES)

    # Calculate the lead_time
    df_res.lead_time = df_res.t_AFP .- df_res.t_ES
    df_fil = filter(x -> !isnan.(x.lead_time), df_res)

    # Use "extend" option for cut, and assin Inf as "ES only".
    #bin_edges = [-Inf, -10_000, -150, -100, -50, 0, 50, 100, 150, 10_000, Inf]
    #bin_labels = ["AFP surv. only",
    #    "<-150 LT", "-150 ~ -101 LT", "-100 ~ -49 LT", "-50 ~ -1 LT",
    #    "0 ~ 49 LT", "50 ~ 99 LT",
    #    "100 ~ 149 LT", ">=150 LT",
    #    "ES only"]
    #bin_edges = [Inf, 10_000, 60, 0, -60, -10_000, -Inf]
    #bin_labels = ["ES only", "≥60 LT", "0 ~ 59 LT", "-60 ~ -1 LT", "<-60 LT", "AFP only"]
    bin_edges = [-Inf, -10_000, -60, 0, 60, 10_000, Inf]
    bin_labels = ["AFP only", "<-60 LT", "-60 ~ -1 LT", "0 ~ 59 LT", "≥60 LT", "ES only"]
    df_fil.lead_time_category= cut(df_fil.lead_time, bin_edges, extend=true, labels=bin_labels)
    return (df_fil, bin_labels)
end

function proportion_each_cate_by_group(df_fil::DataFrame, col::Symbol)::Tuple{Matrix, DataFrame}
    # Group by "ind_site" and calculate the proportion
    grp_prop = DataFrames.combine(groupby(df_fil, [col, :lead_time_category]), nrow => :count)
    df_fil, bin_labels = create_lead_time_category(df_fil)
    uni = grp_prop[:, col] |> unique
    grp_p = DataFrame()
    grp_50 = DataFrame()
    for i in uni
        dfM = filter(x -> x[col] == i, grp_prop)
        dfM[!, :prop_cate] = dfM[:, :count]/ sum(dfM[:,:count])*100
        grp_p = vcat(grp_p, dfM)

        dfM = filter(x -> x[col] == i, df_fil)
        cond = dfM[:, :lead_time] .< 0
        prop_50 = nrow(dfM[cond, :])/nrow(dfM)*100
        grp_50 = vcat(grp_50, DataFrame(col=>i, :prop_50=>prop_50))
    end
    grp_p = sort(grp_p, col)
    grp_50 = sort(grp_50, col)
    df_ind_per = DataFrame(col => df_fil[:, col] |> unique )
    #grp_p = leftjoin(grp_p, df_ind_per, on=col)
    grp_p_unstack = unstack(grp_p, col, :lead_time_category, :prop_cate)
    # Reorder the heatmap
    order = bin_labels[end:-1:begin]
    grp_50[:, :prop_50] .= 100 .- grp_50[:, :prop_50]
    y = @pipe grp_p_unstack[:, order] |> Matrix .|> coalesce(_, 0.0)
    return (y, grp_50)
end

function discretise_balance_color(bin_labels::Vector)
    colors = cgrad(:balance, length(bin_labels), categorical = true, rev=false)
    colors = reshape([colors[i] for i in 1:length(colors)], 1,:)
    colors
end

function visualise_heatmap_include_only(x, y, grp_50, bin_labels)
    colors = discretise_balance_color(bin_labels)
    pl = areaplot(x, y, size=(800,600), label=nothing, fillalpha = [0.85 0.85], color=colors)
    plot!(x, grp_50[:, :prop_50], color="black", label=nothing, lw=1.5,
        marker=:circle, markersize=3, markerstrokewidth = 0,
    )
    #ls=:dashdot)
    pl
end

function add_reverse_order_legend!(pl::Plots.Plot, bin_labels::Vector)
    # for color adjusting.
    colors = discretise_balance_color(bin_labels)
    areaplot!(pl, 1, reshape([0 for i in 1:length(bin_labels)], 1, :),
        label=reshape(bin_labels, 1, :),
        color=reverse(colors))
end

function df_to_heatmap(df_res::DataFrame, x, col;
        xticks=:none, add_zero::Bool=false, pc=1.0
        )
    df_fil, bin_labels = create_lead_time_category(df_res)
    colors = discretise_balance_color(bin_labels)
    y, grp_50 = proportion_each_cate_by_group(df_fil, col)
    if add_zero == true
        n_cate = size(y)[2]
        y0 = hcat(fill(0, 1, n_cate-1), [100])
        y = vcat(y0, y)
        grp_50 = vcat(DataFrame(Dict(col=>0, :prop_50 =>0)), grp_50)
        x = vcat([0], x)
    end

    pl = visualise_heatmap_include_only(x, y, grp_50, bin_labels)
    ticks= [0, 20, 40 ,60, 80, 100]
    hline!(pl, ticks, color=:white, alpha=0.5, ls=:dash, label=:none)
    if (col == :ind_site)  & (maximum(x) <= 100)
        ticks_scaled = Int64.(pc.*ticks)
        plot!(pl, xticks=(ticks, ticks_scaled), yticks=(ticks, ticks),
            xlim=[0,100], ylim=[0,100],
        )
        vline!(pl, ticks, color=:white, alpha=0.5, ls=:dash, label=:none)
    elseif xticks != :none
        plot!(pl, xticks=(xticks, xticks), yticks=(ticks, ticks),
            ylim=[0, 100],
        )
        vline!(pl, xticks, color=:white, alpha=0.5, ls=:dash, label=:none)
    else
        xticks = [0,20,40,60,80]
        plot!(pl, xticks=(xticks, xticks), yticks=(ticks, ticks),
            ylim=[0, 100],
        )
        vline!(pl, xticks, color=:white, alpha=0.5, ls=:dash, label=:none)
    end
    return pl
end

@doc (@doc create_lead_time_category) proportion_each_cate_by_group,
    discretise_balance_color,
    visualise_heatmap_include_only,
    create_reverse_order_legend!,
    df_to_heatmap

function single_stacked_heatmap(path_res;
        x_var="coverage", xlim=:none, ylim=:none,
        legend=true, vis_kwds=(),
        xlabel="", ylabel="",
    )
    @unpack ES_pattern, inc_prop, sim_res, path_trans, pars =
        load(path_res)

    sp_pars = read_spatial_params_file(ES_pattern)
    x = obtain_x_values_for_figure(sp_pars, inc_prop, x_var)

    df_fil, bin_labels = create_lead_time_category(sim_res)

    if length(vis_kwds) == 0
        vis_kwds = (
            left_margin=5Plots.mm, right_margin=40Plots.mm,
            xlabelfontsize=12, ylabelfontsize=12, tickfontsize=12,
        )
    end

    pl = df_to_heatmap(sim_res, x, :ind_site; add_zero=true, pc=pars.pc)
    plot!(pl,
        xlabel=xlabel,
        ylabel=ylabel, # "Probability of each pattern (%)",
        xlim=xlim
        ;vis_kwds...
    )
    if legend == true
        add_reverse_order_legend!(pl, bin_labels)
        plot!(pl, legend=(1.13, 0.9))
    end

    return pl
end

"""From simulated results, extract the early detection probabilities
over ES sensitivity analysis.
"""
function fetch_early_det_50(path_res)::DataFrame
    @unpack sim_res = load(path_res)
    df_fil, bin_labels = create_lead_time_category(sim_res)
    y, grp_50 = proportion_each_cate_by_group(df_fil, :ind_site)
    return grp_50
end


"""
    obtain_x_values_for_figure(sp_pars, inc_prop, x_var)

Args:
- `x_var`: Taks `coverage` or `site`.
"""
function obtain_x_values_for_figure(sp_pars, inc_prop, x_var)
    per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
    sens_index = obtain_ES_sensitivity_index(sp_pars.pop, inc_prop)
    if x_var == "coverage"
        return per_pop[sens_index]
    elseif x_var == "site"
        return sens_index
    else
        error("Specify: coverage or site")
    end
end

function plot_sens_adaptor(path_lis, single_vis!;
        xlabel="Number of ES covered site",
        ylabel="Early detection probability (%)",
        legendtitle="", labels=:none,
        xlabelfontsize=12, ylabelfontsize=12,
        tickfontsize=12, legend=:topleft,
        lw=:none, ls=:none, color=:none,
    )
    ticks= [0, 20, 40 ,60, 80]
    pl = plot(
        xlabel=xlabel, ylabel=ylabel,
        legendtitle = legendtitle,
        legendtitlefontsize = 8,
        legend =legend,
        xlabelfontsize=xlabelfontsize,
        ylabelfontsize=ylabelfontsize,
        tickfontsize=tickfontsize,
        xticks=(ticks, ticks),
        foreground_color_legend=nothing,
        background_color_legend=nothing,
    )

    for (i, path) in enumerate(path_lis)
        label = labels == :none ? :none : labels[i]
        vis_kwds = (lw=lw[i], ls=ls[i], color=color[i])
        if labels != :none
            vis_kwds = merge(vis_kwds, (label=labels[i],))
        end
        single_vis!(pl, path;
            vis_kwds=vis_kwds
        )
    end
    return pl
end

function single_vis_pc_site!(pl, path; vis_kwds=:none)
    x_var="site"
    grp_50 = fetch_early_det_50(path)
    @unpack pars, inc_prop, ES_pattern, par_ES = load(path)
    sp_pars = read_spatial_params_file(ES_pattern)
    x = obtain_x_values_for_figure(sp_pars, inc_prop, x_var)

    plot!(pl, x, grp_50[:, :prop_50],
        xlim=[0,90], label=f"{par_ES.pc:.2f}";
        vis_kwds...
    )
end

function single_vis_pc_coverage!(pl, path; vis_kwds=:none)
    x_var="coverage"
    grp_50 = fetch_early_det_50(path)
    @unpack pars, inc_prop, ES_pattern, par_ES = load(path)
    sp_pars = read_spatial_params_file(ES_pattern)
    x = obtain_x_values_for_figure(sp_pars, inc_prop, x_var)

    plot!(pl, x*par_ES.pc, grp_50[:, :prop_50],
        xlim=[0,100], label=f"{par_ES.pc:.2f}";
        vis_kwds...
    )

end

function single_vis_R0!(pl, path;
        vis_kwds=:none
    )
    x_var="site"
    grp_50 = fetch_early_det_50(path)
    @unpack pars, inc_prop, ES_pattern = load(path)
    sp_pars = read_spatial_params_file(ES_pattern)
    x = obtain_x_values_for_figure(sp_pars, inc_prop, x_var)

    plot!(pl, x, grp_50[:, :prop_50],
        xlim=[0,90], label=f"{pars.R0:d}";
        vis_kwds...
    )
end

function single_vis_α!(pl, path;
        vis_kwds=:none
    )
    x_var="site"
    grp_50 = fetch_early_det_50(path)
    @unpack pars, inc_prop, ES_pattern = load(path)
    sp_pars = read_spatial_params_file(ES_pattern)
    x = obtain_x_values_for_figure(sp_pars, inc_prop, x_var)

    plot!(pl, x, grp_50[:, :prop_50],
        xlim=[0,90], label=f"{pars.α:.3f}";
        vis_kwds...
    )
end

function single_vis_sampling_freq!(pl, path;
        vis_kwds=:none,
    )
    x_var="site"
    grp_50 = fetch_early_det_50(path)

    @unpack pars, inc_prop, ES_pattern, par_ES = load(path)
    sp_pars = read_spatial_params_file(ES_pattern)
    x = obtain_x_values_for_figure(sp_pars, inc_prop, x_var)

    plot!(pl, x, grp_50[:, :prop_50],
        xlim=[0,90], label=f"{par_ES.n_freq:d} day";
        vis_kwds...
    )
end

function single_vis_ES_det!(pl, path;
        vis_kwds=:none,
    )
    x_var="site"
    grp_50 = fetch_early_det_50(path)

    @unpack pars, inc_prop, ES_pattern, par_ES = load(path)
    sp_pars = read_spatial_params_file(ES_pattern)
    x = obtain_x_values_for_figure(sp_pars, inc_prop, x_var)

    plot!(pl, x, grp_50[:, :prop_50],
        xlim=[0,90];
        vis_kwds...
    )
end

function check_single_percentage(path_res)
    @unpack sim_res = load(path_res)
    df_res = sim_res # ES population coverage

    # Summarise the proportion
    dfM = filter(x -> x.ind_site .== 31, df_res)
    dfM_fil, bin_labels = create_lead_time_category(dfM)
    tab1 = countmap(dfM_fil[:, :lead_time_category])
    ks = keys(tab1) .|> String
    vs = values(tab1) |> collect
    tab = DataFrame(label=ks, value=vs)
    tab[:, :prop] = tab[:, :value]/sum(tab[:,:value])*100
    n_nan = size(dfM)[1] - sum(tab[:,:value])
    println("No detection simulations: $(n_nan)")
    return tab
end

function calculate_d_mat(mat::Matrix)::Matrix
    n_point = size(mat)[1]
    d_mat = fill(0., n_point, n_point)
    for i in 1:n_point, j in 1:n_point
        ϕ1, λ1 = mat[i, :]
        ϕ2, λ2 = mat[j, :]
        d_mat[i,j] = harversine_dist(ϕ1, λ1, ϕ2, λ2)
    end
    return d_mat
end

"""

Args:
- `p_imp`: Importation probability vectors.
"""
function calculate_d_ave_over_index(
    df_zaf, p_imp::Vector{Float64}
)
    # to ensure the order is the same as other files.
    #sort!(df_zaf, "value", rev=true)
    pop = df_zaf[:, :value]
    mat = df_zaf[:, [:lat, :lon]] |> Matrix
    d_mat = calculate_d_mat(mat)

    sens_index = obtain_ES_sensitivity_index(pop, 0.01)
    # Calculate the average dist.
    d_sens = []
    for ind_col in sens_index
        d_m = minimum(d_mat[:, begin:ind_col] , dims=2)[:,1]
        d_ave = sum(p_imp .* d_m)
        push!(d_sens, d_ave)
    end
    return d_sens
end

function calculate_distance_and_prop_50(df_zaf, p_imp, path_res)
    d_sens = calculate_d_ave_over_index(df_zaf, p_imp)
    prop_50 = fetch_early_det_50(path_res)
    return d_sens, prop_50
end