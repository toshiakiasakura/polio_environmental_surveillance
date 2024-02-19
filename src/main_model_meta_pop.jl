include("utils.jl")
include("model_meta_pop.jl")

function baseline_results(
    par_AFP, par_ES;
    path_trans=""
    )
    n_sim = fetch_sim_paths(path_trans) |> length |> println
    paths = fetch_sim_paths(path_trans)
    res = deserialize(paths[1])

    sim_res = @time collect_summary_statistics(path_trans, par_AFP, par_ES)

    # Visualisation.
    df = DataFrame(sim_res)
    n_sim = size(df)[1]
    #h1 = histogram(df[:, :R_final_num], bins=20, title="Final size")
    h2 = histogram(df[:, :R_final_site], bins=20,
        xlabel="Number of sites having >=1 infections", ylabel="Density",
        norm=true, legend=:none
    )
    h3 = histogram(df[:, :R_final_AFP], bins=20,
        xlabel="Number of final AFP cases",
        norm=true, legend=:none,
        ylabel="Density"
    )

    max_f = maximum(df[: ,:R_final_AFP])
    df_fil = filter(x -> x.R_final_AFP != 0, df)
    histogram!(h3, df_fil[:, :R_final_AFP], bins=20, legend=:none,
        inset = bbox(0.5Plots.w, 0.1Plots.h, 0.5Plots.w, 0.5Plots.h),
        subplot=2,
        xlabel="Number of final AFP cases excluding 0",
        ylabel="Density",
        norm=true,
    )
    l = @layout [a{0.4w} b]
    plot(h2, h3,  size=(900,400), fmt=:png, dpi=300,
        left_margin=5Plots.mm,
        bottom_margin=5Plots.mm,
        layout=l,
    ) |> display

    days = res.pars.days
    ts_ES = df[:, "t_ES"]
    ts_AFP = df[:, "t_AFP"]
    t_extinct = df[:, "t_extinct"]


    pl1 = plot(
        xlabel="Day", ylabel="Cumulative probability of first detection",
        legend=(0.6, 0.7),
    )
    annotate!((0.07, 0.95), "(A)")
    cum_ES = cumulative_counts(ts_ES, days; prop=true)
    cum_AFP = cumulative_counts(ts_AFP, days; prop=true)
    cum_ES_cond = conditional_cumulative_prob(ts_ES, t_extinct, days)
    cum_AFP_cond = conditional_cumulative_prob(ts_AFP, t_extinct, days)

    plot!(pl1, 1:days, cum_ES, label="Prob. via ES", color=1)
    plot!(pl1, 1:days, cum_ES_cond, label="Conditional Prob. via ES", color=1, linestyle=:dashdot)
    plot!(pl1, 1:days, cum_AFP, label="Prob. via AFP surv.", color=2)
    plot!(pl1, 1:days, cum_AFP_cond, label="Conditional Prob. via AFP surv.", color=2, linestyle=:dashdot)

    dif = leadtime_diff(df)
    x = [1 for i in 1:length(dif)]
    pl2 = violin(x, dif, xticks=:none, ylabel="Lead time of ES (day)", legend=:none)
    boxplot!(pl2, x, dif, fillalpha=0.75)
    annotate!((0.15, 0.95), "(B)")
    l = @layout [a{0.75w} b]
    pl = plot(pl1, pl2,
        fmt=:png, dpi=300, layout=l,
        size=(800,500), left_margin=5Plots.mm,
    )
    display(pl)
    savefig(pl, "../res/fig_baseline_cum_lead.png")

    outcome_num_prop(df) |> display
    detect_pattern(df) |> countmap |> display

    dfM = @pipe filter(x -> isnan(x["t_ES"]) == false, df)
    m_ES  = @pipe mean(dfM[:, "t_ES"]) |> round(_, digits=2)
    dfM = @pipe filter(x -> isnan(x["t_AFP"]) == false, df)
    m_AFP = @pipe mean(dfM[:, "t_AFP"])  |> round(_, digits=2)
    println("Mean ES: $m_ES days, Mean AFP: $m_AFP days")

    #vis_cumulative_prob(df, pars.days; title="sim =$n_sim, R0=$(res.pars.R0)")
    dif = leadtime_diff(df)
    leadtime_diff_statistics(dif)
end

function sensitivity_all_summary(
        par_AFP, par_ES;
        path_trans=""
    )
    n_sim = fetch_sim_paths(path_trans) |> length
    res_all = sensitivity_ana_all(path_trans, par_AFP, par_ES)
    now_str = get_today_time()
    path_save = "../dt_tmp_res/$(now_str).ser"
    res_all_path = (path=path_trans, res_all=res_all)
    println(path_save)
    serialize(path_save, res_all_path)
    return path_save
end

"""To validate the implementation of meta-population model,
Calculate the proportion of at least one importation in each grid.
"""
function validate_meta_population_model(path_trans, path_spatial)
    path_objs = fetch_sim_paths(path_trans)
    n_sim = length(path_objs)
    p = path_objs[1]
    res = deserialize(p)
    bool_mat = fill(0, res.rec.n_site)
    @showprogress for (i, p) in enumerate(path_objs)
        res = deserialize(p)
        bool_mat[:] += sum(res.rec.I, dims=2) .>= 1
    end
    prop_mat = bool_mat ./ n_sim .+ 1/n_sim
    sp_pars = deserialize(path_spatial)
    x = log10.(sp_pars.pop)
    y = log10.(prop_mat)
    plot(x, y, xlabel=:log10_pop, ylabel=:log10_prob)
    scatter!(x, y)
end