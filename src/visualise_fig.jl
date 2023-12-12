function check_single_percentage(path_res)
    path_params, res_all = deserialize(path_res)
    df_res1, df_res2, df_res3 = res_all
    df_res = df_res3 # ES population coverage 
    
    # Summarise the proportion
    dfM = filter(x -> x.ind_site .== 33, df_res)
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

function single_figure(path_spatial, path_res; 
        x_var="coverage", xlim=:none, ylim=:none,
    )
    sp_pars = deserialize(path_spatial)
    per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
    sens_index = obtain_ES_sensitivity_index(sp_pars.pop, 0.01)
    x_per_pop =  per_pop[sens_index]
    
    path_params, res_all = deserialize(path_res)
    df_res1, df_res2, df_res3 = res_all
    df_res = df_res3 # ES population coverage 
    df_fil, bin_labels = create_lead_time_category(df_res)
    
    # Visualise.
    paths = fetch_sim_paths(path_params)
    res = deserialize(paths[1])
    pc = res.pars.pc

    if x_var == "coverage"
        x = x_per_pop
        xlabel = "ES population coverage (%)"
    elseif x_var == "site"
        x = sens_index
        xlabel = "Number of ES covered sites"
    end
    
    pl = df_to_heatmap(df_res, x, :ind_site; add_zero=true, pc=pc)
    add_reverse_order_legend!(pl, bin_labels)
    plot!(pl, 
        left_margin=5Plots.mm, right_margin=40Plots.mm, 
        xlabel=xlabel,
        ylabel="Probability of each pattern (%)", 
        xlabelfontsize=14, ylabelfontsize=14, tickfontsize=12,
        xlim=xlim, 
    )
    plot!(pl, legend=(1.1, 0.9), dpi=300, fmt=:png)
    display(pl)
end

function three_scenario_results(path_spatial, path_res1, path_res2, path_res3;
        x_var="coverage", xlim=:none, 
        airport_order="population",
    )

    sp_pars = deserialize(path_spatial)
    per_pop = cumsum(sp_pars.pop)/sum(sp_pars.pop)*100
    sens_index = obtain_ES_sensitivity_index(sp_pars.pop, 0.01)
    x_per_pop =  per_pop[sens_index]
    if x_var == "coverage"
        x = x_per_pop
        xlabel = "ES population coverage (%)"
    elseif x_var == "site"
        x = sens_index
        xlabel = "Number of ES covered sites"
    end

    path_params1, res_all1 = deserialize(path_res1)
    path_params2, res_all2 = deserialize(path_res2)
    path_params3, res_all3 = deserialize(path_res3)
    
    paths = fetch_sim_paths(path_params1)
    res = deserialize(paths[1])
    pc = res.pars.pc
    
    pl1 = df_to_heatmap(res_all1[3], x, :ind_site; 
        add_zero=true, pc=pc)
    plot!(pl1, 
        xlabel=xlabel,
        ylabel="Probability of each pattern (%)", 
        title="Population size scenario",
        left_margin=5Plots.mm, 
        xlim=xlim,
    )
    pl2 = df_to_heatmap(res_all2[3], x, :ind_site; 
        add_zero=true, pc=pc)
    plot!(pl2, 
        xlabel=xlabel,
        title="Airport scenario", 
        tmargin=50Plots.mm,
        xlim=xlim,
    )
    if airport_order == "population"
        airport_cov = per_pop[[11, 7, 62]]
    elseif airport_order == "mozambique"
        airport_cov = per_pop[[187,65, 435]] 
    end
    for cov in airport_cov
        continue
        #annotate!(pl2, cov, 100, text("↓", :bottom, 20, :black))
    end
    pl3 = df_to_heatmap(res_all3[3], x, :ind_site; 
        add_zero=true, pc=pc)
    plot!(pl3, 
        xlabel=xlabel,
        #right_margin=40Plots.mm, 
        ylabel="Probability of each pattern (%)", 
        title="Mozambique scenario",
        xlim=xlim,
    )
    pl4 = plot(showaxis = false, foreground_color_grid=:white,
        legend=(0.1, 0.9),
    )
    df_fil, bin_labels = create_lead_time_category(res_all3[3])
    add_reverse_order_legend!(pl4, bin_labels)
    pls = [pl1, pl2, pl3, pl4]
    l = @layout [a b; c d]
    pl = plot(pls..., 
        layout=l, dpi=300,
        bottom_margin=5Plots.mm,
        xtickfontsize=10, ytickfontsize=10,
    )
    plot!(pl, size=(1000, 700), fmt=:png)
end

function calculate_d_mat(pop, mat)
    n_point = length(pop)
    d_mat = fill(0., n_point, n_point)
    for i in 1:n_point, j in 1:n_point
        ϕ1, λ1 = mat[i, :]
        ϕ2, λ2 = mat[j, :]
        d_mat[i,j] = harversine_dist(ϕ1, λ1, ϕ2, λ2)
    end
    return d_mat
end

function calculate_d_ave_over_index(df_zaf, p_imp)
    # to ensure the order is the same as other files.
    #sort!(df_zaf, "value", rev=true) 
    pop = df_zaf[:, :value]
    mat = df_zaf[:, [:lat, :lon]] |> Matrix
    d_mat = calculate_d_mat(pop, mat)

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

"""From simulated results, extract the early detection probabilities
over ES sensitivity analysis.
"""
function fetch_early_det_50(path)
    path_params, res_all1 = deserialize(path)
    df_res1, df_res2, df_res3 = res_all1
    df_res = df_res3
    df_fil, bin_labels = create_lead_time_category(df_res)
    y, grp_50 = proportion_each_cate_by_group(df_fil, :ind_site) 
    prop_50 = 100 .- grp_50[:, :prop_0]
    return prop_50
end

function calculate_distance_and_prop_50(df_zaf, p_imp, path_res)
    d_sens = calculate_d_ave_over_index(df_zaf, p_imp)
    prop_50 = fetch_early_det_50(path_res)
    return d_sens, prop_50
end