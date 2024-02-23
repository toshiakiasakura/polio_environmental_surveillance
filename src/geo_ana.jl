using ArchGDAL
using GeoDataFrames
import GeoDataFrames as GDF
using Rasters
using Shapefile

include("utils.jl")

"""
    harvesine_dist(ϕ1, λ1, ϕ2, λ2)::Float64

Calculate distance between two geographical points using Vincenty Formula.

# Arguments
- `ϕ1::Float64`: latitude 1 [deg].
- `λ1::Float64`: longitude 1 [deg].
- `ϕ2::Float64`: latitude 2 [deg].
- `λ2::Float64`: longitude 2 [deg].
"""
function harversine_dist(ϕ1::Float64, λ1::Float64, ϕ2::Float64, λ2::Float64)::Float64
    r = 6371
    Δϕ = ϕ1 - ϕ2
    Δλ = λ1 - λ2
    a = sind(Δϕ / 2.0)^2 + cosd(ϕ1) * cosd(ϕ2) * sind(Δλ / 2.0)^2
    d = 2 * r * atand(√a, √(1.0 - a)) * 2π / 360 # convert degree to radian
    return d
end

function get_gridsize(ras::Raster)
    lon = lookup(ras, X)
    lat = lookup(ras, Y)
    Δlat = @pipe harversine_dist(lat[1], lon[1], lat[2], lon[1]) |> round(_; digits=2)
    Δlon = @pipe harversine_dist(lat[1], lon[1], lat[1], lon[2]) |> round(_; digits=2)
    println("Latitude diff: $Δlat km, Longitude diff: $Δlon km")
end

function raster_to_df(ras::Raster)::DataFrame
    coords = []
    lon = lookup(ras, X)
    lat = lookup(ras, Y)
    for lo in lon, la in lat
        push!(coords, (lon=lo, lat=la, value=ras[At(lo), At(la)]))
    end
    df = DataFrame(coords)
    return df
end

function borders!(p, poly)
    for i in 1:length(p)
        plot!(p, poly; subplot=i, fillalpha=0, linewidth=0.6)
    end
    return p
end

function add_zaf_borders!(pl)
    path_shape = "../dt_geoBoundaries-ZAF-ADM2-all/geoBoundaries-ZAF-ADM2.shp"
    shapes = Shapefile.Handle(path_shape)
    borders!(pl, shapes)
end

function add_zaf_ADM0_borders!(pl)
    path_shape = "../dt_geoBoundaries-ZAF-ADM0-all/geoBoundaries-ZAF-ADM0.shp"
    shapes = Shapefile.Handle(path_shape)
    borders!(pl, shapes)
end

function add_moz_borders!(pl)
    path_shape = "../dt_geoBoundaries-MOZ-ADM0-all/geoBoundaries-MOZ-ADM0.shp"
    shapes = Shapefile.Handle(path_shape)
    borders!(pl, shapes)
end

"""

# Arguments
- `pop`: Population size vector.
- `mat`: Matrix of latitude and longitude.

# Returns
- Matrix: Movement of travellers from columns to rows.
"""
function calculate_probability_of_pi(pop::Vector, mat::Matrix)::Matrix
    n_point = length(pop)
    d_mat = fill(0.0, n_point, n_point)
    for i in 1:n_point, j in 1:n_point
        ϕ1, λ1 = mat[i, :]
        ϕ2, λ2 = mat[j, :]
        d_mat[i, j] = harversine_dist(ϕ1, λ1, ϕ2, λ2)
    end

    # Calculate sij
    s_mat = fill(0.0, n_point, n_point)
    for i in 1:n_point, j in 1:n_point
        d_i = d_mat[i, :]
        cond = d_i .<= d_i[j]
        # exclude the source and destination area.
        s_mat[i, j] = sum(pop[cond]) - pop[i] - pop[j]
    end
    s_mat[diagind(s_mat)] .= 0

    # Calculate πij
    π_mat = fill(0.0, n_point, n_point)
    for i in 1:n_point, j in 1:n_point
        π_mat[i, j] = pop[i] * pop[j] / (pop[i] + s_mat[i, j]) / (pop[i] + pop[j] + s_mat[i, j])
    end
    π_mat[diagind(π_mat)] .= 0
    return π_mat
end

function create_vaccination_coverage_map(df_geo::DataFrame)
    l = @layout [a{0.95w} b]
    pl = plot(axis=nothing, border=:none, dpi=300)
    lower = 0.75
    cm = reverse(ColorSchemes.matter)
    for r in eachrow(df_geo)
        sc = (r.EVP - lower) / (1 - lower)
        plot!(pl, r.geometry, color=cm[sc])
    end

    # International airport stars
    #scatter!(pl, [28.2420, 18.6002, 31.1154], [-26.1282, -33.9705, -29.6087], marker=:star,
    #    markersize=6, color=["yellow", "lawngreen", "darkslategray1"], label=:none,
    #    markerstrokewidth=0.1)
    pl2 = heatmap(rand(2, 2), clims=(lower * 100, 100), framestyle=:none,
        c=cgrad(cm), cbar=true, lims=(-1, 0),
    )
    pl = plot(pl, pl2, layout=l, right_margin=10Plots.mm, fmt=:png)
    display(pl)
    savefig(pl, "../res/fig_vaccine_coverage.png")
end

function save_vaccination_coverage_data(df_vac)
    df_save = copy(df_vac)
    cols = ["shapeName", "sample_size", "OPV0", "OPV1",
        "HEXA1", "HEXA2", "HEXA3", "HEXA4", "EVP",]
    covs = ["OPV0", "OPV1", "HEXA1", "HEXA2", "HEXA3", "HEXA4", "EVP",]
    df_save = df_save[:, cols]
    df_save[:, covs] = df_save[:, covs] .* 100
    df_save[:, "EVP"] = round.(df_save[:, "EVP"], digits=3)
    CSV.write("../res/table_vaccine_coverage.csv", df_save)
end

function read_agg_map(path; agg_scale=230)::Raster
    map_tmp = read(Raster(path))
    map_tmp = replace_missing(map_tmp, 0)
    map_agg = Rasters.aggregate(sum, map_tmp, agg_scale; skipmissingval=true)
    return map_agg
end

function print_basic_map_info(map::Raster)::Nothing
    nr, nc = size(map)
    println("row: $nr, col: $nc, grid num: $(nr*nc)")
    get_gridsize(map)
    map |> sum |> println
end

function merge_two_map_data(path1, path2; agg_scale=230)
    #agg_scale = 110 # Latitude diff: 10.19 km, Longitude diff: 9.44 km
    #agg_scale = 230 # Latitude diff: 21.31 km, Longitude diff: 19.74 km
    f_5_zaf = read(Raster(path1))
    f_5_zaf = replace_missing(f_5_zaf, 0)
    plot(f_5_zaf) |> display

    m_5_zaf = read(Raster(path2))
    m_5_zaf = replace_missing(m_5_zaf, 0)
    plot(m_5_zaf) |> display

    m_5_zaf_agg = Rasters.aggregate(sum, m_5_zaf, agg_scale; skipmissingval=true)
    f_5_zaf_agg = Rasters.aggregate(sum, f_5_zaf, agg_scale; skipmissingval=true)
    zaf_5_agg = m_5_zaf_agg .+ f_5_zaf_agg

    # Print basic info.
    nr, nc = size(zaf_5_agg)
    println("row: $nr, col: $nc, grid num: $(nr*nc)")
    get_gridsize(zaf_5_agg)
    # Check population size change.
    m_5_zaf |> sum |> println
    f_5_zaf |> sum |> println
    zaf_5_agg |> sum |> println

    return zaf_5_agg
end

function cut_validate_raster_dataframe!(df_ras::DataFrame, cut_off_pop)
    println("Original size: ", size(df_ras))
    df_ras[!, :value] = @pipe df_ras[:, :value] .|> round(_; digits=0)
    filter!(x -> x.value > 0.0, df_ras)

    n_bf = df_ras[:, :value] |> sum
    println("Total population size before removing: ", n_bf)
    println("Before removing: ", size(df_ras))
    filter!(x -> x.value > cut_off_pop, df_ras)
    n_af = df_ras[:, :value] |> sum
    println("Total population size after removing: ", n_af)
    println("After removing: ", size(df_ras))
    println("% or removal: ", (n_bf - n_af) / n_bf * 100)
end

function calculate_minimum_dist_given_polygon(lon, lat, pol)
    pol_sim = ArchGDAL.simplify(pol, 0.05)
    boundary = ArchGDAL.boundary(pol_sim)
    n = ArchGDAL.ngeom(boundary)
    dist = Inf
    for i in 1:n
        lon2, lat2, z = ArchGDAL.getpoint(boundary, i - 1)
        dist_tmp = harversine_dist(lat, lon, lat2, lon2)
        dist = minimum([dist, dist_tmp])
    end
    return dist
end

function relate_df_ras_to_district_info(df_ras::DataFrame, df_geo::DataFrame)
    dist = []
    for r in eachrow(df_ras)
        point = ArchGDAL.createpoint(r.lon, r.lat)
        flag = false
        for g in eachrow(df_geo)
            pol = g.geometry
            if ArchGDAL.within(point, pol)
                push!(dist, g.shapeName)
                flag = true
                break
            end
        end
        if flag == false
            push!(dist, "not determined")
        end
    end
    df_ras[!, "shapeName"] = dist

    # some of points are not classified.
    n = (df_ras[:, "shapeName"] .== "not determined") |> sum
    println("Number of unclassified points: ", n)

    # Check those points
    df_fil = filter(x -> x.shapeName == "not determined", df_ras)
    pl = plot()
    scatter!(pl, df_fil[:, :lon], df_fil[:, :lat],)
    add_zaf_borders!(pl)
    display(pl)

    # Classify those points to nearest districts.
    for (i, r) in enumerate(eachrow(df_ras))
        if r.shapeName == "not determined"
            dist = Inf
            dist_name = "not determined"
            for r_geo in eachrow(df_geo)
                dist_tmp = calculate_minimum_dist_given_polygon(r.lon, r.lat, r_geo.geometry)
                if dist > dist_tmp
                    dist = dist_tmp
                    dist_name = r_geo.shapeName
                end
            end
            df_ras[i, "shapeName"] = dist_name
        end
    end
    return df_ras
end

function visualise_population(df_mer::DataFrame, col::Symbol,
    map, title::String
)
    map = copy(map)
    map[:] .= 0
    n = size(df_mer)[1]
    for i in 1:n
        x = df_mer[i, :lon]
        y = df_mer[i, :lat]
        v = df_mer[i, col] #|> Int64
        map[At(x), At(y)] = v
    end
    pl = plot()
    plot!(map,
        xlim=[15, 35], ylim=[-35, -21],
        axis=nothing, border=:none,
        right_margins=9Plots.mm,
        #colorbar_title="log10(πij)",
        title=title,
    )
    add_zaf_borders!(pl)
    return pl
end