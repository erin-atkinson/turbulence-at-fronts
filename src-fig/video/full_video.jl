function full_video(
        foldername,
        filename,
        frames,
        z;
        fig_kw=(; ), 
        ax_kw=(; ),
        axh_kw=(; ),
        ax_series_kw=(; ),
        ht_u_kw=(; ),
        ht_Vq_kw=(; ),
        ln_Vq_kw=(; ),
        ln_VSP_kw=(; ),
        ct_b_kw=(; ),
        ct_v_kw=(; ),
        record_kw=(; ),
        σ=0,
        σh=0,
        background=false
    )
    full_iterations, full_times = iterations_times(foldername)
    sp = simulation_parameters(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    inds = centre_indices(foldername)
    colormap_u = to_colormap(:balance)
    colormap_Vq = to_colormap(:magma)
    z_indᶜ = zᶜbounds(foldername, z)

    iterations = full_iterations[frames]
    times = full_times[frames]

    n = Observable(1)
    frame = @lift frames[$n]
    
    iteration = @lift iterations[$n]
    t = @lift times[$n]
    
    u_title = @lift let time_string = prettytime($t / 3600; digits=1)
        L"t=%$time_string ~\text{hr}"
    end
    
    U = [-sp.α * x for x in xsᶠ, y in 1:1] .* background
    V = [sp.α * y for x in 1:1, y in ysᶠ] .* background
    
    fig = Figure(; 
        size=(960, 540),
        fig_kw...
    )
    
    DFM = jldopen(joinpath(foldername, "DFM.jld2"))
    PV = jldopen(joinpath(foldername, "PV.jld2"))
    OUTPUT = jldopen(joinpath(foldername, "output.jld2"))
    
    colorrange_u = (-0.1, 0.1)
    colorrange_Vq = (0, 1)
    v_levels = range(-0.5, 0.5, 0.05)

    u = @lift get_field(DFM, "u_dfm", $iteration) .+ U
    v = @lift get_field(DFM, "v_dfm", $iteration)
    b = @lift get_field(DFM, "b_dfm", $iteration)
    
    uh = @lift get_field(a->a[:, :, z_indᶜ], OUTPUT, "u", $iteration) .+ U
    vh = @lift get_field(a->a[:, :, z_indᶜ], OUTPUT, "v", $iteration) .+ V
    bh = @lift get_field(a->a[:, :, z_indᶜ], OUTPUT, "b", $iteration)

    Vq = @lift get_field(PV, "Vq_dfm", $iteration)
    
    Vq_series = timeseries_of(joinpath(foldername, "PV.jld2"), "Vq_dfm", full_iterations) do field
        mean(field[inds, :]) * sp.Lh * sp.Ly * sp.Lz
    end

    VSP_series = timeseries_of(joinpath(foldername, "TKE.jld2"), "VSP", full_iterations) do field
        mean(field[inds, :]) * sp.Lh * sp.Ly * sp.Lz
    end

    Vq_point = @lift [Point2D($t / 3600, Vq_series[$frame])]
    VSP_point = @lift [Point2D($t / 3600, VSP_series[$frame])]
    
    ax_u = Axis(fig[2, 1];
        limits=(-sp.Lh/2000, sp.Lh/2000, -sp.H, 0),
        xlabel=L"x / \text{km}",
        ylabel=L"z / \text{m}",
        xticks=[-2, -1, 0, 1, 2],
        ax_kw...
    )
    ax_uh = Axis(fig[1, 1];
        limits=(-sp.Lh/2000, sp.Lh/2000, -sp.Ly/2000, sp.Ly/2000),
        xlabel=L"x / \text{km}",
        ylabel=L"y / \text{km}",
        title=u_title,
        axh_kw...
    )

    ax_Vq = Axis(fig[2, 2];
        limits=(-sp.Lh/2000, sp.Lh/2000, -sp.H, 0),
        xlabel=L"x / \text{km}",
        ylabel=L"z / \text{m}",
        xticks=[-2, -1, 0, 1, 2],
        ax_kw...
    )

    ax_Vq_series = Axis(fig[1, 2];
        limits=(0, full_times[end] / 3600, nothing, nothing),
        xlabel=L"t / \text{hr}",
        ylabel=L"Vq",
        ax_series_kw...
    )
    ax_VSP_series = Axis(fig[1, 2];
        limits=(0, full_times[end] / 3600, nothing, nothing),
        xlabel=L"t / \text{hr}",
        ylabel=L"VSP / \text{kW}",
        yaxisposition = :right,
        ax_series_kw...
    )
    
    ht_u_kw = (;
        colorrange=colorrange_u,
        lowclip=colormap_u[1],
        highclip=colormap_u[end],
        colormap,
        ht_u_kw...
    )
    ct_b_kw = (;
        color=(:black, 0.5),
        levels=b_levels,
        ct_b_kw...
    )
    ct_v_kw = (;
        color=(:blue, 0.5),
        levels=v_levels,
        ct_v_kw...
    )

    ht_Vq_kw = (;
        colorrange=colorrange_Vq,
        colormap,
        ht_Vq_kw...
    )
    ln_Vq_kw = (;
        color=:blue,
        ln_Vq_kw...
    )
    ln_VSP_kw = (;
        color=:green,
        ln_VSP_kw...
    )

    # Mean and slice of u, v and b contours
    ht_u = heatmap!(ax_u, xsᶠ ./ 1000, zsᶜ, u; ht_u_kw...)
    ht_uh = heatmap!(ax_uh, xsᶠ ./ 1000, ysᶜ / 1000, uh; ht_u_kw...)
    contour!(ax_u, xsᶜ ./ 1000, zsᶜ, b; ct_b_kw...)
    contour!(ax_uh, xsᶜ ./ 1000, ysᶜ ./ 1000, bh; ct_b_kw...)
    contour!(ax_u, xsᶜ ./ 1000, zsᶜ, v; ct_v_kw...)
    contour!(ax_uh, xsᶜ ./ 1000, ysᶜ ./ 1000, vh; ct_v_kw...)

    # PV and VSP lines and point
    ln_Vq = lines(ax_Vq_series, full_times ./ 3600, Vq_series; ln_Vq_kw)
    ln_VSP = lines(ax_VSP_series, full_times ./ 3600, VSP_series; ln_VSP_kw)
    scatter!(ax_Vq_series, Vq_point; color=:blue, marker=:+)
    scatter!(ax_VSP_series, VSP_point; color=:green, marker=:+)

    # PV heatmap
    ht_Vq = heatmap!(ax_Vq, xsᶜ ./ 1000, zsᶜ, Vq; ht_Vq_kw...)

    Colorbar(fig[3, 1], ht_u; vertical=false, flipaxis=false, label=L"u / \text{ms}^{-1}")
    Colorbar(fig[3, 2], ht_Vq; vertical=false, flipaxis=false, label=L"V_{q < 0}")

    hidespines!(ax_VSP_series)
    hidexdecorations!(ax_VSP_series)

    hidexdecorations!(ax_uh; ticks=true)
    hideydecorations!(ax_Vq; ticks=false)

    rowgap!(fig.layout, 1, Relative(0.02))

    subfig_label!(fig[1, 1], 1)
    subfig_label!(fig[1, 2], 2)
    subfig_label!(fig[2, 1], 3)
    subfig_label!(fig[2, 2], 4)
    
    lines!(ax_u, [-sp.Lx / 2, sp.Lx / 2], [z, z]; color=(:red, 0.5), linestyle=:dash)

    if length(frames) > 1
        record(fig, filename, 1:length(frames); record_kw...) do i
            n[] = i
            print("$i / $(length(frames)) \r")
        end
        println("")
    end
    close(DFM)
    close(OUTPUT)
    close(PV)
    close(TKE)

    fig
end