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
        ln_TKE_kw=(; ),
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
    colormap_Vq = to_colormap(:diff)
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
    
    U = [-variable_strain_rate(t, sp) * x for t in times, x in xsᶠ, y in 1:1, z in 1:1] .* background
    V = [variable_strain_rate(t, sp) * y for t in times, x in 1:1, y in ysᶠ, z in 1:1] .* background
    
    fig = Figure(; 
        size=(960, 540),
        fig_kw...
    )
    
    DFM = jldopen(joinpath(foldername, "DFM.jld2"))
    PV = jldopen(joinpath(foldername, "PV.jld2"))
    OUTPUT = jldopen(joinpath(foldername, "output.jld2"))
    TKE = jldopen(joinpath(foldername, "TKE.jld2"))
    
    colorrange_u = (-6, 6)
    colorrange_Vq = (0, 1)
    v_levels = range(-0.2, 0.2, 12)

    u = @lift 100 * (get_field(DFM, "u_dfm", $iteration) .+ U[$n, :, 1, :])
    v = @lift get_field(DFM, "v_dfm", $iteration)
    b = @lift get_field(DFM, "b_dfm", $iteration)
    MLD = @lift $b .- ($b[:, end:end] .- (b_levels[3] - b_levels[1]))
    
    uh = @lift 100 * (get_field(a->a[:, :, z_indᶜ], OUTPUT, "u", $iteration) .+ U[$n, :, :, 1])
    vh = @lift get_field(a->a[:, :, z_indᶜ], OUTPUT, "v", $iteration) .+ V[$n, :, :, 1]
    bh = @lift get_field(a->a[:, :, z_indᶜ], OUTPUT, "b", $iteration)

    Vq = @lift get_field(PV, "Vq_dfm", $iteration)
    
    Vq_series = timeseries_of(PV, "Vq_dfm", full_iterations) do field
        mean(field[inds, :])
    end

    TKE_series = timeseries_of(TKE, "TKE", full_iterations) do field
        mean(field[inds, :]) * 1e-6 * 1037 * sp.Lh * sp.Ly * sp.Lz
    end

    Vq_point = @lift [Point2($t / 3600, Vq_series[$frame])]
    TKE_point = @lift [Point2($t / 3600, TKE_series[$frame])]
    
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
        limits=(0, full_times[end] / 3600, 0.20, 0.40),
        xlabel=L"t / \text{hr}",
        ylabel="Negative PV fraction",
        yticklabelcolor=:blue,
        ylabelcolor=:blue,
        xaxisposition=:top,
        yticks=range(0.24, 0.40, 3),
        ax_series_kw...
    )
    ax_TKE_series = Axis(fig[1, 2];
        limits=(0, full_times[end] / 3600, 20, 40),
        xlabel=L"t / \text{hr}",
        ylabel=L"TKE / \text{MJ}",
        yticklabelcolor=:green,
        ylabelcolor=:green,
        yaxisposition = :right,
        xaxisposition=:top,
        yticks=range(24, 40, 3),
        ax_series_kw...
    )
    
    ht_u_kw = (;
        colorrange=colorrange_u,
        lowclip=colormap_u[1],
        highclip=colormap_u[end],
        colormap=colormap_u,
        ht_u_kw...
    )
    ct_b_kw = (;
        color=(:black, 0.5),
        levels=b_levels,
        ct_b_kw...
    )
    ct_v_kw = (;
        colormap=colormap_u,
        levels=v_levels,
        ct_v_kw...
    )

    ht_Vq_kw = (;
        colorrange=colorrange_Vq,
        colormap=colormap_Vq,
        ht_Vq_kw...
    )
    ln_Vq_kw = (;
        color=ax_Vq_series.yticklabelcolor,
        ln_Vq_kw...
    )
    ln_TKE_kw = (;
        color=ax_TKE_series.yticklabelcolor,
        ln_TKE_kw...
    )

    # Mean and slice of u, v and b contours
    ht_u = heatmap!(ax_u, xsᶠ ./ 1000, zsᶜ, u; ht_u_kw...)
    ht_uh = heatmap!(ax_uh, xsᶠ ./ 1000, ysᶜ / 1000, uh; ht_u_kw...)
    contour!(ax_u, xsᶜ ./ 1000, zsᶜ, b; ct_b_kw...)
    contour!(ax_uh, xsᶜ ./ 1000, ysᶜ ./ 1000, bh; ct_b_kw...)
    #contour!(ax_u, xsᶜ ./ 1000, zsᶜ, v; ct_v_kw...)
    #contour!(ax_u, xsᶜ ./ 1000, zsᶜ, MLD; levels=[0], color=:blue, linestyle=:dash)
    #contour!(ax_uh, xsᶜ ./ 1000, ysᶜ ./ 1000, vh; ct_v_kw...)

    # PV and TKE lines and point
    ln_Vq = lines!(ax_Vq_series, full_times ./ 3600, Vq_series; ln_Vq_kw...)
    ln_TKE = lines!(ax_TKE_series, full_times ./ 3600, TKE_series; ln_TKE_kw...)
    scatter!(ax_Vq_series, Vq_point; ln_Vq_kw.color, markersize=14)
    scatter!(ax_TKE_series, TKE_point; ln_TKE_kw.color, markersize=14)
    scatter!(ax_Vq_series, Vq_point; color=:white, markersize=8)
    scatter!(ax_TKE_series, TKE_point; color=:white, markersize=8)

    # PV heatmap
    ht_Vq = heatmap!(ax_Vq, xsᶜ ./ 1000, zsᶜ, Vq; ht_Vq_kw...)

    Colorbar(fig[3, 1], ht_u; vertical=false, flipaxis=false, label=L"u / \text{cm}\,{s}^{-1}")
    Colorbar(fig[3, 2], ht_Vq; vertical=false, flipaxis=false, label=L"\overline{\chi}_{q < 0}")

    hidespines!(ax_TKE_series)
    hidexdecorations!(ax_TKE_series)

    hidexdecorations!(ax_uh; ticks=true)
    #hideydecorations!(ax_Vq; ticks=false)

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