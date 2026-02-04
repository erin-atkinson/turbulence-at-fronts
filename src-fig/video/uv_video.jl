function uv_video(
        foldername,
        filename,
        frames;
        fig_kw=(; ), 
        ax_kw=(; ),
        axh_kw=(; ),
        ht_kw=(; ),
        ct_kw=(; ),
        record_kw=(; ),
        σ=0,
        σh=0,
        background=false
    )
    
    SLICES = joinpath(foldername, "SLICES.jld2")
    DFM = joinpath(foldername, "DFM.jld2")

    iterations, times = iterations_times(DFM)
    sp = simulation_parameters(DFM)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(DFM)
    inds = center_indices(DFM)
    
    colormap = to_colormap(:balance)

    iterations = iterations[frames]
    times = times[frames]

    n = Observable(1)
    iteration = @lift iterations[$n]
    t = @lift times[$n]
    u_title = @lift let t = $t
        t_val = @sprintf "%03.1f" sp.f * t / 2π
        hr_val = @sprintf "%03.0f" t / 3600
        
        L"ft / 2\pi = %$(t_val) \quad t = %$(hr_val)~\text{hr}"
    end
    v_title = @lift let t = $t
        t_val = @sprintf "%03.1f" sp.f * t / 2π
        hr_val = @sprintf "%03.0f" t / 3600
        
        L"ft / 2\pi = %$(t_val) \quad t = %$(hr_val)~\text{hr}"
    end
    
    U = [-variable_strain_rate(t, sp) * x for t in times, x in xsᶠ, y in 1:1, z in 1:1] .* background
    
    fig = Figure(; 
        size=(960, 540),
        fig_kw...
    )
    
    colorrange_u = (-6, 6)
    colorrange_v = (-20, 20)

    SLICES = jldopen(joinpath(foldername, "SLICES.jld2"))
    DFM = jldopen(joinpath(foldername, "DFM.jld2"))

    b0 = get_field(a->mean(a), DFM, "b_dfm", iterations[1])
    u = @lift 100 * (get_field(a->filt(a, σ), DFM, "u_dfm", $iteration) .+ U[$n, :, 1, :])
    v = @lift 100 * get_field(a->filt(a, σ), DFM, "v_dfm", $iteration)
    b = @lift get_field(a->filt(a .- (mean(a) - b0), σ), DFM, "b_dfm", $iteration)
    
    uh = @lift 100 * (get_field(a -> filt(a, σh), SLICES, "z050_u", $iteration) .+ U[$n, :, :, 1])
    vh = @lift 100 * get_field(a -> filt(a, σh), SLICES, "z050_v", $iteration)
    bh = @lift get_field(a -> filt(a .- (mean(a) - b0), σh), SLICES, "z050_b", $iteration)
    
    ax_u = Axis(fig[2, 1];
        limits=(-sp.Lh/2000, sp.Lh/2000, -sp.H, 0),
        xlabel=L"x / \text{km}",
        ylabel=L"z / \text{m}",
        xticks=[-2, -1, 0, 1, 2],
        ax_kw...
    )
    ax_v = Axis(fig[2, 2];
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
    ax_vh = Axis(fig[1, 2];
        limits=(-sp.Lh/2000, sp.Lh/2000, -sp.Ly/2000, sp.Ly/2000),
        xlabel=L"x / \text{km}",
        ylabel=L"y / \text{km}",
        title=v_title,
        axh_kw...
    )

    ht_u_kw = (;
        colorrange=colorrange_u,
        lowclip=colormap[1],
        highclip=colormap[end],
        colormap,
        ht_kw...
    )
    ht_v_kw = (;
        colorrange=colorrange_v,
        lowclip=colormap[1],
        highclip=colormap[end],
        colormap,
        ht_kw...
    )
    ct_kw = (;
        color=(:black, 0.5),
        levels=b_levels,
        ct_kw...
    )

    ht_u = heatmap!(ax_u, xsᶠ ./ 1000, zsᶜ, u; ht_u_kw...)
    ht_v = heatmap!(ax_v, xsᶜ ./ 1000, zsᶜ, v; ht_v_kw...)

    contour!(ax_u, xsᶜ ./ 1000, zsᶜ, b; ct_kw...)
    contour!(ax_v, xsᶜ ./ 1000, zsᶜ, b; ct_kw...)
    
    #contour!(ax_u, xsᶜ ./ 1000, zsᶜ, MLD; levels=[0], color=:blue, linestyle=:dash)
    #contour!(ax_v, xsᶜ ./ 1000, zsᶜ, MLD; levels=[0], color=:blue, linestyle=:dash)
    
    ht_uh = heatmap!(ax_uh, xsᶠ ./ 1000, ysᶜ / 1000, uh; ht_u_kw...)
    ht_vh = heatmap!(ax_vh, xsᶜ ./ 1000, ysᶠ / 1000, vh; ht_v_kw...)

    contour!(ax_uh, xsᶜ ./ 1000, ysᶜ ./ 1000, bh; ct_kw...)
    contour!(ax_vh, xsᶜ ./ 1000, ysᶜ ./ 1000, bh; ct_kw...) 
    
    Colorbar(fig[3, 1], ht_u; vertical=false, flipaxis=false, label=L"u / \text{cm}\,\text{s}^{-1}")
    Colorbar(fig[3, 2], ht_v; vertical=false, flipaxis=false, label=L"v / \text{cm}\,\text{s}^{-1}")

    hidexdecorations!(ax_uh; ticks=true)
    hidexdecorations!(ax_vh; ticks=true)

    hideydecorations!(ax_v; ticks=false)
    hideydecorations!(ax_vh; ticks=false)

    rowgap!(fig.layout, 1, Relative(0.02))

    subfig_label!(fig[1, 1], 1)
    subfig_label!(fig[1, 2], 2)
    subfig_label!(fig[2, 1], 3)
    subfig_label!(fig[2, 2], 4)
    
    for ax in [ax_u, ax_v]
        lines!(ax, [-sp.Lx / 2000, sp.Lx / 2000], [-0.5sp.H, -0.5sp.H]; color=(:red, 0.5), linestyle=:dash)
    end
    if length(frames) > 1
        record(fig, filename, 1:length(frames); record_kw...) do i
            n[] = i
            print("$i / $(length(frames)) \r")
        end
        println("")
    end
    close(DFM)
    close(SLICES)

    fig
end