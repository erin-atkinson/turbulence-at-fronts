function uv_video(
        foldername,
        filename,
        frames,
        z;
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
    iterations, times = iterations_times(foldername)
    sp = simulation_parameters(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    inds = centre_indices(foldername)
    colormap = to_colormap(:balance)
    z_indᶜ = zᶜbounds(foldername, z)

    iterations = iterations[frames]
    times = times[frames]

    n = Observable(1)
    iteration = @lift iterations[$n]
    t = @lift times[$n]
    u_title = @lift let time_string = prettytime($t / 3600; digits=1)
        L"u, t=%$time_string ~\text{hr}"
    end
    v_title = @lift let time_string = prettytime($t / 3600; digits=1)
        L"v, t=%$time_string ~\text{hr}"
    end
    
    U = [-sp.α * x for x in xsᶠ, y in 1:1] .* background
    V = [sp.α * y for x in 1:1, y in ysᶠ] .* background
    
    fig = Figure(; 
        size=(960, 540),
        fig_kw...
    )
    
    DFM = jldopen(joinpath(foldername, "DFM.jld2"))
    OUTPUT = jldopen(joinpath(foldername, "output.jld2"))
    
    colorrange_u = (-0.1, 0.1)
    colorrange_v = (-0.3, 0.3)

    u = @lift get_field(DFM, "u_dfm", $iteration) .+ U
    v = @lift get_field(DFM, "v_dfm", $iteration)
    
    uh = @lift get_field(a->a[:, :, z_indᶜ], OUTPUT, "u", $iteration) .+ U
    vh = @lift get_field(a->a[:, :, z_indᶜ], OUTPUT, "v", $iteration) .+ V
    
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

    ht_u = heatmap!(ax_u, xsᶠ ./ 1000, zsᶜ, u; ht_u_kw...)
    ht_v = heatmap!(ax_v, xsᶜ ./ 1000, zsᶜ, v; ht_v_kw...)
    
    ht_uh = heatmap!(ax_uh, xsᶠ ./ 1000, ysᶜ / 1000, uh; ht_u_kw...)
    ht_vh = heatmap!(ax_vh, xsᶜ ./ 1000, ysᶠ / 1000, vh; ht_v_kw...)
    
    Colorbar(fig[3, 1], ht_u; vertical=false, flipaxis=false, label=L"u / \text{ms}^{-1}")
    Colorbar(fig[3, 2], ht_v; vertical=false, flipaxis=false, label=L"v / \text{ms}^{-1}")

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
        lines!(ax, [-sp.Lx / 2, sp.Lx / 2], [z, z]; color=(:red, 0.5), linestyle=:dash)
    end
    if length(frames) > 1
        record(fig, filename, 1:length(frames); record_kw...) do i
            n[] = i
            print("$i / $(length(frames)) \r")
        end
        println("")
    end
    close(DFM)
    close(OUTPUT)

    fig
end