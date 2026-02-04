function full_video(
        foldername,
        filename,
        frames,
        z;
        fig_kw=(; ), 
        ax_kw=(; ),
        axh_kw=(; ),
        ax_hist_kw=(; ),
        ht_u_kw=(; ),
        ht_q_kw=(; ),
        hist_kw=(; ),
        ct_b_kw=(; ),
        record_kw=(; ),
        σ=0,
        σh=0,
        σq=0,
        background=false
    )
    DFM = jldopen(joinpath(foldername, "DFM.jld2"))
    PV = jldopen(joinpath(foldername, "PV.jld2"))
    SLICES = jldopen(joinpath(foldername, "SLICES.jld2"))
    RORI = jldopen(joinpath(foldername, "RORI.jld2"))
    
    full_iterations, full_times = iterations_times(DFM)
    sp = simulation_parameters(DFM)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(DFM)
    inds = center_indices(DFM)
    colormap_u = to_colormap(:balance)
    colormap_q = to_colormap(:curl)
    z_indᶜ = zᶜbounds(DFM, z)
    z1, z2 = zᶜbounds(DFM, -0.95sp.H, -0.05sp.H)

    iterations = full_iterations[frames]
    times = full_times[frames]

    n = Observable(1)
    frame = @lift frames[$n]
    
    iteration = @lift iterations[$n]
    t = @lift times[$n]
    
    u_title = @lift let t = $t
        t_val = @sprintf "%03.1f" sp.f * t / 2π
        hr_val = @sprintf "%03.0f" t / 3600
        L"t=%$hr_val ~\text{hr}"
    end
    
    U = [-variable_strain_rate(t, sp) * x for t in times, x in xsᶠ, y in 1:1, z in 1:1] .* background
    
    fig = Figure(; 
        size=(960, 540),
        fig_kw...
    )
    
    colorrange_u = (-6, 6)
    colorrange_q = (-0.25, 0.25)

    u = @lift 100 * (get_field(a->filt(a, σ), DFM, "u_dfm", $iteration) .+ U[$n, :, 1, :])
    b = @lift get_field(a->filt(a, σ), DFM, "b_dfm", $iteration)
    
    uh = @lift 100 * (get_field(a->filt(a, σh), SLICES, "z050_u", $iteration) .+ U[$n, :, :, 1])
    bh = @lift get_field(a->filt(a, σh), SLICES, "z050_b", $iteration)

    q = @lift get_field(a->filt(a, σq), PV, "q", $iteration) ./ (sp.f * sp.N₀²)

    mask = [maskfromlines(x, z, regions.arrest) for x in xsᶜ, z in zsᶜ]
    Ri = @lift get_field(RORI, "N²", $iteration) ./ get_field(RORI, "S²", $iteration)
    tanhRi_bins = range(-1, 1, 100)
    tanhRis = (tanhRi_bins[1:end-1] .+ tanhRi_bins[2:end]) ./ 2
    ρ = @lift bin_counts(tanh.($Ri), tanhRi_bins)
    ρ_mask = @lift bin_counts(tanh.($Ri[mask]), tanhRi_bins)
    
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

    ax_q = Axis(fig[2, 2];
        limits=(-sp.Lh/2000, sp.Lh/2000, -sp.H, 0),
        xlabel=L"x / \text{km}",
        ylabel=L"z / \text{m}",
        xticks=[-2, -1, 0, 1, 2],
        ax_kw...
    )

    ax_hist = Axis(fig[1, 2];
        limits=(-1, 1, nothing, nothing),
        xlabel=L"\text{Ri} = N^2 / S^2",
        xticks = (tanh.([-Inf, -1, -0.5, 0, 0.5, 1.0, Inf]), [L"-\infty", "-1", "-0.5", "0", "0.5", "1", L"\infty"]),
        xminorticks = tanh.([-0.75, -0.25, 0.25, 0.75]),
        xminorticksvisible = true,
        ylabel=L"\rho",
        xaxisposition=:top,
        #yticks=range(0.24, 0.40, 3),
        ax_hist_kw...
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

    ht_q_kw = (;
        colorrange=colorrange_q,
        colormap=colormap_q,
        ht_q_kw...
    )

    # Mean and slice of u, v and b contours
    ht_u = heatmap!(ax_u, xsᶠ ./ 1000, zsᶜ, u; ht_u_kw...)
    ht_uh = heatmap!(ax_uh, xsᶠ ./ 1000, ysᶜ / 1000, uh; ht_u_kw...)
    contour!(ax_u, xsᶜ ./ 1000, zsᶜ, b; ct_b_kw...)
    contour!(ax_uh, xsᶜ ./ 1000, ysᶜ ./ 1000, bh; ct_b_kw...)
    #contour!(ax_u, xsᶜ ./ 1000, zsᶜ, v; ct_v_kw...)
    #contour!(ax_u, xsᶜ ./ 1000, zsᶜ, MLD; levels=[0], color=:blue, linestyle=:dash)
    #contour!(ax_uh, xsᶜ ./ 1000, ysᶜ ./ 1000, vh; ct_v_kw...)

    # Ri histogram
    lines!(ax_hist, tanhRis, ρ)
    lines!(ax_hist, tanhRis, ρ_mask; color=:magenta)

    # PV heatmap
    ht_q = heatmap!(ax_q, xsᶜ ./ 1000, zsᶜ, q; ht_q_kw...)

    Colorbar(fig[3, 1], ht_u; vertical=false, flipaxis=false, label=L"u / \text{cm}\,{s}^{-1}")
    Colorbar(fig[3, 2], ht_q; vertical=false, flipaxis=false, label=L"q / fN_0^2")

    #hidespines!(ax_TKE_series)
    #hidexdecorations!(ax_TKE_series)

    hidexdecorations!(ax_uh; ticks=true)
    #hideydecorations!(ax_q; ticks=false)

    rowgap!(fig.layout, 1, Relative(0.02))

    subfig_label!(fig[1, 1], 1)
    subfig_label!(fig[1, 2], 2)
    subfig_label!(fig[2, 1], 3)
    subfig_label!(fig[2, 2], 4)
    
    lines!(ax_u, [-sp.Lx / 2, sp.Lx / 2], [z, z]; color=(:red, 0.5), linestyle=:dash)
    lines!(ax_u, [-sp.Lx / 2, sp.Lx / 2], [z, z]; color=(:red, 0.5), linestyle=:dash)

    if length(frames) > 1
        record(fig, filename, 1:length(frames); record_kw...) do i
            n[] = i
            print("$i / $(length(frames)) \r")
        end
        println("")
    end
    close(DFM)
    close(SLICES)
    close(PV)
    close(RORI)

    fig
end