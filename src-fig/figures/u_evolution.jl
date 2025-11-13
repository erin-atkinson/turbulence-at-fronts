function u_evolution(
        foldername,
        frames,
        z,
        region=nothing; 
        fig_kw=(; ), 
        ax_kw=(; ),
        axh_kw=(; ),
        ht_kw=(; ),
        ct_kw=(; ),
        σ=0,
        σh=0,
        background=true,
        mixed_depth=false,
    )

    OUTPUT = joinpath(foldername, "OUTPUT.jld2")
    DFM = joinpath(foldername, "DFM.jld2")

    iterations, times = iterations_times(DFM)
    sp = simulation_parameters(DFM)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(DFM)
    inds = centre_indices(DFM)
    colormap = to_colormap(:balance)
    z_indᶜ = zᶜbounds(DFM, z)
    
    fig = Figure(; 
        size=(1000, 400),
        fig_kw...
    )

    iterations = iterations[frames]
    times = times[frames]

    U = [-variable_strain_rate(t, sp) * x for t in times, x in xsᶠ, yz in 1:1] .* background

    u_dfm = 100 * (timeseries_of(a->filt(a, σ), DFM, "u_dfm", iterations) .+ U)
    uh = 100 * (timeseries_of(a->filt(a[:, :, z_indᶜ], σh), OUTPUT, "u", iterations) .+ U)
    
    b_dfm = timeseries_of(a->filt(a, σ), DFM, "b_dfm", iterations)
    bh = timeseries_of(a->filt(a[:, :, z_indᶜ], 1, 1), OUTPUT, "b", iterations)

    MLD = b_dfm .- (b_dfm[:, :, end:end] .- (b_levels[3] - b_levels[1])) 
    
    u_max = max(maximum(abs, u_dfm[:, inds, :]), maximum(abs, uh[:, inds, :]))
    
    titles = map(times) do t
        t_val = prettytime(sp.f * t / 2π; l=1, digits=1)
        hr_val = prettytime(t / 3600; l=3, digits=0)
        
        L"ft / 2\pi = %$(t_val) \quad t = %$(hr_val)~\text{hr}"
    end

    ax_kw = (;
        xlabel=L"x/\text{km}",
        ylabel=L"z/\text{m}",
        xticks=[-2, -1, 0, 1, 2],
        limits=(-sp.ℓ / 1000, sp.ℓ / 1000, -sp.H, 0),
        ax_kw...
    )
    
    axh_kw = (;
        xlabel=L"x/\text{km}",
        ylabel=L"y/\text{km}",
        xticks=[-2, -1, 0, 1, 2],
        limits=(-sp.ℓ / 1000, sp.ℓ / 1000, -sp.Ly / 2000, sp.Ly / 2000),
        axh_kw...
    )

    axs = map(1:length(frames), frames, titles) do i, frame, title
        Axis(fig[i, 1]; ax_kw..., title)
    end

    axsh = map(1:length(frames), frames) do i, frame
        Axis(fig[i, 2]; axh_kw...)
    end
    
    hidexdecorations!.(axs[1:end-1]; ticks=false)
    hidexdecorations!.(axsh[1:end-1]; ticks=false)
    hideydecorations!.(axsh; ticks=false)

    ht_kw = (;
        colormap,
        colorrange=(-u_max, u_max),
        ht_kw...
    )
    
    hts = map(1:length(frames), axs) do i, ax
        heatmap!(ax, xsᶠ / 1000, zsᶜ, u_dfm[i, :, :]; ht_kw...)
    end

    htsh = map(1:length(frames), axsh) do i, ax
        heatmap!(ax, xsᶠ / 1000, ysᶜ / 1000, uh[i, :, :]; ht_kw...)
    end

    ct_kw = (;
        color=(:black, 0.5),
        levels=b_levels,
        ct_kw...
    )
    
    cts = map(1:length(frames), axs) do i, ax
        contour!(ax, xsᶜ / 1000, zsᶜ, b_dfm[i, :, :]; ct_kw...)
    end
    
    map(1:length(frames), axs) do i, ax
        mixed_depth && contour!(ax, xsᶜ / 1000, zsᶜ, MLD[i, :, :]; levels=[0], color=:blue, linestyle=:dash)
    end

    ctsh = map(1:length(frames), axsh) do i, ax
        contour!(ax, xsᶜ / 1000, ysᶠ / 1000, bh[i, :, :]; ct_kw...)
    end

    Colorbar(fig[length(frames) + 1, 1:2], hts[1]; label=L"u / \text{cm}\,{s}^{-1}", vertical=false, flipaxis=false)

    for i in 1:length(frames)
        subfig_label!(fig[i, 1], 2(i-1) +1)
        subfig_label!(fig[i, 2], 2i)
    end
    for ax in axs
        lines!(ax, [ax_kw.limits[1], ax_kw.limits[2]], [z, z]; color=(:red, 0.5), linestyle=:dash)
    end
    if region !== nothing
        mask = [maskfromlines(1000x, z, region) for x in range(-sp.Lh/2000, sp.Lh/2000, 1000), z in range(-sp.Lz, 0, 1000)]
        contour!(axs[end], range(-sp.Lh/2000, sp.Lh/2000, 1000), range(-sp.Lz, 0, 1000), mask, levels=[0.5]; color=:magenta, linestyle=:dashdot, linewidth=1)
    end
    fig
end