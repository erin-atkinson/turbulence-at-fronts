function u_evolution(
        foldername,
        frames,
        x_slice = 0,
        region = nothing; 
        fig_kw = (; ), 
        ax_kw = (; ),
        axh_kw = (; ),
        ax_slice_kw = (; ),
        ht_kw = (; ),
        ct_kw = (; ),
        arr_kw = (; ),
        σ = 0,
        σh = 0,
        background = true,
        mixed_depth = false,
        arrow = false
    )

    SLICES = joinpath(foldername, "SLICES.jld2")
    DFM = joinpath(foldername, "DFM.jld2")
    OUTPUT = joinpath(foldername, "output.jld2")

    full_iterations, full_times = iterations_times(DFM)
    sp = simulation_parameters(DFM)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(DFM)
    inds = center_indices(DFM)
    iᶠ = xᶠbounds(DFM, x_slice)
    iᶜ = xᶜbounds(DFM, x_slice)

    colormap = to_colormap(:balance)
    
    fig = Figure(; 
        size=(1000, 400),
        fig_kw...
    )

    iterations = full_iterations[frames]
    times = full_times[frames]

    U = if background
        [-variable_strain_rate(t, sp) * x for t in times, x in xsᶠ, yz in 1:1]
    else
        0
    end

    u_dfm = 100 * (timeseries_of(a->filt(a, σ), DFM, "u_dfm", iterations) .+ U)
    uh = 100 * (timeseries_of(a->filt(a, σh), SLICES, "z050_u", iterations) .+ U)

    u_slice = 100 * (timeseries_of(a->filt(a[iᶠ, :, :], σ), OUTPUT, "u", iterations) .+ U[:, iᶠ:iᶠ, 1:1])

    b_dfm = timeseries_of(a->filt(a, σ), DFM, "b_dfm", iterations)
    bh = timeseries_of(a->filt(a, 1, 1), SLICES, "z050_b", iterations)
    
    b_slice = timeseries_of(a->filt(a[iᶜ, :, :], σ), OUTPUT, "b", iterations)

    u = time_average_of(a->filt(a, σ), DFM, "u_dfm", full_iterations[end-100:end]) .+ [-sp.α * x for x in xsᶠ, z in 1:1]
    u = (u[1:end-1, :] .+ u[2:end, :]) ./ 2
    w = time_average_of(a->filt(a, σ), DFM, "w_dfm", full_iterations[end-100:end])
    w = (w[:, 1:end-1] .+ w[:, 2:end]) ./ 2

    # Arrow plots
    xs_u = xsᶜ[1:60:end] ./ 1000
    zs_u = zsᶜ[1:20:end]
    u = u[1:60:end, 1:20:end] ./ 1000
    w = w[1:60:end, 1:20:end]
    
    u_max = max(maximum(abs, u_dfm[:, inds, :]), maximum(abs, uh[:, inds, :]))
    
    titles = map(times) do t
        t_val = @sprintf "%.1f" sp.f * t / 2π
        hr_val = @sprintf "%.0f" t / 3600
        
        L"ft / 2\pi = %$(t_val) \quad t = %$(hr_val)~\text{hr}"
    end

    ax_kw = (;
        xlabel=L"x/\text{km}",
        ylabel=L"z/\text{m}",
        xticks=[-2, -1, 0, 1, 2],
        limits=(-sp.Lh/ 2000, sp.Lh / 2000, -sp.H, 0),
        ax_kw...
    )
    
    axh_kw = (;
        xlabel=L"x/\text{km}",
        ylabel=L"y/\text{km}",
        xticks=[-2, -1, 0, 1, 2],
        limits=(-sp.Lh / 2000, sp.Lh / 2000, -sp.Ly / 2000, sp.Ly / 2000),
        axh_kw...
    )

    ax_slice_kw = (;
        xlabel=L"y/\text{km}",
        ylabel=L"z/\text{m}",
        xticks=[-0.1, 0, 0.1],
        limits=(-sp.Ly / 2000, sp.Ly / 2000, -sp.H, 0),
        ax_slice_kw...
    )

    axs = map(1:length(frames), frames, titles) do i, frame, title
        Axis(fig[i, 1]; ax_kw..., title)
    end

    axsh = map(1:length(frames), frames) do i, frame
        Axis(fig[i, 3]; axh_kw...)
    end

    axs_slice = map(1:length(frames), frames) do i, frame
        Axis(fig[i, 2]; ax_slice_kw...)
    end
    
    hidexdecorations!.(axs[1:end-1]; ticks=false)
    hidexdecorations!.(axsh[1:end-1]; ticks=false)
    hidexdecorations!.(axs_slice[1:end-1]; ticks=false)
    hideydecorations!.(axs_slice; ticks=false)

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

    htsh = map(1:length(frames), axs_slice) do i, ax
        heatmap!(ax, ysᶜ / 1000, zsᶜ, u_slice[i, :, :]; ht_kw...)
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
        #mixed_depth && contour!(ax, xsᶜ / 1000, zsᶜ, MLD[i, :, :]; levels=[0], color=:blue, linestyle=:dash)
    end

    ctsh = map(1:length(frames), axsh) do i, ax
        contour!(ax, xsᶜ / 1000, ysᶜ / 1000, bh[i, :, :]; ct_kw...)
    end

    ctsh = map(1:length(frames), axs_slice) do i, ax
        contour!(ax, ysᶜ / 1000, zsᶜ, b_slice[i, :, :]; ct_kw...)
    end

    Colorbar(fig[length(frames) + 1, 1:3], hts[1]; label=L"(u + U) / \text{cm}\,{s}^{-1}", vertical=false, flipaxis=false)

    for i in 1:length(frames)
        subfig_label!(fig[i, 1], 3(i-1) + 1)
        subfig_label!(fig[i, 2], 3(i-1) + 2)
        subfig_label!(fig[i, 3], 3i)
    end
    for ax in axs
        lines!(ax, [ax_kw.limits[1], ax_kw.limits[2]], [-50, -50]; color=(:red, 1.0), linestyle=:dash, linewidth=1)
        lines!(ax, [x_slice, x_slice] ./ 1000, [ax_kw.limits[3], ax_kw.limits[4]]; color=(:red, 1.0), linestyle=:dash, linewidth=1)
    end
    if region !== nothing
        mask = [maskfromlines(1000x, z, region) for x in range(-sp.Lh/2000, sp.Lh/2000, 1000), z in range(-sp.Lz, 0, 1000)]
        contour!(axs[end], range(-sp.Lh/2000, sp.Lh/2000, 1000), range(-sp.Lz, 0, 1000), mask, levels=[0.5]; color=:magenta, linestyle=:dashdot, linewidth=1)
    end

    arr_kw = (; 
        lengthscale=5000, 
        color = abs.(w)[:], 
        colormap = [RGBA(0, 0, 0, 0), RGBA(0, 0, 0, 1)],
        #colorrange = (0, maximum(abs, w) / 2),
        colorrange = (0, 100 / (24 * 3600)),
        arr_kw...
    )
    #@info maximum(abs, w) ./ 2
    if arrow
        arrows2d!(axs[end], xs_u, zs_u, u, w; arr_kw...)
    end
    
    fig
end
