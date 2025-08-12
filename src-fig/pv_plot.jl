# pv_video.jl

# Potential vorticity and proportion of negative pv

function instability_plot(
        foldername,
        frames,
        region=nothing; 
        fig_kw=(; ), 
        ax_kw=(; ),
        ht_kw=(; ),
        σ=0,
    )
    
    iterations, times = iterations_times(foldername)
    sp = simulation_parameters(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    inds = centre_indices(foldername)
    colormap = to_colormap(:diff)
    #colormap = [repeat(colormap[1:1], 2*length(colormap)÷4); colormap]
    fig_kw = (;
        size=(800, 300),
        fig_kw...
    )
    
    fig = Figure(; fig_kw...)
    
    iterations = iterations[frames]
    times = times[frames]

    Vq = timeseries_of(a->filt(a, σ), joinpath(foldername, "PV.jld2"), "Vq_dfm", iterations)
    
    # Set each title

    titles = map(times) do t
        t_val = prettytime(sp.f * t / 2π; l=1, digits=1)
        hr_val = prettytime(t / 3600; l=3, digits=0)
        L"ft / 2\pi = %$(t_val) \quad t = %$(hr_val)~\text{hr}"
    end
    
    ax_kw = (;
        xlabel=L"x/\text{km}",
        ylabel=L"z/\text{m}",
        limits=(-sp.ℓ / 1000, sp.ℓ / 1000, -sp.H, 0),
        ax_kw...
    )
    
    # Heatmaps
    ht_kw=(;
        colorrange=(0, 1),
        colormap,
        ht_kw...
    )

    axes = map(enumerate(titles)) do (i, title)
        Axis(fig[1, i]; ax_kw..., title)
    end

    hts = map(enumerate(axes)) do (i, ax)
        ht = heatmap!(ax, xsᶜ/1000, zsᶜ, Vq[i, :, :]; ht_kw...)
        subfig_label!(fig[1, i], i)
        ht
    end
    hideydecorations!.(axes[2:end])
    Colorbar(fig[1, length(axes)+1], hts[1]; label=L"\overline{\chi}_{q < 0}")

    if region !== nothing
        mask = [maskfromlines(1000x, z, region) for x in range(-sp.Lh/2000, sp.Lh/2000, 1000), z in range(-sp.Lz, 0, 1000)]
        contour!(axes[end], range(-sp.Lh/2000, sp.Lh/2000, 1000), range(-sp.Lz, 0, 1000), mask, levels=[0.5]; color=:magenta, linestyle=:dashdot, linewidth=1)
    end
    fig
end
