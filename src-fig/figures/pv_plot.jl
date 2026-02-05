# pv_video.jl

# Potential vorticity and proportion of negative pv

function pv_plot(
        foldername,
        frames; 
        fig_kw = (; ), 
        ax_kw = (; ),
        ht_kw = (; ),
        arr_kw = (; ),
        σ = 0,
        arrow = false,
        region = nothing
    )

    PV = joinpath(foldername, "PV.jld2")
    RORI = joinpath(foldername, "RORI.jld2")
    
    iterations, times = iterations_times(PV)
    sp = simulation_parameters(PV)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(PV)
    inds = center_indices(PV)
    
    colormap = to_colormap(:curl)
    
    fig_kw = (;
        size=(800, 300),
        fig_kw...
    )
    
    fig = Figure(; fig_kw...)
    
    u = time_average_of(a->filt(a, σ), DFM, "u_dfm", iterations[end-100:end]) .+ [-sp.α * x for x in xsᶠ, z in 1:1]
    u = (u[1:end-1, :] .+ u[2:end, :]) ./ 2
    w = time_average_of(a->filt(a, σ), DFM, "w_dfm", iterations[end-100:end])
    w = (w[:, 1:end-1] .+ w[:, 2:end]) ./ 2

    # Arrow plots
    xs_u = xsᶜ[1:60:end] ./ 1000
    zs_u = zsᶜ[1:20:end]
    u = u[1:60:end, 1:20:end] ./ 1000
    w = w[1:60:end, 1:20:end]
    
    iterations = iterations[frames]
    times = times[frames]

    qs = timeseries_of(a->filt(a, σ), PV, "q", iterations) ./ (sp.f * sp.N₀²)
    Ris = timeseries_of(a->filt(a, σ), RORI, "N²", iterations) ./ (timeseries_of(a->filt(a, σ), RORI, "S²", iterations))
    #for i in 1:length(frames)
    #    Ris[i, :, :] .= filt(Ris[i, :, :], σ)
    #end
    # Set each title
    titles = map(times) do t
        t_val = @sprintf "%03.1f" sp.f * t / 2π
        hr_val = @sprintf "%03.0f" t / 3600
        
        L"t = %$(hr_val)~\text{hr}"
    end
    
    ax_kw = (;
        xlabel = L"x/\text{km}",
        ylabel = L"z/\text{m}",
        limits = (-sp.Lh / 2000, sp.Lh / 2000, -sp.H, 0),
        ax_kw...
    )
    
    # Heatmaps
    ht_kw = (;
        colorrange = (-0.5, 0.5),
        colormap,
        highclip = colormap[end],
        lowclip = colormap[1],
        ht_kw...
    )

    axes = map(enumerate(titles)) do (i, title)
        Axis(fig[1, i]; ax_kw..., title)
    end

    axes_Ri = map(enumerate(titles)) do (i, title)
        Axis(fig[2, i]; ax_kw...)
    end

    hts = map(enumerate(axes)) do (i, ax)
        ht = heatmap!(ax, xsᶜ/1000, zsᶜ, qs[i, :, :]; ht_kw...)
        subfig_label!(fig[1, i], i)
        ht
    end

    hts_Ri = map(enumerate(axes_Ri)) do (i, ax)
        cmap = to_colormap(:amp)
        #ht = heatmap!(ax, xsᶜ/1000, zsᶜ, tanh.(Ris[i, :, :]); 
        #    colorrange=(0, 1.0),
        #    colormap = cmap,
        #    lowclip = cmap[1]
        #)

        ht = contourf!(ax, xsᶜ/1000, zsᶜ, 0.999 .* tanh.(Ris[i, :, :]); 
            levels=range(0, 1, 11),
            colormap = cmap,
            extendlow = cmap[1]
        )
        subfig_label!(fig[2, i], i+3)
        ht
    end

    hideydecorations!.(axes[2:end]; ticks=false)
    hidexdecorations!.(axes; ticks=false)
    hideydecorations!.(axes_Ri[2:end]; ticks=false)
    
    Colorbar(fig[1, length(axes)+1], hts[1]; label=L"q / fN_0^2")
    
    ticks = (tanh.([-1, -0.5, 0, 0.5, 1.0, Inf]), ["-1", "-0.5", "0", "0.5", "1", L"\infty"])
    minorticks = tanh.([-0.75, -0.25, 0.25, 0.75])
    minorticksvisible = true
    Colorbar(fig[2, length(axes)+1], hts_Ri[1]; label=L"\text{Ri}", ticks, minorticks, minorticksvisible)

    arr_kw = (; 
        lengthscale=5000, 
        color = abs.(w)[:], 
        colormap = [RGBA(0, 0, 0, 0), RGBA(0, 0, 0, 1)],
        #colorrange = (0, maximum(abs, w) / 2),
        colorrange = (0, 100 / (24 * 3600)),
        arr_kw...
    )
    @info maximum(abs, w) ./ 2
    if arrow
        arrows2d!(axes[end], xs_u, zs_u, u, w; arr_kw...)
    end
    
    if region !== nothing
        mask = [maskfromlines(1000x, z, region) for x in range(-sp.Lh/2000, sp.Lh/2000, 1000), z in range(-sp.Lz, 0, 1000)]
        contour!(axes[end], range(-sp.Lh/2000, sp.Lh/2000, 1000), range(-sp.Lz, 0, 1000), mask, levels=[0.5]; color=:magenta, linestyle=:dashdot, linewidth=1)
    end
    fig
end
