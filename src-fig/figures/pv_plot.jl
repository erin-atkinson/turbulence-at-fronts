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
        arrows = false,
        region = nothing
    )

    PV = joinpath(foldername, "PV.jld2")
    
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
    
    iterations = iterations[frames]
    times = times[frames]

    qs = timeseries_of(a->filt(a, σ), PV, "q", iterations) ./ (sp.f * sp.N₀²)
    
    u = get_field(a->filt(a, σ), DFM, "u_dfm", iterations[end]) .+ [-sp.α * x for x in xsᶠ, z in 1:1]
    u = (u[1:end-1, :] .+ u[2:end, :]) ./ 2
    w = get_field(a->filt(a, σ), DFM, "w_dfm", iterations[end])
    w = (w[:, 1:end-1] .+ w[:, 2:end]) ./ 2
    
    # Arrow plots
    arrow_i = 300:100:2048
    arrow_k = 150:30:256
    xs_u = xsᶜ[arrow_i] ./ 1000
    zs_u = zsᶜ[arrow_k]
    u = u[arrow_i, arrow_k] ./ 1000
    w = w[arrow_i, arrow_k]

    u = u .* qs[end, arrow_i, arrow_k]
    w = w .* qs[end, arrow_i, arrow_k]

    #u_norm = u ./ sqrt.(u.^2 .+ w.^2)
    #w_norm = w ./ sqrt.(u.^2 .+ w.^2)
    
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

    # Arrows
    #arrow_color = qs[end, arrow_i, arrow_k] .* 
    arr_kw = (; 
        lengthscale=5e5, 
        color = (:black, 0.3), 
#        colormap = [RGBA(0, 0, 0, 0), RGBA(0, 0, 0, 1)],
#        colorrange = (0, maximum(abs, w) ./ 2),
        arr_kw...
    )

    axes = map(enumerate(titles)) do (i, title)
        Axis(fig[1, i]; ax_kw..., title)
    end

    hts = map(enumerate(axes)) do (i, ax)
        ht = heatmap!(ax, xsᶜ/1000, zsᶜ, qs[i, :, :]; ht_kw...)
        subfig_label!(fig[1, i], i)
        ht
    end

    hideydecorations!.(axes[2:end]; ticks=false)
    
    Colorbar(fig[1, length(axes)+1], hts[1]; label=L"q / fN_0^2")

    if region !== nothing
        mask = [maskfromlines(1000x, z, region) for x in range(-sp.Lh/2000, sp.Lh/2000, 1000), z in range(-sp.Lz, 0, 1000)]
        contour!(axes[end], range(-sp.Lh/2000, sp.Lh/2000, 1000), range(-sp.Lz, 0, 1000), mask, levels=[0.5]; color=:magenta, linestyle=:dashdot, linewidth=1)
    end
    if arrows
        arrows2d!(axes[end], xs_u, zs_u, u, w; arr_kw...)
    end
    fig
end
