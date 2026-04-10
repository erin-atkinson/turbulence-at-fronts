# pv_video.jl

# Potential vorticity and Richardson number

function pv_plot(
        foldername,
        f1, f2, f3; 
        fig_kw = (; ), 
        ax_kw = (; ),
        ht_kw = (; ),
        arr_kw = (; ),
        σ = 0,
        arrow = false,
        region = nothing
    )

    frames = [f1, f2, f3]
    
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
    zs_u = zsᶜ[1:30:end]
    u = u[1:60:end, 1:30:end] ./ 1000
    w = w[1:60:end, 1:30:end]

    # Suppress vertical velocities larger than 1mm/s
    s = abs.(w)
    s_max = 0.001
    s = map(s) do S
        S > s_max ? S/s_max : 1
    end
    
    u = u ./ s
    w = w ./ s
    
    iterations = iterations[frames]
    times = times[frames]

    qs = timeseries_of(a->filt(a, σ), PV, "q", iterations) ./ (sp.f * sp.N₀²)
    Ri = get_field(a->filt(a, σ), RORI, "N²", iterations[end]) ./ get_field(a->filt(a, σ), RORI, "S²", iterations[end])

    # Set each title
    hr_vals = map(times) do t
        @sprintf "%.0f" t / 3600
    end
    titles = [
        L"t = %$(hr_vals[1])~\text{hr}",
        L"t = %$(hr_vals[2])~\text{hr}",
        L"q / fN_0^2, \quad t = %$(hr_vals[3])~\text{hr}",
        L"\text{Ri}, \quad t = %$(hr_vals[3])~\text{hr}",
    ]
    
    ax_kw = (;
        xlabel = L"x/\text{km}",
        ylabel = L"z/\text{m}",
        limits = (-sp.Lh / 2000, sp.Lh / 2000, -sp.H, 0),
        ax_kw...
    )
    
    ax_f1 = Axis(fig[1, 1]; ax_kw..., title=titles[1])
    hidexdecorations!(ax_f1; ticks=false)
    
    ax_f2 = Axis(fig[1, 2]; ax_kw..., title=titles[2])
    hidexdecorations!(ax_f2; ticks=false)
    hideydecorations!(ax_f2; ticks=false)
    
    ax_f3_q = Axis(fig[2, 1]; ax_kw..., title=titles[3])
    
    ax_f3_Ri = Axis(fig[2, 2]; ax_kw..., title=titles[4])
    hideydecorations!(ax_f3_Ri; ticks=false)
    
    # Heatmaps
    ht_kw = (;
        colorrange = (-0.5, 0.5),
        colormap,
        highclip = colormap[end],
        lowclip = colormap[1],
        ht_kw...
    )

    ht = heatmap!(ax_f1, xsᶜ/1000, zsᶜ, qs[1, :, :]; ht_kw...)
    subfig_label!(fig[1, 1], 1)

    ht = heatmap!(ax_f2, xsᶜ/1000, zsᶜ, qs[2, :, :]; ht_kw...)
    subfig_label!(fig[1, 2], 2)

    ht = heatmap!(ax_f3_q, xsᶜ/1000, zsᶜ, qs[3, :, :]; ht_kw...)
    subfig_label!(fig[2, 1], 3)

    cmap = to_colormap(:amp)
    ctf = contourf!(ax_f3_Ri, xsᶜ/1000, zsᶜ, 0.999 .* tanh.(Ri); 
            levels=range(0, 1, 11),
            colormap = cmap,
            extendlow = cmap[1]
        )
    subfig_label!(fig[2, 2], 4)
    
    
    Colorbar(fig[1, 3], ht; label=L"q / fN_0^2")
    
    ticks = (tanh.([-1, -0.5, 0, 0.5, 1.0, Inf]), ["-1", "-0.5", "0", "0.5", "1", L"\infty"])
    minorticks = tanh.([-0.75, -0.25, 0.25, 0.75])
    minorticksvisible = true
    Colorbar(fig[2, 3], ctf; label=L"\text{Ri}", ticks, minorticks, minorticksvisible)

    arr_kw = (; 
        lengthscale=10000, 
        color = abs.(w)[:], 
        colormap = [RGBA(0, 0, 0, 0.0), RGBA(0, 0, 0, 1), RGBA(0, 0, 0, 1)],
        #colorrange = (0, maximum(abs, w) / 2),
        colorrange = (0, 80 / (24 * 3600)),
        arr_kw...
    )

    if arrow
        arrows2d!(ax_f3_q, xs_u, zs_u, u, w; arr_kw...)
    end
    
    fig
end
