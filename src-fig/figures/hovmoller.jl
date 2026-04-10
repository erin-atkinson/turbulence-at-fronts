function hovmoller(
        foldername,
        z;
        fig_kw=(; ), 
        ax_kw=(; ),
        ht_kw=(; ),
        ct_kw=(; ),
        σ=0,
        background=true,
    )

    DFM = joinpath(foldername, "DFM.jld2")
    
    iterations, times = iterations_times(DFM)
    sp = simulation_parameters(DFM)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(DFM)
    inds = center_indices(DFM)
    colormap = to_colormap(:balance)
    
    zᶜ = zᶜbounds(DFM, z)
    zᶠ = zᶠbounds(DFM, z)
    zᶜ_surface = zᶜbounds(DFM, -5)
    zᶠ_surface = zᶠbounds(DFM, -5)
    
    fig = Figure(; 
        size=(800, 400),
        fig_kw...
    )

    U = [-variable_strain_rate(t, sp) * x for t in times, x in xsᶠ] .* background
    
    u = 100 * (timeseries_of(a->filt(a[:, zᶠ], σ), DFM, "u_dfm", iterations) .+ U)
    b = timeseries_of(a->a[:, zᶜ], DFM, "b_dfm", iterations)
    
    u_surface = 100 * (timeseries_of(a->filt(a[:, zᶠ_surface], σ), DFM, "u_dfm", iterations) .+ U)
    b_surface = timeseries_of(a->a[:, zᶜ_surface], DFM, "b_dfm", iterations)

    b_offset = filt(b .+ mean(b[1:1, inds] .- b[:, inds]; dims=2), σ)
    b_surface_offset = filt(b_surface .+ mean(b_surface[1:1, inds] .- b_surface[:, inds]; dims=2), σ)

    b = filt(b, σ)
    b_surface = filt(b_surface, σ)
    
    u_max = max(maximum(abs, u[:, inds]), maximum(abs, u_surface[:, inds]))

    ax_kw = (;
        ylabel=L"t/\text{hr}",
        xlabel=L"x/\text{km}",
        xticks=[-2, -1, 0, 1, 2],
        limits=(-sp.Lh / 2000, sp.Lh / 2000, 0, 200),
        ax_kw...
    )
    ax = Axis(fig[1, 1]; ax_kw..., title=L"z=%$z\,\text{m}")
    ax_surface = Axis(fig[1, 2]; ax_kw..., title=L"z=-5\,\text{m}")
    
    hideydecorations!.(ax_surface; ticks=false)

    ht_kw = (;
        colormap,
        colorrange=(-u_max, u_max),
        ht_kw...
    )
    
    ht = heatmap!(ax, xsᶠ / 1000, times / 3600, transpose(u); ht_kw...)
    ht_surface = heatmap!(ax_surface, xsᶠ / 1000, times / 3600, transpose(u_surface); ht_kw...)
    
    ct_kw = (;
        color=(:black, 0.5),
        levels=b_levels,
        ct_kw...
    )
    
    contour!(ax, xsᶜ / 1000, times / 3600, transpose(b); ct_kw...)
    contour!(ax_surface, xsᶜ / 1000, times / 3600, transpose(b_surface); ct_kw...)

    #contour!(ax, times / 3600, xsᶜ / 1000, b_offset; ct_kw...)
    #contour!(ax_surface, times / 3600, xsᶜ / 1000, b_surface_offset; ct_kw...)

    Colorbar(fig[1, 3], ht; label=L"(u + U) / \text{cm}\,\text{s}^{-1}")

    subfig_label!(fig[1, 1], 1)
    subfig_label!(fig[1, 2], 2)
    fig
end
