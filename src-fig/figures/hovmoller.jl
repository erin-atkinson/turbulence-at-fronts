function hovmoller(
        foldername,
        z;
        fig_kw=(; ), 
        ax_kw=(; ),
        ht_kw=(; ),
        ct_kw=(; ),
        σ=0,
        background=false,
    )

    iterations, times = iterations_times(foldername)
    sp = simulation_parameters(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    inds = centre_indices(foldername)
    colormap = to_colormap(:balance)
    z_indᶜ = zᶜbounds(foldername, z)
    
    fig = Figure(; 
        size=(800, 400),
        fig_kw...
    )

    U = [-variable_strain_rate(t, sp) * x for t in times, x in xsᶠ] .* background
    
    u = 100 * (timeseries_of(a->filt(a[:, z_indᶜ], σ), joinpath(foldername, "DFM.jld2"), "u_dfm", iterations) .+ U)
    b = timeseries_of(a->a[:, z_indᶜ], joinpath(foldername, "DFM.jld2"), "b_dfm", iterations)
    
    u_surface = 100 * (timeseries_of(a->filt(a[:, end], σ), joinpath(foldername, "DFM.jld2"), "u_dfm", iterations))
    b_surface = timeseries_of(a->a[:, end], joinpath(foldername, "DFM.jld2"), "b_dfm", iterations)

    b_offset = b[1:1, inds[1:1]] .- b[:, inds[1:1]]
    
    b = filt(b .+ b_offset, σ)
    b_surface = filt(b_surface .+ b_offset, σ)
    
    u_max = max(maximum(abs, u[:, inds]), maximum(abs, u_surface[:, inds]))

    ax_kw = (;
        xlabel=L"t/\text{hr}",
        ylabel=L"x/\text{km}",
        yticks=[-2, -1, 0, 1, 2],
        limits=(0, 200, -sp.ℓ / 1000, sp.ℓ / 1000),
        ax_kw...
    )
    ax = Axis(fig[1, 1]; ax_kw..., title=L"z=%$z~\text{m}")
    ax_surface = Axis(fig[1, 2]; ax_kw..., title="Surface")
    
    hideydecorations!.(ax_surface; ticks=false)

    ht_kw = (;
        colormap,
        colorrange=(-u_max, u_max),
        ht_kw...
    )
    
    ht = heatmap!(ax, times / 3600, xsᶠ / 1000, u; ht_kw...)
    ht_surface = heatmap!(ax_surface, times / 3600, xsᶠ / 1000, u_surface; ht_kw...)
    
    ct_kw = (;
        color=(:black, 0.5),
        levels=b_levels,
        ct_kw...
    )
    
    contour!(ax, times / 3600, xsᶜ / 1000, b; ct_kw...)
    contour!(ax_surface, times / 3600, xsᶜ / 1000, b_surface; ct_kw...)

    Colorbar(fig[1, 3], ht; label=L"u / \text{cm}\,\text{s}^{-1}")

    subfig_label!(fig[1, 1], 1)
    subfig_label!(fig[1, 2], 2)
    fig
end