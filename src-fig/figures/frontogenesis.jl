function frontogenesis(
        foldername;
        fig_kw=(; ), 
        ax_kw=(; ),
        ht_kw=(; ),
        ct_kw=(; ),
        σ=0,
    )

    iterations, times = iterations_times(foldername)
    sp = simulation_parameters(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    inds = centre_indices(foldername)
    colormap = to_colormap(:balance)
    z_indᶜ = zᶜbounds(foldername, -sp.H)
    
    fig = Figure(; 
        size=(800, 400),
        fig_kw...
    )

    v_dfm = timeseries_of(a->a[:, z_indᶜ:end], joinpath(foldername, "DFM.jld2"), "v_dfm", iterations)
    ζ_dfm = zeros(Float64, length(times), length(xsᶠ), length(zsᶜ[z_indᶜ:end]))
    for i in 2:(size(ζ_dfm, 2)-1)
        ζ_dfm[:, i, :] .= (v_dfm[:, i, :] - v_dfm[:, i-1, :]) / (xsᶜ[i] - xsᶜ[i-1])
    end
    ζ_dfm .= filt(ζ_dfm ./ sp.f, σ)
    
    b_dfm = timeseries_of(a->a[:, z_indᶜ:end], joinpath(foldername, "DFM.jld2"), "b_dfm", iterations)
    b_dfm = filt(b_dfm .+ (b_dfm[1, inds[1:1], 1:1] .- b_dfm[:, inds[1:1], 1:1]), σ)
    
    ζ_max = maximum(ζ_dfm[:, :, end])
    ζ_min = minimum(ζ_dfm[:, :, end])

    ax_kw = (;
        xlabel=L"t/\text{hr}",
        ylabel=L"x/\text{km}",
        yticks=[-2, -1, 0, 1, 2],
        limits=(0, 200, -sp.ℓ / 1000, sp.ℓ / 1000),
        ax_kw...
    )
    ax = Axis(fig[1, 1]; ax_kw..., title="Depth average")
    ax_surface = Axis(fig[1, 2]; ax_kw..., title="Surface")
    
    hideydecorations!.(ax_surface; ticks=false)

    ht_kw = (;
        colormap,
        colorrange=(-max(ζ_max, -ζ_min), max(ζ_max, -ζ_min)),
        highclip=ζ_max > ht_kw.colorrange[2] ? colormap[end] : nothing,
        lowclip=nothing, #ζ_min < ht_kw.colorrange[1] ? colormap[end] : nothing,
        ht_kw...
    )
    
    ht = heatmap!(ax, times / 3600, xsᶠ / 1000, mean(ζ_dfm; dims=3)[:, :, 1]; ht_kw...)
    ht_surface = heatmap!(ax_surface, times / 3600, xsᶠ / 1000, ζ_dfm[:, :, end]; ht_kw...)
    
    ct_kw = (;
        color=(:black, 0.5),
        levels=b_levels,
        ct_kw...
    )
    
    contour!(ax, times / 3600, xsᶜ / 1000, mean(b_dfm; dims=3)[:, :, 1]; ct_kw...)
    contour!(ax_surface, times / 3600, xsᶜ / 1000, b_dfm[:, :, end]; ct_kw...)

    Colorbar(fig[1, 3], ht; label=L"\zeta / f")

    subfig_label!(fig[1, 1], 1)
    subfig_label!(fig[1, 2], 2)
    fig
end