# pv_video.jl

# Potential vorticity and proportion of negative pv

function front_balance(
        foldername,
        frames; 
        fig_kw = (; ), 
        ax_kw = (; ),
        ht_kw = (; ),
        ct_kw = (; ),
        arr_kw = (; ),
        σ = 0
    )

    DFM = joinpath(foldername, "DFM.jld2")
    MBALANCE = joinpath(foldername, "MBALANCE.jld2")
    
    iterations, times = iterations_times(DFM)
    sp = simulation_parameters(DFM)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(DFM)
    inds = center_indices(DFM)
    
    colormap = to_colormap(:balance)
    
    fig_kw = (;
        size=(1000, 500),
        fig_kw...
    )
    
    fig = Figure(; fig_kw...)

    # Average over final state
    iterations = iterations[frames]
    times = times[frames]

    u = time_average_of(a->filt(a, σ), DFM, "u_dfm", iterations) .+ [-sp.α * x for x in xsᶠ, z in 1:1]
    u = (u[1:end-1, :] .+ u[2:end, :]) ./ 2
    w = time_average_of(a->filt(a, σ), DFM, "w_dfm", iterations)
    w = (w[:, 1:end-1] .+ w[:, 2:end]) ./ 2
    
    b = time_average_of(a->filt(a, σ), DFM, "b_dfm", iterations)

    # Arrow plots
    xs_u = xsᶜ[1:60:end] ./ 1000
    zs_u = zsᶜ[1:15:end]
    u = u[1:60:end, 1:15:end] ./ 1000
    w = w[1:60:end, 1:15:end]

    αΓ_arrest = time_average_of(a->filt(a, σ), MBALANCE, "M²", iterations)[500, 240] * sp.α

    αΓ = time_average_of(a->filt(a, σ), MBALANCE, "M²", iterations) .* sp.α ./ αΓ_arrest
    
    DΓDt = time_average_of(a->filt(a, σ), MBALANCE, "DM²Dt", iterations) ./ αΓ_arrest
    dΓ = get_field(a->filt(a, σ), MBALANCE, "M²", iterations[end]) ./ αΓ_arrest
    dΓ -= get_field(a->filt(a, σ), MBALANCE, "M²", iterations[1]) ./ αΓ_arrest
    ∂Γ∂t = dΓ ./ (times[end] - times[1])
    
    SC = -time_average_of(a->filt(a, σ), MBALANCE, "∂u∂xM²", iterations) ./ αΓ_arrest
    SC -= time_average_of(a->filt(a, σ), MBALANCE, "∂w∂xN²", iterations) ./ αΓ_arrest

    Fx = time_average_of(a->filt(a, σ), MBALANCE, "Fx", iterations) ./ αΓ_arrest
    Fz = time_average_of(a->filt(a, σ), MBALANCE, "Fz", iterations) ./ αΓ_arrest

    titles = [
        L"(\overline{u} - \alpha x)\frac{\partial M^2}{\partial x} + \overline{w}\frac{\partial M^2}{\partial z}",
        L"-\frac{\partial \overline{u}}{\partial x}M^2 - \frac{\partial \overline{w}}{\partial x}N^2",
        L"10\alpha M^2",
        L"F_{M^2}"
    ]
    
    ax_kw = (;
        xlabel = L"x/\text{km}",
        ylabel = L"z/\text{m}",
        limits = (-sp.Lh / 2000, sp.Lh / 2000, -sp.H, 0),
        ax_kw...
    )
    
    # Heatmaps
    ht_kw = (;
        colorrange = (-0.1, 0.1),
        colormap,
        highclip = colormap[end],
        lowclip = colormap[1],
        ht_kw...
    )

    ct_kw = (;
        color = :red,
        levels = 10,
        ct_kw...
    )

    arr_kw = (; 
        lengthscale=5000, 
        color = abs.(w)[:], 
        colormap = [RGBA(0, 0, 0, 0), RGBA(0, 0, 0, 1)],
        #colorrange = (0, maximum(abs, w) / 2),
        colorrange = (0, 100 / (24 * 3600)),
        arr_kw...
    )
    @info maximum(abs, w) ./ 2
    
    title = titles[1]
    ax = Axis(fig[1, 1]; ax_kw..., title)
    heatmap!(ax, xsᶠ ./ 1000, zsᶜ, DΓDt .- ∂Γ∂t; ht_kw...)
    #contour!(ax, xsᶠ ./ 1000, zsᶜ, ψ; ct_kw...)
    contour!(ax, xsᶜ ./ 1000, zsᶜ, b; levels=40, color=:black)
    arrows2d!(ax, xs_u, zs_u, u, w; arr_kw...)
    
    subfig_label!(fig[1, 1], 1)
    hidexdecorations!(ax; ticks=false)

    title = titles[2]
    ax = Axis(fig[1, 2]; ax_kw..., title)
    heatmap!(ax, xsᶠ ./ 1000, zsᶜ, SC; ht_kw...)
    #contour!(ax, xsᶠ ./ 1000, zsᶜ, ψ; ct_kw...)
    contour!(ax, xsᶜ ./ 1000, zsᶜ, b; levels=40, color=:black)
    arrows2d!(ax, xs_u, zs_u, u, w; arr_kw...)
    
    subfig_label!(fig[1, 2], 2)
    hidexdecorations!(ax; ticks=false)
    hideydecorations!(ax; ticks=false)

    title = titles[3]
    ax = Axis(fig[2, 1]; ax_kw..., title)
    heatmap!(ax, xsᶠ ./ 1000, zsᶜ, 10αΓ; ht_kw...)
    #contour!(ax, xsᶠ ./ 1000, zsᶜ, ψ; ct_kw...)
    contour!(ax, xsᶜ ./ 1000, zsᶜ, b; levels=40, color=:black)
    arrows2d!(ax, xs_u, zs_u, u, w; arr_kw...)
    
    subfig_label!(fig[2, 1], 3)

    title = titles[4]
    ax = Axis(fig[2, 2]; ax_kw..., title)
    ht = heatmap!(ax, xsᶠ ./ 1000, zsᶜ, Fx .+ Fz; ht_kw...)
    #contour!(ax, xsᶠ ./ 1000, zsᶜ, ψ; ct_kw...)
    contour!(ax, xsᶜ ./ 1000, zsᶜ, b; levels=40, color=:black)
    arrows2d!(ax, xs_u, zs_u, u, w; arr_kw...)
    
    subfig_label!(fig[2, 2], 4)
    hideydecorations!(ax; ticks=false)
    
    Colorbar(fig[1:2, 3], ht)

    fig
end
