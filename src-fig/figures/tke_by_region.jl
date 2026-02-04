function plot_terms(scene, ax_kw, ln_kw, x, ys)
    ax = Axis(scene; ax_kw...)
    lns = map(ys) do y
        lines!(ax, x, y; ln_kw...)
    end
    return ax, lns
end

function tke_by_region(
        foldername,
        frames=[];
        fig_kw=(; ), 
        ax_kw=(; ),
        ln_kw=(; )
    )

    TKE = joinpath(foldername, "TKE.jld2")
    
    iterations, times = iterations_times(TKE)
    sp = simulation_parameters(TKE)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(TKE)
    inds = center_indices(TKE)
    
    term_labels = [L"\text{VSP}", L"\text{GSP}", L"\text{LSP}", L"\text{BFLUX}", L"\text{DSP}'", L"\varepsilon"]
    colors = Makie.wong_colors()
    
    # Plot in front and in arrest region
    
    fig = Figure(; size=(800, 250), fig_kw...)

    Δm = 1.027 * sp.Lh * sp.Lz * sp.Ly / (sp.Nh * sp.Nz)
    Δt = 3600 * 24

    # Reference energy 
    u = get_field(v->v.^2 / 2, joinpath(foldername, "DFM.jld2"), "u_dfm", iterations[21])
    v = get_field(v->v.^2 / 2, joinpath(foldername, "DFM.jld2"), "v_dfm", iterations[21])
    w = get_field(v->v.^2 / 2, joinpath(foldername, "DFM.jld2"), "w_dfm", iterations[21])
    tke = get_field(joinpath(foldername, "TKE.jld2"), "TKE", iterations[21])

    ΔE = sum(u) + sum(v) + sum(w) + sum(tke)
    
    #ΔE = get_field(v->sum(v), joinpath(foldername, "TKE.jld2"), "TKE", iterations[21])
    
    mask = [maskfromlines(x, z, regions.total) for x in xsᶜ, z in zsᶜ]
    VSP = timeseries_of(a->sum(mask .* a), TKE, "VSP", iterations) * Δt / ΔE
    GSP = timeseries_of(a->sum(mask .* a), TKE, "GSP", iterations) * Δt / ΔE
    LSP = timeseries_of(a->sum(mask .* a), TKE, "LSP", iterations) * Δt / ΔE
    BFLUX = timeseries_of(a->sum(mask .* a), TKE, "BFLUX", iterations) * Δt / ΔE
    DSP′ = timeseries_of(a->sum(mask .* a), TKE, "DSP", iterations) * Δt / ΔE

    tke = timeseries_of(a->sum(mask .* a), joinpath(foldername, "TKE.jld2"), "TKE", iterations) * Δt / ΔE
    ε = diff([0; tke]) ./ diff([times[1] - (times[2] - times[1]); times]) .- (VSP .+ LSP .+ BFLUX .+ DSP′)

    mask = [maskfromlines(x, z, regions.arrest) for x in xsᶜ, z in zsᶜ]
    VSP_arrest = timeseries_of(a->sum(mask .* a), TKE, "VSP", iterations) * Δt / ΔE
    GSP_arrest = timeseries_of(a->sum(mask .* a), TKE, "GSP", iterations) * Δt / ΔE
    LSP_arrest = timeseries_of(a->sum(mask .* a), TKE, "LSP", iterations) * Δt / ΔE
    BFLUX_arrest = timeseries_of(a->sum(mask .* a), TKE, "BFLUX", iterations) * Δt / ΔE
    DSP′_arrest = timeseries_of(a->sum(mask .* a), TKE, "DSP", iterations) * Δt / ΔE
    #ε_arrest = timeseries_of(a->sum(mask .* a), TKE, "ε", iterations) * Δt / ΔE
    
    ax_kw = (;
        xlabel=L"t / \text{hr}",
        ylabel=L"\Delta P / \text{day}^{-1}",
        title = L"\text{Full domain}",
        limits=(0, 160, -2, 2),
        ax_kw...
    )

    ax = Axis(fig[1, 1]; ax_kw...)

    lns = [
        lines!(ax, times ./ 3600, VSP; ln_kw..., color=(colors[1], 0.5)),
        lines!(ax, times ./ 3600, GSP; ln_kw..., color=colors[1], linestyle=:dash),
        lines!(ax, times ./ 3600, LSP; ln_kw..., color=colors[2]),
        lines!(ax, times ./ 3600, BFLUX; ln_kw..., color=colors[3]),
        lines!(ax, times ./ 3600, DSP′; ln_kw..., color=colors[4]),
        lines!(ax, times ./ 3600, ε; ln_kw..., color=colors[5]),
    ]

    ax_kw = (;
        xlabel=L"t / \text{hr}",
        ax_kw...,
        limits=(0, 160, -0.5, 0.5),
        title = L"\text{Arrest region}",
        ylabel = ""
    )
    
    ax_arrest = Axis(fig[1, 2]; ax_kw..., )
    
    lines!(ax_arrest, times ./ 3600, VSP_arrest; ln_kw..., color=(colors[1], 0.5))
    lines!(ax_arrest, times ./ 3600, GSP_arrest; ln_kw..., color=colors[1], linestyle=:dash)
    lines!(ax_arrest, times ./ 3600, LSP_arrest; ln_kw..., color=colors[2])
    lines!(ax_arrest, times ./ 3600, BFLUX_arrest; ln_kw..., color=colors[3])
    lines!(ax_arrest, times ./ 3600, DSP′_arrest; ln_kw..., color=colors[4])
    #lines!(ax_arrest, times ./ 3600, ε_arrest; ln_kw..., color=colors[5])
    
    Legend(fig[1, 3], lns, term_labels, L"E_0 \Delta P")
    
    #hideydecorations!(ax_arrest; ticks=false, grid=false)
    subfig_label!(fig[1, 1], 1)
    subfig_label!(fig[1, 2], 2)

    for frame in frames
        lines!(ax, [times[frame], times[frame]] / 3600, [-10, 10]; linestyle=:dash, color=:grey)
        lines!(ax_arrest, [times[frame], times[frame]] / 3600, [-10, 10]; linestyle=:dash, color=:grey)
    end
    fig
end
