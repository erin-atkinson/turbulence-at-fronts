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

    iterations, times = iterations_times(foldername)
    sp = simulation_parameters(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    inds = centre_indices(foldername)
    
    term_names = ["VSP", "LSP", "BFLUX", "DSP", "ε"]
    term_labels = [L"\text{VSP}", L"\text{LSP}", L"\text{BFLUX}", L"\text{DSP}'", L"-\varepsilon"]

    # Plot in front and in arrest region
    
    fig = Figure(; size=(800, 250), fig_kw...)

    region_scenes = (; arrest=fig[1, 2], total=fig[1, 1], )

    TKE = joinpath(foldername, "TKE.jld2")

    Δm = 1.027 * sp.Lh * sp.Lz * sp.Ly / (sp.Nh * sp.Nz)
    Δt = 3600 * 24

    # Reference energy 
    ΔE = sum(get_field(v->Δm * v.^2 / 2, joinpath(foldername, "DFM.jld2"), "v_dfm", iterations[1]))

    terms = map(regions) do region
        mask = [maskfromlines(x, z, region) for x in xsᶜ, z in zsᶜ]
        terms = map(term_names[1:end-1]) do term_name
            timeseries_of(a->sum(mask .* a), TKE, term_name, iterations) * Δm * Δt / ΔE
        end
        ε = (timeseries_of(a->sum(mask .* a), TKE, "DTKEDt", iterations) * Δm * Δt / ΔE) .- sum(terms)
        [terms; [-ε]]
    end
    termmax = mapreduce(max, terms.total) do term
        maximum(abs, term)
    end
    
    ax_kw = (;
        xlabel=L"t / \text{hr}",
        ylabel=L"\Delta P / \text{day}^{-1}",
        limits=(0, nothing, -termmax, termmax),
        ax_kw...
    )
    
    axeslns = map(region_names, region_scenes, terms) do title, scene, ys
        ax, lns = plot_terms(scene, (; title, ax_kw...), ln_kw, times/3600, ys)
        (; ax, lns)
    end

    Legend(fig[1, 3], axeslns.arrest.lns, term_labels, L"E_0 \Delta P")
    
    #hidexdecorations!(axeslns.arrest.ax; ticks=false, grid=false)
    #hidexdecorations!(axeslns.top.ax; ticks=false, grid=false)
    #hideydecorations!(axeslns.top.ax; ticks=false, grid=false)
    hideydecorations!(axeslns.arrest.ax; ticks=false, grid=false)
    subfig_label!(fig[1, 1], 1)
    subfig_label!(fig[1, 2], 2)

    
    for frame in frames
        lines!(axeslns.total.ax, [times[frame], times[frame]] / 3600, [-termmax, termmax]; linestyle=:dash, color=:grey)
        lines!(axeslns.arrest.ax, [times[frame], times[frame]] / 3600, [-termmax, termmax]; linestyle=:dash, color=:grey)
    end
    fig
end

function frontogenesis_by_region(
        foldername;
        fig_kw=(; ), 
        ax_kw=(; ),
        ln_kw=(; )
    )

    iterations, times = iterations_times(foldername)
    sp = simulation_parameters(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    inds = centre_indices(foldername)
    
    term_names = [
        "∇bD∇bDt_dfm",
        "background_strain_dfm",
        "strain_dfm",
        "shear_dfm",
        "divergence_dfm",
        "tilting_dfm",
        "subgrid_dfm"
    ]
    term_labels = [
        L"\frac{D}{Dt}\frac{1}{2}|\nabla b|^2",
        L"\alpha_0 (b_{,x}^2 - b_{,y}^2)",
        L"\alpha (b_{,x}^2 - b_{,y}^2)",
        L"-2\tau (b_{,x}b_{,y})",
        L"-\delta (b_{,x}^2 + b_{,y}^2)",
        L"-(w_{,x}b_{,x} + w_{,y}b_{,y})b_{,z}",
        L"\text{sgs}"
    ]

    # Plot in front and in arrest region
    
    fig = Figure(; size=(800, 250), fig_kw...)

    region_scenes = (; arrest=fig[1, 2], total=fig[1, 1], )

    TKE = joinpath(foldername, "FRONTOGENESIS.jld2")
    terms = map(regions) do region
        mask = [maskfromlines(x, z, region) for x in xsᶜ, z in zsᶜ]
        map(term_names) do term_name
            timeseries_of(a->1e12 * sum(mask .* a), TKE, term_name, iterations)
        end
    end

    ax_kw = (;
        xlabel=L"t / \text{hr}",
        ylabel=L"10^12\Delta / \text{m}^2\text{s}^{-3}",
        limits=(0, nothing, nothing, nothing),
        ax_kw...
    )
    
    axeslns = map(region_names, region_scenes, terms) do title, scene, ys
        ax, lns = plot_terms(scene, (; title, ax_kw...), ln_kw, times/3600, ys)
        (; ax, lns)
    end

    Legend(fig[1, 3], axeslns.arrest.lns, term_labels; title=L"\Delta")
    
    #hidexdecorations!(axeslns.arrest.ax; ticks=false, grid=false)
    #hidexdecorations!(axeslns.top.ax; ticks=false, grid=false)
    #hideydecorations!(axeslns.top.ax; ticks=false, grid=false)
    #hideydecorations!(axeslns.underneath.ax; ticks=false, grid=false)
    hideydecorations!(axeslns.arrest.ax; ticks=false, grid=false, ticklabels=false)
    
    fig
end

function mean_frontogenesis_by_region(
        foldername;
        fig_kw=(; ), 
        ax_kw=(; ),
        ln_kw=(; )
    )

    iterations, times = iterations_times(foldername)
    sp = simulation_parameters(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    inds = centre_indices(foldername)
    
    term_names = [
        "background_strain_dfm",
        "divergence_dfm",
        "tilting_dfm",
        "turbulence_h_dfm",
        "turbulence_z_dfm",
        "subgrid_dfm"
    ]
    term_labels = [
        L"\alpha (b_{,x}^2 - b_{,y}^2)",
        L"-\delta (b_{,x}^2 + b_{,y}^2)",
        L"-(w_{,x}b_{,x} + w_{,y}b_{,y})b_{,z}",
        L"\text{turbulence_x}",
        L"\text{turbulence_z}",
        L"\text{sgs}"
    ]

    # Plot in front and in arrest region
    
    fig = Figure(; size=(800, 250), fig_kw...)

    region_scenes = (; arrest=fig[1, 2], total=fig[1, 1], )

    TKE = joinpath(foldername, "MEANFRONTOGENESIS.jld2")
    terms = map(regions) do region
        mask = [maskfromlines(x, z, region) for x in xsᶜ, z in zsᶜ]
        map(term_names) do term_name
            filt(timeseries_of(a->1e14 * sum(mask .* a), TKE, term_name, iterations), 16)
        end
    end

    ax_kw = (;
        xlabel=L"t / \text{hr}",
        ylabel=L"10^14\Delta / \text{m}^2\text{s}^{-3}",
        limits=(0, 200, -5, 5),
        ax_kw...
    )
    
    axeslns = map(region_names, region_scenes, terms) do title, scene, ys
        ax, lns = plot_terms(scene, (; title, ax_kw...), ln_kw, times/3600, ys)
        (; ax, lns)
    end

    Legend(fig[1, 3], axeslns.arrest.lns, term_labels; title=L"\Delta")
    
    #hidexdecorations!(axeslns.arrest.ax; ticks=false, grid=false)
    #hidexdecorations!(axeslns.top.ax; ticks=false, grid=false)
    #hideydecorations!(axeslns.top.ax; ticks=false, grid=false)
    #hideydecorations!(axeslns.underneath.ax; ticks=false, grid=false)
    hideydecorations!(axeslns.arrest.ax; ticks=false, grid=false, ticklabels=false)
    
    fig
end

TW(px, fv, u) = u .* (px + fv)
TTW(px, fv, w′u′, u) = u .* (px + fv + w′u′)

function ubalance_by_region(
        foldername;
        fig_kw=(; ), 
        ax_kw=(; ),
        ln_kw=(; )
    )

    iterations, times = iterations_times(foldername)
    sp = simulation_parameters(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    inds = centre_indices(foldername)
    
    term_names = [
        "background_strain_dfm",
        "pressure_dfm",
        "coriolis_dfm",
        "turbulence_h_dfm",
        "turbulence_z_dfm",
        "subgrid_dfm"
    ]
    term_labels = [
        L"\alpha u^2",
        L"fvu-p_x u",
        L"\text{turbulence_x}",
        L"\text{turbulence_z}",
        L"\text{sgs}"
    ]

    # Plot in front and in arrest region
    
    fig = Figure(; size=(800, 250), fig_kw...)

    region_scenes = (; arrest=fig[1, 2], total=fig[1, 1], )

    UBALANCE = joinpath(foldername, "UBALANCE.jld2")
    terms = map(regions) do region
        mask = [maskfromlines(x, z, region) for x in xsᶠ, z in zsᶜ]
        terms = map(term_names) do term_name
            filt(timeseries_of(a->sum(mask .* a) * 1e9, UBALANCE, term_name, iterations), 1) ./ sum(mask)
        end
        [terms[1], filt(terms[2] + terms[3], 20), terms[4], terms[5], terms[6]]
    end
    
    ax_kw = (;
        xlabel=L"t / \text{hr}",
        ylabel=L"10^9 \Delta/ \text{m}\text{s}^{-2}",
        limits=(0, 100, nothing, nothing),
        ax_kw...
    )
    
    axeslns = map(region_names, region_scenes, terms) do title, scene, ys
        ax, lns = plot_terms(scene, (; title, ax_kw...), ln_kw, times/3600, ys)
        (; ax, lns)
    end

    Legend(fig[1, 3], axeslns.arrest.lns, term_labels; title=L"\Delta")
    
    #hidexdecorations!(axeslns.arrest.ax; ticks=false, grid=false)
    #hidexdecorations!(axeslns.top.ax; ticks=false, grid=false)
    #hideydecorations!(axeslns.top.ax; ticks=false, grid=false)
    #hideydecorations!(axeslns.underneath.ax; ticks=false, grid=false)
    hideydecorations!(axeslns.arrest.ax; ticks=false, grid=false, ticklabels=false)
    
    fig
end

function tke_by_term(
        foldername;
        fig_kw=(; ), 
        ax_kw=(; ),
        ln_kw=(; )
    )

    iterations, times = iterations_times(foldername)
    sp = simulation_parameters(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    inds = centre_indices(foldername)
    
    term_names = ["VSP", "LSP", "BFLUX", "ε"]
    term_labels = [L"\text{VSP}", L"\text{LSP}", L"\text{BFLUX}", L"\varepsilon"]

    fig = Figure(; size=(1000, 500), fig_kw...)

    term_scenes = [fig[1, 1], fig[1, 2], fig[2, 1], fig[2, 2]]
    
    TKE = joinpath(foldername, "TKE.jld2")
    terms = map(term_names) do term_name
        map(regions) do region
            mask = [maskfromlines(x, z, region) for x in xsᶜ, z in zsᶜ]
            timeseries_of(a->sum(mask .* a) ./ sum(mask), TKE, term_name, iterations)
        end
    end

    ax_kw = (;
        xlabel=L"t / \text{hr}",
        ylabel=L"\Delta E / \text{m}^2\text{s}^{-3}",
        limits=(0, nothing, -1e-7, 1e-7)
    )
    
    axeslns = map(term_names, term_scenes, terms) do title, scene, ys
        ax, lns = plot_terms(scene, (; title, ax_kw...), ln_kw, times/3600, ys)
        (; ax, lns)
    end

    Legend(fig[1:2, 3], [ln for ln in axeslns[1].lns], [region for region in region_names])
    
    hidexdecorations!(axeslns[1].ax; ticks=false, grid=false)
    hidexdecorations!(axeslns[2].ax; ticks=false, grid=false)
    #hideydecorations!(axeslns.top.ax; ticks=false, grid=false)
    #hideydecorations!(axeslns.underneath.ax; ticks=false, grid=false)
    
    fig
end

function pv_by_term(
        foldername;
        fig_kw=(; ), 
        ax_kw=(; ),
        ln_kw=(; )
    )

    iterations, times = iterations_times(foldername)
    sp = simulation_parameters(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
    inds = centre_indices(foldername)
    
    term_names = ["q_dfm", "uq_dfm", "wq_dfm", "div_Uq_dfm"]
    term_labels = [L"q", L"\overline{uq}", L"\overline{wq}", L"\overline{\nabla \cdot \vec{u}q}"]
    
    fig = Figure(; size=(1000, 500), fig_kw...)

    term_scenes = [fig[1, 1], fig[1, 2], fig[2, 1], fig[2, 2]]
    
    PV = joinpath(foldername, "PV.jld2")
    terms = map(term_names) do term_name
        map(regions) do region
            mask = [maskfromlines(x, z, region) for x in xsᶜ, z in zsᶜ]
            filt(timeseries_of(a->sum(mask .* a[1:1024, 1:128]) ./ sum(mask), PV, term_name, iterations), 1)
        end
    end

    ax_kw = (;
        xlabel=L"t / \text{hr}",
        ylabel=L"A"
    )
    
    axeslns = map(term_labels, term_scenes, terms) do title, scene, ys
        ax, lns = plot_terms(scene, (; title, ax_kw...), ln_kw, times/3600, ys)
        (; ax, lns)
    end

    Legend(fig[1:2, 3], [ln for ln in axeslns[1].lns], [region for region in region_names])
    
    hidexdecorations!(axeslns[1].ax; ticks=false, grid=false)
    hidexdecorations!(axeslns[2].ax; ticks=false, grid=false)
    #hideydecorations!(axeslns.top.ax; ticks=false, grid=false)
    #hideydecorations!(axeslns.underneath.ax; ticks=false, grid=false)
    
    fig
end
