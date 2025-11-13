using Oceananigans
include("../src-simulation/grid_faces.jl")
include("../src-simulation/base_state.jl")

@inline function σ(z, sp)
    s = min((z+sp.Lz) / (sp.Lz-sp.H), 1)
    return sp.σ * (1 - abs(s))^2
end

function preinitial_conditions(sp)
    # Pre-initial conditions
    xs, _, _ = get_grid_faces(sp)
    zs = range(-sp.Lz, 0, sp.Nz)

    vs = [approximate_front_velocity(x, z, sp) for x in xs, z in zs]
    bs = [front_buoyancy(x, z, sp) for x in xs, z in zs]
    σs = [σ(z, sp) for x in xs, z in zs] .* 3600
    hs = [-sp.H for x in xs]

    levels = minimum(bs):(sp.Δb / 5):maximum(bs)
    
    fig = Figure(; size=(600, 160))

    ax = Axis(fig[1, 1]; 
        xlabel=L"x / \text{km}", 
        ylabel=L"z / \text{m}", 
        limits=(xs[1] / 1000, xs[end] / 1000, -sp.Lz, 0),
        xminorticks=xs[1:32:end] ./ 1000,
        xticks=round.([-sp.Lx/2, -sp.Lh/2, sp.Lh/2, sp.Lx/2] ./ 1000; digits=1),
        xminorticksvisible=true
    )
    σ_color = to_colormap(:tempo)[200]
    σ_colormap = [RGBA(σ_color.r, σ_color.g, σ_color.b, 0), σ_color]

    ct = contourf!(ax, xs ./ 1000, zs, vs; levels=range(0, maximum(abs, vs), 5), colormap=to_colormap(:amp)[1:200])
    contour!(ax, xs ./ 1000, zs, bs; levels, color=(:black, 0.5))
    #lines!(ax, xs ./ 1000, hs; color=(:blue, 0.5), linestyle=:dash)
    lines!(ax, [-sp.Lh/2000, -sp.Lh/2000, sp.Lh/2000, sp.Lh/2000, -sp.Lh/2000], [0, -sp.H, -sp.H, 0, 0]; color=:magenta, linestyle=:dashdot)
    ht = heatmap!(ax, xs ./ 1000, zs, σs; colormap=σ_colormap)
    
    Colorbar(fig[1, 2], ct; label=L"v / \text{m}\,\text{s}^{-1}")
    Colorbar(fig[1, 3], ht; label=L"\sigma / \text{hr}^{-1}")
    fig
end
