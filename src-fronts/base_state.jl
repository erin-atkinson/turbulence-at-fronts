using SpecialFunctions
using OffsetArrays: no_offset_view

@inline g(s) = s < 0 ? 0 : sqrt(1 + s^2) - 1
#@inline g(s) = log(1 + exp(s))
@inline G(s) = s < 0 ? 0 : 0.5 * (s * (-2 + sqrt(1 + s^2)) + asinh(s))

@inline γ(s, δ) = -1 + δ * (erf(s) + 1) / 2
@inline γ′(s, δ) = δ * exp(-s^2) / sqrt(π)
@inline A(δ) = maximum(range(-3, 3, 1000)) do s
    γ(s, δ) * (γ′(s+5e-7, δ) - γ′(s-5e-7, δ)) / 1e-6 + γ′(s, δ)^2
end

@inline B(H, a, ℓ, δ) = let h(x) = H * γ(x/ℓ, δ)
    integrated_h(x) = h(x) - h(x + a*h(x+a*h(x+a*h(x+a*h(x+a*h(x+a*h(x))))))) # This is overkill
    integrated_h′(x) = (integrated_h(x+5e-5ℓ) - integrated_h(x-5e-5ℓ)) / 1e-4ℓ
    maximum(integrated_h′, range(-2ℓ, 2ℓ, 100))
end
#= 
The base/reference state is a front above a mixed layer, with the varying mixed
layer height providing the 
=#
#@inline function get_base_state_filament(simulation_parameters) end
#=
@inline function get_base_state_front(simulation_parameters)
    sp = simulation_parameters
    # Create the front
    
    @inline h(x) = sp.H * γ(x/sp.ℓ, sp.δ)
    @inline h′(x) = sp.H * γ′(x/sp.ℓ, sp.δ) / sp.ℓ
    
    @inline function b(x, z)
        sp.N₀² * z + sp.ΔN² * sp.λ * sp.H * g((z - h(x))/(sp.λ * sp.H))
    end
    @inline function v(x, z)
        -sp.ΔN² * sp.λ * sp.H * g((z - h(x))/(sp.λ * sp.H)) * h′(x) / sp.f
    end
    
    @inline w₀(x, y, z) = 0
    @inline b₀(x, y, z) = b(x, z)
    @inline u₀(x, y, z) = 0
    @inline v₀(x, y, z) = v(x, z)
    
    (; b=b₀, u=u₀, v=v₀, w=w₀)
end
=#
# This one is more complicated, so it's going to return an array
@inline function get_base_state_front(grid, simulation_parameters; with_halos=false)
    sp = simulation_parameters
    @inline h(x) = sp.H * γ(x/sp.ℓ, sp.δ)
    @inline h′(x) = sp.H * γ′(x/sp.ℓ, sp.δ) / sp.ℓ
    @inline function b(x, z)
        sp.N₀² * (z -(1-sp.ε) * sp.λ * sp.H * g((z - h(x + sp.a * (z + sp.H/2)))/(sp.λ * sp.H)))
    end
    
    N²(x, z) = (b(x, z+5e-8) - b(x, z-5e-8)) / 1e-7
    M²(x, z) = (b(x+5e-8, z) - b(x-5e-8, z)) / 1e-7
    
    # Now compute the appropriate velocity
    xsᶜ, ysᶜ, zsᶜ = nodes(grid, Center(), Center(), Center(); with_halos)
    xsᶠ, ysᶠ, zsᶠ = nodes(grid, Face(), Face(), Face(); with_halos)
    (xsᶜ, ysᶜ, zsᶜ, xsᶠ, ysᶠ, zsᶠ) = no_offset_view.((xsᶜ, ysᶜ, zsᶜ, xsᶠ, ysᶠ, zsᶠ))
    vs = cumsum([(zsᶜ[i] - zsᶜ[max(i-1, 1)]) * M²(x, z) / sp.f for x in xsᶜ, y in ysᶜ, (i, z) in enumerate(zsᶜ)]; dims=3)
    bs = [b(x, z) for x in xsᶜ, y in ysᶜ, z in zsᶜ]
    return (; u=zeros(length(xsᶠ), length(ysᶜ), length(zsᶜ)), w=zeros(length(xsᶜ), length(ysᶜ), length(zsᶠ)), v=vs, b=bs)
end

@inline function get_base_state_front(simulation_parameters)
    sp = simulation_parameters
    @inline h(x) = sp.H * γ(x/sp.ℓ, sp.δ)
    @inline h′(x) = sp.H * γ′(x/sp.ℓ, sp.δ) / sp.ℓ
    @inline function b(x, z)
        sp.N₀² * (z -(1-sp.ε) * sp.λ * sp.H * g((z - h(x + sp.a * (z + sp.H/2)))/(sp.λ * sp.H)))
    end
    
    N²(x, z) = (b(x, z+5e-8) - b(x, z-5e-8)) / 1e-7
    M²(x, z) = (b(x+5e-8, z) - b(x-5e-8, z)) / 1e-7
    
    # Now compute the appropriate velocity
    xs = range(-5sp.ℓ, 5sp.ℓ, 500)
    zs = range(-sp.H, 0, 500)
    
    vs = cumsum([(zs[i] - zs[max(i-1, 1)]) * M²(x, z) / sp.f for x in xs, (i, z) in enumerate(zs)]; dims=2)
    bs = [b(x, z) for x in xs, z in zs]
    
    return (; xs, zs, u=(x, y, z)->0, w=(x, y, z)->0, v=vs, b=bs)
end

@inline get_base_state(simulation_parameters) = get_base_state_front(simulation_parameters)
@inline get_base_state(grid, simulation_parameters) = get_base_state_front(grid, simulation_parameters)