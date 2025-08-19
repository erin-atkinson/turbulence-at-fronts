# base_state.jl
# Functions describing the initial state of the simulations

using SpecialFunctions
using Oceananigans: fill_halo_regions!

#@inline g(s) = s < 0 ? 0 : sqrt(1 + s^2) - 1
@inline g(s) = log(1 + exp(s))

@inline γ(s, δ) = -1 + δ * (erf(s) + 1) / 2
@inline γ′(s, δ) = δ * exp(-s^2) / sqrt(π)
@inline A(δ) = maximum(range(-3, 3, 1000)) do s
    γ(s, δ) * (γ′(s+5e-7, δ) - γ′(s-5e-7, δ)) / 1e-6 + γ′(s, δ)^2
end

@inline mixed_layer_depth(x, sp) = sp.H * γ(x/sp.ℓ, sp.δ)

# Actual bottom boundary is given by fixed point
@inline h₀(x, sp) = h₀(x, sp, Val(10))
@inline h₀(x, sp, ::Val{0}) = -sp.H
@inline h₀(x, sp, ::Val{T}) where T = mixed_layer_depth(x + sp.a * (h₀(x, sp, Val(T-1)) + sp.H/2), sp)

@inline function front_buoyancy(x, z, sp)
    # Well-mixed, Ri achieved by tilting isopycnals
    MLD = mixed_layer_depth(x + sp.a * (z + sp.H/2), sp)
    return sp.N₀² * (z - sp.λ * sp.H * g((z - MLD)/(sp.λ * sp.H)))
end

@inline function approximate_front_velocity(x, z, sp)
    # We can approximate the velocity far above the thermocline (-(λ + δ)H << z)
    MLD = mixed_layer_depth(x + sp.a * (z + sp.H/2), sp)
    return sp.N₀² * max(MLD - h₀(x, sp), 0) / sp.f / sp.a
end

@inline function approximate_front_vorticity(x, z, sp)
    return (approximate_front_velocity(x + 5e-6sp.ℓ, 0, sp) - approximate_front_velocity(x - 5e-6sp.ℓ, 0, sp)) / (1e-5sp.ℓ)
end

@inline function front_M²(x, z, sp)
    return (front_buoyancy(x + 5e-6sp.ℓ, z, sp) - front_buoyancy(x - 5e-6sp.ℓ, z, sp)) / (1e-5sp.ℓ)
end

@inline front_shear(x, z, sp) = front_M²(x, z, sp) / sp.f
    
@inline function front_N²(x, z, sp)
    return (front_buoyancy(x, z + 5e-6sp.H, sp) - front_buoyancy(x, z - 5e-6sp.H, sp)) / 1e-5sp.H
end

@inline front_Ri(x, z, sp) = front_N²(x, z, sp) / front_shear(x, z, sp)^2

@inline function minimum_Ri(sp)
    # Minimum Ri occurs in the centre of the front throughout the mixed layer
    return minimum(x->front_Ri(x, 0, sp), range(-3sp.ℓ, 3sp.ℓ, 1000))
end

@inline function minimum_Ro(sp)
    # Maximum Rossby number occurs at the surface
    return minimum(range(-3sp.ℓ, 3sp.ℓ, 1000)) do x
        approximate_front_vorticity(x, 0, sp) / sp.f
    end
end

@inline function front_initial_conditions(grid::RectilinearGrid, sp)
    # Use Oceananigans fields to setup the initial thermal wind properly
    
    b = Field{Center, Center, Center}(grid)
    set!(b, (x, y, z)->front_buoyancy(x, z, sp))
    fill_halo_regions!(b)

    # Compute the thermal wind shear
    S_op = @at (Center, Face, Center) ∂x(b) / sp.f
    S = compute!(Field(B_op))

    # Integrate
    V_op = CumulativeIntegral(S; dims=3)
    v = compute!(Field(V_op))
    fill_halo_regions!(v)
    
    # Random secondary circulation
    u(x, y, z) = 1e-14 * randn()
    w(x, y, z) = 1e-14 * randn()
    
    return (; u, v, w, b)
end
