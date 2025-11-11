# base_state.jl
# Functions describing the initial state of the simulations

using Oceananigans: fill_halo_regions!



@inline f(s) = 1 + tanh(s)
@inline f′(s) = sech(s)^2
const f′′_max = maximum(s->-2tanh(s) * sech(s)^2, range(-10, 10, 1000))

@inline G(s) = log(1 + exp(s))
@inline g(s) = 1 / (1 + exp(-s))
@inline g′(s) = exp(-s) / (1 + exp(-s))^2

# Background stratification
@inline b∞(z, sp) = -sp.λ * sp.H * sp.N₀² * G(-(z + sp.H) / (sp.λ * sp.H))

# Buoyancy
@inline function front_buoyancy(x, z, sp)
    x₁ = x / sp.ℓ + sp.a * (z + sp.H / 2) / sp.H
    z₁ = (z + sp.H) / (sp.λ * sp.H)
    
    return (sp.Δb / 2) * f(x₁) * g(z₁) + b∞(z, sp)
end

# Stratification
@inline function front_N²(x, z, sp)
    x₁ = x / sp.ℓ + sp.a * (z + sp.H / 2) / sp.H
    z₁ = (z + sp.H) / (sp.λ * sp.H)
    
    return (
          (sp.a * sp.Δb / 2sp.H) * f′(x₁) * g(z₁)
        + (sp.Δb / (2sp.λ * sp.H)) * f(x₁) * g′(z₁)
        + sp.N₀² * g(-z₁)
    )
end

# Horizontal buoyancy gradient
@inline function front_M²(x, z, sp)
    x₁ = x / sp.ℓ + sp.a * (z + sp.H / 2) / sp.H
    z₁ = (z + sp.H) / (sp.λ * sp.H)
    
    return (sp.Δb / 2sp.ℓ) * f′(x₁) * g(z₁)
end

# Thermal wind shear
@inline front_S(x, z, sp) = front_M²(x, z, sp) / sp.f

@inline function approximate_front_velocity(x, z, sp)
    x₁ = x / sp.ℓ + sp.a * (z + sp.H / 2) / sp.H
    x₂ = x / sp.ℓ + sp.a * (-sp.H + sp.H / 2) / sp.H
    z₁ = (z + sp.H) / (sp.λ * sp.H)
    
    return (sp.a * sp.H * sp.Δb / (2sp.ℓ * sp.f)) * (f(x₁) - f(x₂)) * g(z₁)
end

@inline front_Ri(x, z, sp) = front_N²(x, z, sp) / front_shear(x, z, sp)^2

@inline function minimum_Ri(sp)
    return 2sp.a * sp.ℓ^2 * sp.f^2 / sp.H / sp.Δb
end

@inline function maximum_front_velocity(sp)
    return 2sp.a^2 * sp.H * sp.Δb / sp.ℓ / sp.f 
end

@inline function maximum_Ro(sp)
    # Maximum Rossby number occurs at the surface
    f′′_max * maximum_front_velocity(sp) / sp.ℓ / sp.f 
end

@inline function create_front_parameters(ip)
    λ = 0.1
    a = (ip.Ro * ip.Ri / A_Ro)^(1/3)'
    ℓ = sqrt(f′′_max * a^2 * ip.H * ip.Δb / 2ip.f^2 / ip.Ro)
    return merge(ip, (; λ, a, ℓ))
end

@inline function front_initial_conditions(grid::RectilinearGrid, sp)
    # Use Oceananigans fields to setup the initial thermal wind properly
    
    b = Field{Center, Center, Center}(grid)
    set!(b, (x, y, z)->front_buoyancy(x, z, sp))
    fill_halo_regions!(b)

    # Compute the thermal wind shear
    S_op = @at (Center, Face, Center) ∂x(b) / sp.f
    S = compute!(Field(S_op))

    # Integrate
    V_op = CumulativeIntegral(S; dims=3)
    v = compute!(Field(V_op))
    fill_halo_regions!(v)
    
    # Random secondary circulation
    u(x, y, z) = 1e-14 * randn()
    w(x, y, z) = 1e-14 * randn()
    
    return (; u, v, w, b)
end
