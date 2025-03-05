using Oceananigans
using Oceananigans.Operators

# Here we define fields that store the mean buoyancy profiles at the west and east, for forcing

#b_west = Field{Nothing, Nothing, Center}(grid)
#b_east = Field{Nothing, Nothing, Center}(grid)

# Kernel functions that compute the above fields

@inline function b_west_func(i, j, k, grid, b, sp)
    # This is a reduction over i and j
    # mean of b in the left of the outer domain
    j_inds = 1:grid.Ny
    
    total = 0
    
    for _j in j_inds
        total += @inbounds b[1, _j, k]
    end
    
    return total / grid.Ny
end

@inline function b_east_func(i, j, k, grid, b, sp)
    # This is a reduction over i and j
    # mean of b in the right of the outer domain
    j_inds = 1:grid.Ny
    
    total = 0
    
    for _j in j_inds
        total += @inbounds b[grid.Nx, _j, k]
    end
    
    return total / grid.Ny
end

# Damping mask profile
@inline σ₁(x) = abs(x) > 1 ? 0 : (1 - abs(x))^2

@inline damping_shape_x(x) = (1 - σ₁((x-sp.Lx/2) / sp.damping_width)) * (1 - σ₁((x+sp.Lx/2) / sp.damping_width))
@inline velocity_mask(x, y, z) = abs(1 - damping_shape_x(x) * (1 - σ₁((z+sp.Lz) / (sp.Lz-sp.H))))

@inline b_west_mask(x) = σ₁((x+sp.Lx/2) / sp.damping_width)
@inline b_east_mask(x) = σ₁((x-sp.Lx/2) / sp.damping_width)
@inline b_bottom_mask(z) = abs(σ₁((z+sp.Lz) / (sp.Lz-sp.H)))

# Strain turns on slowly
@inline strain_function(t, x, α, f, ℓ) = -α * max(1 - exp(-f * t / 15), 0, 1) * (tanh((x-4ℓ) / ℓ) + tanh((-x-4ℓ) / ℓ)) / 2

@inline function strain_forcing_u(i, j, k, grid, clock, model_fields, p)
    
    strain_rate = @inbounds strain_function(clock.time, grid.xᶠᵃᵃ[i], p.α, p.f, p.ℓ)
    
    ∂xu = ℑxᶠᵃᵃ(i, j, k, grid, ∂xᶜᶜᶜ, model_fields.u)
    ∂yu = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, model_fields.u)
    
    return @inbounds strain_rate * (grid.xᶠᵃᵃ[i, j, k] * ∂xu - grid.yᵃᶜᵃ[i, j, k] * ∂yu + model_fields.u[i, j, k])
end

@inline function strain_forcing_v(i, j, k, grid, clock, model_fields, p)
    
    strain_rate = @inbounds strain_function(clock.time, grid.xᶜᵃᵃ[i], p.α, p.f, p.ℓ)
    
    ∂xv = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶠᶜ, model_fields.v)
    ∂yv = ℑyᵃᶠᵃ(i, j, k, grid, ∂yᶜᶜᶜ, model_fields.v)
    
    return @inbounds strain_rate * (grid.xᶜᵃᵃ[i, j, k] * ∂xv - grid.yᵃᶠᵃ[i, j, k] * ∂yv - model_fields.v[i, j, k])
end

@inline function strain_forcing_w(i, j, k, grid, clock, model_fields, p)
    
    strain_rate = @inbounds strain_function(clock.time, grid.xᶜᵃᵃ[i], p.α, p.f, p.ℓ)
    
    ∂xw = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶠ, model_fields.w)
    ∂yw = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶠ, model_fields.w)
    
    return @inbounds strain_rate * (grid.xᶜᵃᵃ[i, j, k] * ∂xw - grid.yᵃᶜᵃ[i, j, k] * ∂yw)
end

@inline function strain_forcing_b(i, j, k, grid, clock, model_fields, p)
    
    strain_rate = @inbounds strain_function(clock.time, grid.xᶜᵃᵃ[i], p.α, p.f, p.ℓ)
    
    ∂xb = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, model_fields.b)
    ∂yb = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, model_fields.b)
    
    return @inbounds strain_rate * (grid.xᶜᵃᵃ[i, j, k] * ∂xb - grid.yᵃᶜᵃ[i, j, k] * ∂yb)
end

@inline function strain_forcing_ξ(i, j, k, grid, clock, model_fields, p)
    
    strain_rate = @inbounds strain_function(clock.time, grid.xᶜᵃᵃ[i], p.α, p.f, p.ℓ)
    
    ∂xξ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, model_fields.ξ)
    ∂yξ = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, model_fields.ξ)
    
    return @inbounds strain_rate * (grid.xᶜᵃᵃ[i, j, k] * ∂xξ - grid.yᵃᶜᵃ[i, j, k] * ∂yξ)
end

# Velocities get damped to zero at the edge and bottom of the domain
u_forcing = (
    Forcing(strain_forcing_u; discrete_form=true, parameters=(; sp.α, sp.f, sp.ℓ)),
    Relaxation(; rate=sp.damping_rate, mask=velocity_mask, target=(x, y, z, t)->0)
)
v_forcing = (
    Forcing(strain_forcing_v; discrete_form=true, parameters=(; sp.α, sp.f, sp.ℓ)),
    Relaxation(; rate=sp.damping_rate, mask=velocity_mask, target=(x, y, z, t)->0)
)
w_forcing = (
    Forcing(strain_forcing_w; discrete_form=true, parameters=(; sp.α, sp.f, sp.ℓ)),
    Relaxation(; rate=sp.damping_rate, mask=velocity_mask, target=(x, y, z, t)->0)
)

# Buoyancy is damped to the mean profiles, and I need a seperate forcing function

@inline function b_damping_func(i, j, k, grid, clock, model_fields, p)
    # mean of b in the right of the outer domain
    b_east = @inbounds p.b_east[1, 1, k]
    east_rate = @inbounds p.damping_rate * abs(b_east_mask(grid.xᶜᵃᵃ[i]))
    
    b_west = @inbounds p.b_west[1, 1, k]
    west_rate = @inbounds p.damping_rate * abs(b_west_mask(grid.xᶜᵃᵃ[i]))
    
    return @inbounds east_rate * (b_east - model_fields.b[i, j, k]) + west_rate * (b_west - model_fields.b[i, j, k])
end

#b_func = get_base_state_function(sp).b

b_forcing = (
    Forcing(strain_forcing_b; discrete_form=true, parameters=(; sp.α, sp.f, sp.ℓ)),
    Forcing(b_damping_func; discrete_form=true, parameters=(; b_east, b_west, sp.damping_rate)),
    Relaxation(; rate=sp.damping_rate, mask=(x, y, z)->b_bottom_mask(z), target=(x, y, z, t)->base_b_func(x, z, sp))
)

#=
b_forcing = (
    Forcing(strain_forcing_b; discrete_form=true, parameters=(; sp.α, sp.f, sp.ℓ)),
    Relaxation(; rate=sp.damping_rate, mask=(x, y, z)->b_bottom_mask(z), target=(x, y, z, t)->base_b_func(x, z, sp))
)
=#
ξ_forcing = (
    Forcing(strain_forcing_ξ; discrete_form=true, parameters=(; sp.α, sp.f, sp.ℓ)),
)

forcing = (; u=u_forcing, v=v_forcing, w=w_forcing, b=b_forcing, ξ=ξ_forcing)