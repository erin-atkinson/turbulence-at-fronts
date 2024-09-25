using Oceananigans
using Oceananigans.Operators

@inline σ₁(x) = abs(x) > 1 ? 0 : (1 - abs(x))^2

@inline function strain_forcing_u(i, j, k, grid, clock, model_fields, p)
    
    ∂xu = ℑxᶠᵃᵃ(i, j, k, grid, ∂xᶜᶜᶜ, model_fields.u)
    ∂yu = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, model_fields.u)
    
    return @inbounds p.strain_function(clock.time, grid.xᶠᵃᵃ[i, j, k], p.α, p.f, p.ℓ) * (grid.xᶠᵃᵃ[i, j, k] * ∂xu - grid.yᵃᶜᵃ[i, j, k] * ∂yu + model_fields.u[i, j, k])
end

@inline function strain_forcing_v(i, j, k, grid, clock, model_fields, p)
    
    ∂xv = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶠᶜ, model_fields.v)
    ∂yv = ℑyᵃᶠᵃ(i, j, k, grid, ∂yᶜᶜᶜ, model_fields.v)
    
    return @inbounds p.strain_function(clock.time, grid.xᶜᵃᵃ[i, j, k], p.α, p.f, p.ℓ) * (grid.xᶜᵃᵃ[i, j, k] * ∂xv - grid.yᵃᶠᵃ[i, j, k] * ∂yv - model_fields.v[i, j, k])
end

@inline function strain_forcing_w(i, j, k, grid, clock, model_fields, p)
    
    ∂xw = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶠ, model_fields.w)
    ∂yw = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶠ, model_fields.w)
    
    return @inbounds p.strain_function(clock.time, grid.xᶜᵃᵃ[i, j, k], p.α, p.f, p.ℓ) * (grid.xᶜᵃᵃ[i, j, k] * ∂xw - grid.yᵃᶜᵃ[i, j, k] * ∂yw)
end

@inline function strain_forcing_b(i, j, k, grid, clock, model_fields, p)
    
    ∂xb = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, model_fields.b)
    ∂yb = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, model_fields.b)
    
    return @inbounds p.strain_function(clock.time, grid.xᶜᵃᵃ[i, j, k], p.α, p.f, p.ℓ) * (grid.xᶜᵃᵃ[i, j, k] * ∂xb - grid.yᵃᶜᵃ[i, j, k] * ∂yb)
end

# Strain turns on slowly
@inline strain_function(t, x, α, f, ℓ) = -α * (1 - exp(-f * t / 15)) * (tanh((x-4ℓ) / ℓ) + tanh((-x-4ℓ) / ℓ)) / 2

@inline function create_forcings(base_state, simulation_parameters)
    sp = simulation_parameters
    Lx = sp.Lx
    Ly = sp.Ly

    damping_width = sp.damping_frac * sp.ℓ
    rate = sp.c * sp.N₀² / (2π)
    
    @inline damping_shape_x(x) = (1 - σ₁((x-sp.Lx/2) / damping_width)) * (1 - σ₁((x+sp.Lx/2) / damping_width))
    @inline mask(x, y, z) = abs( 1 - damping_shape_x(x) * (1 - σ₁((z+sp.Lz) / (sp.Lz-sp.H))) )
    @inline b_mask(x, y, z) = abs(σ₁((z+sp.Lz) / (sp.Lz-sp.H)))
    
    b_func = get_base_state_function(simulation_parameters).b
    
    u_forcing, v_forcing, w_forcing, b_forcing = if sp.α == 0
        u_forcing = (
            Relaxation(; rate, mask, target=(x, y, z, t)->0),
        )
        v_forcing = (
            Relaxation(; rate, mask, target=(x, y, z, t)->0),
        )
        w_forcing = (
            Relaxation(; rate, mask, target=(x, y, z, t)->0),
        )
        b_forcing = (
            Relaxation(; rate, mask=b_mask, target=(x, y, z, t)->b_func(x, z)),
        )
        u_forcing, v_forcing, w_forcing, b_forcing
    else
        u_forcing = (
            Forcing(strain_forcing_u; discrete_form=true, parameters=(; sp.α, sp.f, sp.ℓ, strain_function)),
            Relaxation(; rate, mask, target=(x, y, z, t)->0)
        )
        v_forcing = (
            Forcing(strain_forcing_v; discrete_form=true, parameters=(; sp.α, sp.f, sp.ℓ, strain_function)),
            Relaxation(; rate, mask, target=(x, y, z, t)->0)
        )
        w_forcing = (
            Forcing(strain_forcing_w; discrete_form=true, parameters=(; sp.α, sp.f, sp.ℓ, strain_function)),
            Relaxation(; rate, mask, target=(x, y, z, t)->0)
        )
        b_forcing = (
            Forcing(strain_forcing_b; discrete_form=true, parameters=(; sp.α, sp.f, sp.ℓ, strain_function)),
            Relaxation(; rate, mask=b_mask, target=(x, y, z, t)->b_func(x, z))
        )
        u_forcing, v_forcing, w_forcing, b_forcing
    end
    
    
    

    @info "Created forcing functions"
    (; u=u_forcing, v=v_forcing, w=w_forcing, b=b_forcing)
end