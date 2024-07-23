using Oceananigans
using Oceananigans.Operators

@inline σ₁(x) = abs(x) > 1 ? 0 : (1 - abs(x))^2

@inline function strain_forcing_u(i, j, k, grid, clock, model_fields, p)
    
    ∂xu = ℑxᶠᵃᵃ(i, j, k, grid, ∂xᶜᶜᶜ, model_fields.u)
    ∂yu = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, model_fields.u)
    
    return @inbounds p.α * (grid.xᶠᵃᵃ[i, j, k] * ∂xu[i, j, k] - grid.yᵃᶜᵃ[i, j, k] * ∂yu[i, j, k] + model_fields.u[i, j, k])
end

@inline function strain_forcing_v(i, j, k, grid, clock, model_fields, p)
    
    ∂xv = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶠᶜ, model_fields.v)
    ∂yv = ℑyᵃᶠᵃ(i, j, k, grid, ∂yᶜᶜᶜ, model_fields.v)
    
    return @inbounds p.α * (grid.xᶜᵃᵃ[i, j, k] * ∂xv[i, j, k] - grid.yᵃᶠᵃ[i, j, k] * ∂yv[i, j, k] - model_fields.v[i, j, k])
end

@inline function strain_forcing_w(i, j, k, grid, clock, model_fields, p)
    
    ∂xw = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶠ, model_fields.w)
    ∂yw = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶠ, model_fields.w)
    
    return @inbounds p.α * (grid.xᶜᵃᵃ[i, j, k] * ∂xw[i, j, k] - grid.yᵃᶜᵃ[i, j, k] * ∂yw[i, j, k])
end

@inline function strain_forcing_b(i, j, k, grid, clock, model_fields, p)
    
    ∂xb = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, model_fields.b)
    ∂yb = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, model_fields.b)
    
    return @inbounds p.α * (grid.xᶜᵃᵃ[i, j, k] * ∂xb[i, j, k] - grid.yᵃᶜᵃ[i, j, k] * ∂yb[i, j, k])
end

@inline function create_forcings(base_state, simulation_parameters; verbose=true)
    α = simulation_parameters.α
    Lx = simulation_parameters.L
    Ly = simulation_parameters.L * (simulation_parameters.Ny / simulation_parameters.Nx) 
    
    damping_width = simulation_parameters.damping_frac * simulation_parameters.a
    rate = simulation_parameters.damping_rate
    
    @inline damping_shape_x(x) = (1 - σ₁((x-Lx/2) / damping_width)) * (1 - σ₁((x+Lx/2) / damping_width))
    @inline damping_shape_y(y) = (1 - σ₁((y-Ly/2) / damping_width)) * (1 - σ₁((y+Ly/2) / damping_width)) 
    @inline mask(x, y, z) = abs((1 - damping_shape_x(x) * damping_shape_y(y)))
    
    u_forcing = (Forcing(strain_forcing_u; discrete_form=true, parameters=(; α)), Relaxation(; rate, mask, target=(x, y, z, t)->base_state.u(x, y, z)))
    v_forcing = (Forcing(strain_forcing_v; discrete_form=true, parameters=(; α)), Relaxation(; rate, mask, target=(x, y, z, t)->base_state.v(x, y, z)))
    w_forcing = (Forcing(strain_forcing_w; discrete_form=true, parameters=(; α)), Relaxation(; rate, mask, target=(x, y, z, t)->base_state.w(x, y, z)))
    b_forcing = (Forcing(strain_forcing_b; discrete_form=true, parameters=(; α)), Relaxation(; rate, mask, target=(x, y, z, t)->base_state.b(x, y, z)))
    
    
    
    @info "Created forcing functions"
    (; u=u_forcing, v=v_forcing, w=w_forcing, b=b_forcing)
end