using Oceananigans
using Oceananigans.Operators
using Oceananigans.Advection
using Oceananigans.Advection: advective_momentum_flux_Uu,
                              advective_momentum_flux_Vu,
                              advective_momentum_flux_Wu,
                              advective_momentum_flux_Uv,
                              advective_momentum_flux_Vv,
                              advective_momentum_flux_Wv,
                              advective_momentum_flux_Uw,
                              advective_momentum_flux_Vw,
                              advective_momentum_flux_Ww

# These are just for a reminder
# div_𝐯u, div_𝐯v, div_𝐯w, div_Uc,
# advective_tracer_flux_x,
# advective_tracer_flux_y,
# advective_tracer_flux_z,

@inline advective_momentum_flux_density_Uu(i, j, k, grid, advection, U, u) = @inbounds advective_momentum_flux_Uu(i, j, k, grid, advection, U, u) / Axᶜᶜᶜ(i, j, k, grid)
@inline advective_momentum_flux_density_Vu(i, j, k, grid, advection, V, u) = @inbounds advective_momentum_flux_Vu(i, j, k, grid, advection, V, u) / Ayᶠᶠᶜ(i, j, k, grid)
@inline advective_momentum_flux_density_Wu(i, j, k, grid, advection, W, u) = @inbounds advective_momentum_flux_Wu(i, j, k, grid, advection, W, u) / Azᶠᶜᶠ(i, j, k, grid)

@inline advective_momentum_flux_density_Uv(i, j, k, grid, advection, U, v) = @inbounds advective_momentum_flux_Uv(i, j, k, grid, advection, U, v) / Axᶠᶠᶜ(i, j, k, grid)
@inline advective_momentum_flux_density_Vv(i, j, k, grid, advection, V, v) = @inbounds advective_momentum_flux_Vv(i, j, k, grid, advection, V, v) / Ayᶜᶜᶜ(i, j, k, grid)
@inline advective_momentum_flux_density_Wv(i, j, k, grid, advection, W, v) = @inbounds advective_momentum_flux_Wv(i, j, k, grid, advection, W, v) / Azᶜᶠᶠ(i, j, k, grid)

@inline advective_momentum_flux_density_Uw(i, j, k, grid, advection, U, w) = @inbounds advective_momentum_flux_Uw(i, j, k, grid, advection, U, w) / Axᶠᶜᶠ(i, j, k, grid)
@inline advective_momentum_flux_density_Vw(i, j, k, grid, advection, V, w) = @inbounds advective_momentum_flux_Vw(i, j, k, grid, advection, V, w) / Ayᶜᶠᶠ(i, j, k, grid)
@inline advective_momentum_flux_density_Ww(i, j, k, grid, advection, W, w) = @inbounds advective_momentum_flux_Ww(i, j, k, grid, advection, W, w) / Azᶜᶜᶜ(i, j, k, grid)

# Derivatives of fluxes
@inline ∂Uu∂x_func(i, j, k, grid, advection, U, u) = @inbounds δxᶠᵃᵃ(i, j, k, grid, advective_momentum_flux_Uu, advection, U, u) / Vᶠᶜᶜ(i, j, k, grid)
@inline ∂Vu∂y_func(i, j, k, grid, advection, V, u) = @inbounds δyᵃᶜᵃ(i, j, k, grid, advective_momentum_flux_Vu, advection, V, u) / Vᶠᶜᶜ(i, j, k, grid)
@inline ∂Wu∂z_func(i, j, k, grid, advection, W, u) = @inbounds δzᵃᵃᶜ(i, j, k, grid, advective_momentum_flux_Wu, advection, W, u) / Vᶠᶜᶜ(i, j, k, grid)

@inline ∂Uv∂x_func(i, j, k, grid, advection, U, v) = @inbounds δxᶜᵃᵃ(i, j, k, grid, advective_momentum_flux_Uv, advection, U, v) / Vᶜᶠᶜ(i, j, k, grid)
@inline ∂Vv∂y_func(i, j, k, grid, advection, V, v) = @inbounds δyᵃᶠᵃ(i, j, k, grid, advective_momentum_flux_Vv, advection, V, v) / Vᶜᶠᶜ(i, j, k, grid)
@inline ∂Wv∂z_func(i, j, k, grid, advection, W, v) = @inbounds δzᵃᵃᶜ(i, j, k, grid, advective_momentum_flux_Wv, advection, W, v) / Vᶜᶠᶜ(i, j, k, grid)

@inline ∂Uw∂x_func(i, j, k, grid, advection, U, w) = @inbounds δxᶜᵃᵃ(i, j, k, grid, advective_momentum_flux_Uw, advection, U, w) / Vᶜᶜᶠ(i, j, k, grid)
@inline ∂Vw∂y_func(i, j, k, grid, advection, V, w) = @inbounds δyᵃᶜᵃ(i, j, k, grid, advective_momentum_flux_Vw, advection, V, w) / Vᶜᶜᶠ(i, j, k, grid)
@inline ∂Ww∂z_func(i, j, k, grid, advection, W, w) = @inbounds δzᵃᵃᶠ(i, j, k, grid, advective_momentum_flux_Ww, advection, W, w) / Vᶜᶜᶠ(i, j, k, grid)

# Tracers
@inline advective_tracer_flux_density_x(i, j, k, grid, advection, U, c) = @inbounds advective_tracer_flux_x(i, j, k, grid, advection, U, c) / Axᶠᶜᶜ(i, j, k, grid)
@inline advective_tracer_flux_density_y(i, j, k, grid, advection, V, c) = @inbounds advective_tracer_flux_y(i, j, k, grid, advection, V, c) / Ayᶜᶠᶜ(i, j, k, grid)
@inline advective_tracer_flux_density_z(i, j, k, grid, advection, W, c) = @inbounds advective_tracer_flux_z(i, j, k, grid, advection, W, c) / Azᶜᶜᶠ(i, j, k, grid)

@inline ∂Uc∂x_func(i, j, k, grid, advection, U, c) = @inbounds δxᶜᵃᵃ(i, j, k, grid, advective_tracer_flux_x, advection, U, c) / Vᶜᶜᶜ(i, j, k, grid)
@inline ∂Vc∂y_func(i, j, k, grid, advection, V, c) = @inbounds δyᵃᶜᵃ(i, j, k, grid, advective_tracer_flux_y, advection, V, c) / Vᶜᶜᶜ(i, j, k, grid)
@inline ∂Wc∂z_func(i, j, k, grid, advection, W, c) = @inbounds δzᵃᵃᶜ(i, j, k, grid, advective_tracer_flux_z, advection, W, c) / Vᶜᶜᶜ(i, j, k, grid)

@inline div_𝐯u′(i, j, k, grid, advection, total_velocities, u, u_dfm) = div_𝐯u(i, j, k, grid, advection, total_velocities, u) - div_𝐯u(i, j, k, grid, advection, total_velocities, u_dfm)
@inline div_𝐯v′(i, j, k, grid, advection, total_velocities, v, v_dfm) = div_𝐯v(i, j, k, grid, advection, total_velocities, v) - div_𝐯v(i, j, k, grid, advection, total_velocities, v_dfm)
@inline div_𝐯w′(i, j, k, grid, advection, total_velocities, w, w_dfm) = div_𝐯w(i, j, k, grid, advection, total_velocities, w) - div_𝐯w(i, j, k, grid, advection, total_velocities, w_dfm)
