using Oceananigans.Operators
using CUDA: @allowscalar

# ---------------------------------------
# Quadratic damping
@inline function sponge_layer(x, y, z)
    s = min((z+sp.Lz) / (sp.Lz-sp.H), 1)
    return (1 - abs(s))^2
end

# Damp b at the bottom towards a linear profile
@inline function b_forcing_func(i, j, k, grid, clock, model_fields)
    x = @inbounds grid.xᶜᵃᵃ[i]
    y = @inbounds grid.yᵃᶜᵃ[j]
    z = @inbounds grid.zᵃᵃᶜ[k]
    
    b = @inbounds model_fields.b[i, j, k]
    tb = @inbounds model_fields.b[i, j, grid.Hz] + sp.N₀² * (z - grid.zᵃᵃᶜ[grid.Hz])
    
    return sp.σ * (tb - b) * sponge_layer(x, y, z)
end
# ---------------------------------------

# ---------------------------------------
# Strain turns on slowly starting at t=0
@inline function variable_strain_rate(t, α, f)
    turnon = max(1-exp(-f * t / 15), 0)
    return α * turnon
end

# Background velocity fields
U = Field{Face, Nothing, Nothing}(grid)
V = Field{Nothing, Face, Nothing}(grid)

# This is probably unnecessary
Xs = Field{Face, Nothing, Nothing}(grid)
Ys = Field{Nothing, Face, Nothing}(grid)
@allowscalar begin 
    Xs.data.parent[:, 1, 1] .= grid.xᶠᵃᵃ.parent
    Ys.data.parent[1, :, 1] .= grid.yᵃᶠᵃ.parent
end

# U and V are calculated directly to avoid issues with boundary conditions
@inline function calculate_UV_callback(simulation, sp)
    t = simulation.model.clock.time
    α = variable_strain_rate(t, sp.α, sp.f)
    
    # Need to bypass periodic halo
    U.data.parent .= -α * Xs.data.parent
    V.data.parent .= α * Ys.data.parent
    
    return nothing
end
# ---------------------------------------

# ---------------------------------------
# Background velocity
@inline αv_func(x, y, z, t, v) = -variable_strain_rate(t, sp.α, sp.f) * v
@inline αu_func(x, y, z, t, u) = variable_strain_rate(t, sp.α, sp.f) * u
# ---------------------------------------

# ---------------------------------------
# Combination of forcings
u_forcing = (
    AdvectiveForcing(; u=U, v=V),
    Relaxation(; rate=sp.σ, mask=sponge_layer, target=0),
    Forcing(αu_func; field_dependencies=(:u, ))
)
v_forcing = (
    AdvectiveForcing(; u=U, v=V),
    Relaxation(; rate=sp.σ, mask=sponge_layer, target=0),
    Forcing(αv_func; field_dependencies=(:v, )),
)
w_forcing = (
    AdvectiveForcing(; u=U, v=V),
    Relaxation(; rate=sp.σ, mask=sponge_layer, target=0),
)
b_forcing = (
    AdvectiveForcing(; u=U, v=V),
    Forcing(b_forcing_func; discrete_form=true)
)
# ---------------------------------------

forcing = (; u=u_forcing, v=v_forcing, w=w_forcing, b=b_forcing)