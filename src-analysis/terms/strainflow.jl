using Oceananigans.Fields: ZeroField

@inline function variable_strain_rate(t)
    turnon = max(1-exp(-sp.f * t / 15), 0)
    return sp.α * turnon
end

# Background velocity fields
U = Field{Face, Nothing, Nothing}(grid)
V = Field{Nothing, Face, Nothing}(grid)
W = ZeroField()

# This is probably unnecessary
const Xs = Field{Face, Nothing, Nothing}(grid)
const Ys = Field{Nothing, Face, Nothing}(grid)

Xs.data.parent[:, 1, 1] .= grid.xᶠᵃᵃ.parent
Ys.data.parent[1, :, 1] .= grid.yᵃᶠᵃ.parent

# U and V are calculated directly to avoid issues with boundary conditions
# Open does NOT work
@inline function compute_background!(U, V, W, clock)
    α = variable_strain_rate(clock.time)
    
    # Need to bypass periodic halo
    U.data.parent .= -α * Xs.data.parent
    V.data.parent .= 0 # α * Ys.data.parent
    
    return nothing
end
