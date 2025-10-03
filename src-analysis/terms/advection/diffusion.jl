using Oceananigans
using Oceananigans.Operators
using Oceananigans.Advection
using Oceananigans.Utils: SumOfArrays
#= 
=#

# Create advection schemes
const weno = WENO(; order=9)
const centered = Centered(; order=10)


# Momentum
for v in (:u, :v, :w)
    for V in (:U, :V, :W)
        f = Symbol(:effective_diffusive_momentum_flux_, V, v)
        adv_func = Symbol(:advective_momentum_flux_, V, v)
        @eval begin
            @inline function $f(i, j, k, grid, $V, $v)
                return (Oceananigans.Advection.$adv_func(i, j, k, grid, weno, $V, $v) - Oceananigans.Advection.$adv_func(i, j, k, grid, centered, $V, $v))
            end
        end
    end
end

# Divergences, just copy from Oceananigans source
@inline function div_effective_Uu(i, j, k, grid, U, u)
    return 1/Vᶠᶜᶜ(i, j, k, grid) * (δxᶠᵃᵃ(i, j, k, grid, effective_diffusive_momentum_flux_Uu, U[1], u) +
                                    δyᵃᶜᵃ(i, j, k, grid, effective_diffusive_momentum_flux_Vu, U[2], u) +
                                    δzᵃᵃᶜ(i, j, k, grid, effective_diffusive_momentum_flux_Wu, U[3], u))
end

@inline function div_effective_Uv(i, j, k, grid, U, v)
    return 1/Vᶜᶠᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, effective_diffusive_momentum_flux_Uv, U[1], v) +
                                    δyᵃᶠᵃ(i, j, k, grid, effective_diffusive_momentum_flux_Vv, U[2], v) +
                                    δzᵃᵃᶜ(i, j, k, grid, effective_diffusive_momentum_flux_Wv, U[3], v))
end

@inline function div_effective_Uw(i, j, k, grid, U, w)
    return 1/Vᶜᶜᶠ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, effective_diffusive_momentum_flux_Uw, U[1], w) +
                                    δyᵃᶜᵃ(i, j, k, grid, effective_diffusive_momentum_flux_Vw, U[2], w) +
                                    δzᵃᵃᶠ(i, j, k, grid, effective_diffusive_momentum_flux_Ww, U[3], w))
end

# Tracers
for (ξ, v) in zip((:x, :y, :z), (:u, :v, :w))
    f = Symbol(:effective_diffusive_tracer_flux_, ξ)
    adv_func = Symbol(:advective_tracer_flux_, ξ)
    @eval begin
        @inline function $f(i, j, k, grid, $v, c)
            return Oceananigans.Advection.$adv_func(i, j, k, grid, weno, $v, c) - Oceananigans.Advection.$adv_func(i, j, k, grid, centered, $v, c)
        end
    end
end

@inline function div_effective_Uc(i, j, k, grid, U, c)
    return 1/Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, effective_diffusive_tracer_flux_x, U.u, c) +
                                    δyᵃᶜᵃ(i, j, k, grid, effective_diffusive_tracer_flux_y, U.v, c) +
                                    δzᵃᵃᶜ(i, j, k, grid, effective_diffusive_tracer_flux_z, U.w, c))
end
