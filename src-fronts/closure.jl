using Oceananigans

@inline function create_closure(xs, ys, zs, base_state, simulation_parameters)
    return AnisotropicMinimumDissipation()
    sp = simulation_parameters
    # Suppress frontogenesis below grid scale
    Δx = if sp.cell_ratio == 1
        (xs[2] - xs[1]) / sp.Nx
    else
        minimum(diff(xs))
    end
    Δz = (zs[2] - zs[1]) / sp.Nz
    νₕ = Δx^2  * sp.α / 2
    νᵥ = Δz^2  * sp.α / 2
    νᵥ = νₕ / 10
    #closure = (HorizontalScalarDiffusivity(; ν=sp.ν_mult*νₕ, κ=sp.ν_mult*νₕ / sp.Pr), VerticalScalarDiffusivity(; ν=sp.ν_mult*νᵥ, κ=sp.ν_mult*νᵥ / sp.Pr))
    #closure = (HorizontalScalarDiffusivity(; ν=sp.ν_mult*νₕ, κ=sp.ν_mult*νₕ / sp.Pr), SmagorinskyLilly(; Pr=sp.Pr, Cb=0))
    @info "Created closure"
    return closure
end