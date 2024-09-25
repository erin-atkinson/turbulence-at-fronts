using Oceananigans

@inline function create_closure(xs, ys, zs, base_state, simulation_parameters)
    sp = simulation_parameters
    #(HorizontalScalarDiffusivity(ν=sp.νₕ, κ=sp.νₕ), VerticalScalarDiffusivity(ν=sp.νᵥ, κ=sp.νᵥ))
    AnisotropicMinimumDissipation()
end