@inline function background_strain_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    α = variable_strain_rate(clock.time)
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, dependency_fields.b_dfm)
    return α * bx
end
strain_dependencies = (:b_dfm, )
