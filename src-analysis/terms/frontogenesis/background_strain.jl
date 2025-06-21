@inline function background_strain_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    α = variable_strain_rate(clock.time)
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, fields.b)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, fields.b)
    return α * (bx^2 - by^2)
end
strain_dependencies = ()
