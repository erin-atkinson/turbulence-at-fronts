@inline function strain_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    α = (-∂xᶜᶜᶜ(i, j, k, grid, fields.u) + ∂yᶜᶜᶜ(i, j, k, grid, fields.v)) / 2
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, fields.b)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, fields.b)
    return α * (bx^2 - by^2)
end

strain_dependencies = ()
