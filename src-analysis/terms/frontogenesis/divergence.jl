@inline function divergence_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    δ = ∂xᶜᶜᶜ(i, j, k, grid, fields.u) + ∂yᶜᶜᶜ(i, j, k, grid, fields.v)
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, fields.b)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, fields.b)
    return -δ * (bx^2 + by^2) / 2
end
divergence_dependencies = ()
