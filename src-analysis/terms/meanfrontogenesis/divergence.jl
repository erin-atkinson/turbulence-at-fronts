@inline function divergence_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    δ = ∂xᶜᶜᶜ(i, j, k, grid, dependency_fields.u_dfm)
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, dependency_fields.b_dfm)
    return -δ * bx
end
divergence_dependencies = (:b_dfm, :u_dfm)
