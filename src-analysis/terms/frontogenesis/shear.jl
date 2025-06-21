@inline function shear_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    uy = ℑxyᶜᶜᵃ(i, j, k, grid, ∂yᶠᶠᶜ, fields.u)
    vx = ℑxyᶜᶜᵃ(i, j, k, grid, ∂xᶠᶠᶜ, fields.v)
    τ = (uy + vx) / 2
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, fields.b)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, fields.b)
    return -2τ * bx * by
end

shear_dependencies = ()
