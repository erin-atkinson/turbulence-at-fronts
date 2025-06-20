@inline function ∇b²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, fields.b)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, fields.b)
    return bx^2 + by^2
end

∇b²_dependencies = ()
