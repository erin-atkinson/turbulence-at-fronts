@inline function ∇b²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, dependency_fields.b_dfm)
    return bx^2
end

∇b²_dependencies = (:b_dfm, )
