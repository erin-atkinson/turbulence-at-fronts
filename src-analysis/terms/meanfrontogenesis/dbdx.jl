@inline function dbdx_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, dependency_fields.b_dfm)
    return bx
end

dbdx_dependencies = (:b_dfm, )
