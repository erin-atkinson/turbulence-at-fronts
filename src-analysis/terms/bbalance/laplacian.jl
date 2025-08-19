@inline function laplacian_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    b_dfm = dependency_fields.b_dfm
    return ∂xᶜᶜᶜ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
end

laplacian_dependencies = (:b_dfm, )
