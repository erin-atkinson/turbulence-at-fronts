function M²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    b = dependency_fields.b_dfm
    return ∂xᶠᶜᶜ(i, j, k, grid, b)
end
