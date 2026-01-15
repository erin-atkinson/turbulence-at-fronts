function N²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    b = dependency_fields.b_dfm
    return ∂zᶜᶜᶠ(i, j, k, grid, b)
end
