@inline function q_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    
    ηx = ℑyzᵃᶜᶜ(i, j, k, grid, ηx_func, clock, fields, dependency_fields, sp)
    ηy = ℑxzᶜᵃᶜ(i, j, k, grid, ηy_func, clock, fields, dependency_fields, sp)
    ηz = ℑxyᶜᶜᵃ(i, j, k, grid, ηz_func, clock, fields, dependency_fields, sp)

    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, fields.b)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, fields.b)
    bz = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, fields.b)
    
    return (bx * ηx
          + by * ηy
          + bz * ηz)
end

q_dependencies = ()
