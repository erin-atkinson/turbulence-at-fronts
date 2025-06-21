@inline function tilting_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    wx = ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, fields.w)
    wy = ℑyzᵃᶜᶜ(i, j, k, grid, ∂yᶜᶠᶠ, fields.w)
    
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, fields.b)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, fields.b)
    bz = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, fields.b)
    return -(bx * wx * bz + by * wy * bz)
end

tilting_dependencies = ()
