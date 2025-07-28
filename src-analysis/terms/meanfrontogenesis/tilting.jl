@inline function tilting_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    wx = ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, dependency_fields.w_dfm)
    bz = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, dependency_fields.b_dfm)
    return -wx * bz
end

tilting_dependencies = (:b_dfm, :w_dfm)
