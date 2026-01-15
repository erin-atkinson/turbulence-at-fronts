# These are calculated with linear interpolation

# Convergence
function ∂w∂zN²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    N² = @inbounds dependency_fields.N²[i, j, k]
    w = dependency_fields.w_dfm

    ∂w∂z = ℑzᵃᵃᶠ(i, j, k, grid, ∂zᶜᶜᶜ, w)

    return ∂w∂z * N²
end

# Tilting
function ∂u∂zM²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    b = dependency_fields.b_dfm
    u = dependency_fields.u_dfm

    M² = ℑxzᶜᵃᶠ(i, j, k, grid, ∂xᶠᶜᶜ, b)
    ∂u∂z = ℑxᶜᵃᵃ(i, j, k, grid, ∂zᶠᶜᶠ, u)

    return ∂u∂z * M²
end
