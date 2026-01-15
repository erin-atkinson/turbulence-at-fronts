# These are calculated with linear interpolation

# Convergence
function ∂u∂xM²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    M² = @inbounds dependency_fields.M²[i, j, k]
    u = dependency_fields.u_dfm

    ∂u∂x = ℑxᶠᵃᵃ(i, j, k, grid, ∂xᶜᶜᶜ, u)

    return ∂u∂x * M²
end

# Tilting
function ∂w∂xN²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    b = dependency_fields.b_dfm
    w = dependency_fields.w_dfm

    N² = ℑxzᶠᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b)
    ∂w∂x = ℑzᵃᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, w)

    return ∂w∂x * N²
end
