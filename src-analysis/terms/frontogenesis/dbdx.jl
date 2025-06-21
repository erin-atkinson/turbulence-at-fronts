@inline function ∇b²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, fields.b)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, fields.b)
    return bx^2 + by^2
end

∇b²_dependencies = (; )

# This should be computed BEFORE ∇b²
@inline function ∇bD∇bDt_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    # Get "next" bx, by
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, fields.b)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, fields.b)

    # Difference with "current"
    ∂∇b²∂t = @inbounds (bx^2 + by^2 - ∇b²[i, j, k]) / clock.last_Δt

    total_u = SumOfArrays{2}(fields.U, fields.u)
    total_v = SumOfArrays{2}(fields.V, fields.v)
    total_w = SumOfArrays{2}(fields.W, fields.w)

    total_velocities = (; u=total_u, v=total_v, w=total_w)

    Fb = div_Uc(i, j, k, grid, centered, total_velocities, ∇b²)

    return (∂∇b²∂t - Fb) / 2
end

∇bD∇bDt_dependencies = ()
