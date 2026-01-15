@inline function F_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    b = fields.b
    div_UN² = @inbounds dependency_fields.div_UN²[i, j, k]
    ∂u∂zM² = @inbounds dependency_fields.∂u∂zM²[i, j, k]
    ∂w∂zN² = @inbounds dependency_fields.∂w∂zN²[i, j, k]
    adv_2D = div_UN² + ∂u∂zM² + ∂w∂zN²

    u = fields.u
    v = fields.v
    w = fields.w

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; 
        u = SumOfArrays{2}(u, U),
        v = SumOfArrays{2}(v, V),
        w = SumOfArrays{2}(w, W)
    )

    ∂z_div_Ub = ∂zᶜᶜᶠ(i, j, k, grid, div_Uc, weno, total_velocities, b) + sp.α * ∂zᶜᶜᶠ(i, j, k, grid, b)

    F = ∂z_div_Ub - adv_2D
    
    return -F
end
