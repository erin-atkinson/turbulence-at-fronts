@inline function F_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    b = fields.b
    div_UM² = @inbounds dependency_fields.div_UM²[i, j, k]
    ∂u∂xM² = @inbounds dependency_fields.∂u∂xM²[i, j, k]
    ∂w∂xN² = @inbounds dependency_fields.∂w∂xN²[i, j, k]
    adv_2D = div_UM² + ∂u∂xM² + ∂w∂xN²

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

    ∂x_div_Ub = ∂xᶠᶜᶜ(i, j, k, grid, div_Uc, weno, total_velocities, b) + sp.α * ∂xᶠᶜᶜ(i, j, k, grid, b)

    F = ∂x_div_Ub - adv_2D
    
    return -F
end

@inline function Fx_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    b = fields.b
    u∂M²∂x = @inbounds dependency_fields.u∂M²∂x[i, j, k]
    ∂u∂xM² = @inbounds dependency_fields.∂u∂xM²[i, j, k]
    adv_2D = u∂M²∂x + ∂u∂xM²

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
    
    ∂x_u∂b∂x = ∂xᶠᶜᶜ(i, j, k, grid, ∂Uc∂x_func, weno, total_velocities.u, b) + sp.α * ∂xᶠᶜᶜ(i, j, k, grid, b)

    Fx = ∂x_u∂b∂x - adv_2D
    
    return -Fx
end

@inline function Fz_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    b = fields.b
    w∂M²∂z = @inbounds dependency_fields.w∂M²∂z[i, j, k]
    ∂w∂xN² = @inbounds dependency_fields.∂w∂xN²[i, j, k]
    adv_2D = w∂M²∂z + ∂w∂xN²

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
    
    ∂x_w∂b∂z = ∂xᶠᶜᶜ(i, j, k, grid, ∂Wc∂z_func, weno, total_velocities.w, b)

    Fz = ∂x_w∂b∂z - adv_2D
    
    return -Fz
end
