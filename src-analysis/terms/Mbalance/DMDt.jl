# Lagrangian tendency
function DM²Dt_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    u∂M²∂x = @inbounds dependency_fields.u∂M²∂x[i, j, k]
    w∂M²∂z = @inbounds dependency_fields.w∂M²∂z[i, j, k]

    b_next = dependency_fields.b_next_dfm
    b = dependency_fields.b_dfm

    ∂M²∂t = (∂xᶠᶜᶜ(i, j, k, grid, b_next) - ∂xᶠᶜᶜ(i, j, k, grid, b)) / clock.last_Δt

    return ∂M²∂t + u∂M²∂x + w∂M²∂z
end

# Advection
function div_UM²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    M² = dependency_fields.M²

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
    
    return @inbounds div_𝐯u(i, j, k, grid, weno, total_velocities, M²) + sp.α * M²[i, j, k]
end

function u∂M²∂x_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    M² = dependency_fields.M²

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
    
    return @inbounds ∂Uu∂x_func(i, j, k, grid, weno, total_velocities.u, M²) + sp.α * M²[i, j, k]
end

function w∂M²∂z_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    M² = dependency_fields.M²

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
    
    return @inbounds ∂Wu∂z_func(i, j, k, grid, weno, total_velocities.w, M²)
end