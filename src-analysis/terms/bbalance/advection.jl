@inline function adv_α_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w

    U = fields.U
    V = fields.V
    W = fields.W

    b_dfm = dependency_fields.b_dfm

    total_velocities = (; u = U,
                        v = V,
                        w = W)
    advection = centered
    
    return -(∂Uc∂x_func(i, j, k, grid, advection, total_velocities.u, b_dfm)+∂Vc∂y_func(i, j, k, grid, advection, total_velocities.v, b_dfm))
end


@inline function adv_x_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w

    U = fields.U
    V = fields.V
    W = fields.W

    b_dfm = dependency_fields.b_dfm

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))
    advection = centered
    
    return -∂Uc∂x_func(i, j, k, grid, advection, total_velocities.u, b_dfm)
end

@inline function adv_y_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w

    U = fields.U
    V = fields.V
    W = fields.W

    b_dfm = dependency_fields.b_dfm

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))
    advection = centered

    return -∂Vc∂y_func(i, j, k, grid, advection, total_velocities.v, b_dfm)
end

@inline function adv_z_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w

    U = fields.U
    V = fields.V
    W = fields.W

    b_dfm = dependency_fields.b_dfm

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))
    advection = centered

    return -∂Wc∂z_func(i, j, k, grid, advection, total_velocities.w, b_dfm)
end

advection_dependencies = (:b_dfm, )
