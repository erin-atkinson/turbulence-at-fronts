# Turbulence acting on front strength

@inline function adv_x_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w

    U = fields.U
    V = fields.V
    W = fields.W

    dbdx = dependency_fields.dbdx

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    return -∂Uc∂x_func(i, j, k, grid, weno, total_velocities.u, dbdx)
end

@inline function adv_y_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    
    U = fields.U
    V = fields.V
    W = fields.W

    dbdx = dependency_fields.dbdx

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    return -∂Vc∂y_func(i, j, k, grid, weno, total_velocities.v, dbdx)
end

@inline function adv_z_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    
    U = fields.U
    V = fields.V
    W = fields.W

    dbdx = dependency_fields.dbdx

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    return -∂Wc∂z_func(i, j, k, grid, weno, total_velocities.w, dbdx)
end
turbulence_dependencies = (:dbdx, )
