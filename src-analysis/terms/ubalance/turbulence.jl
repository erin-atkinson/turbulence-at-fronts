# Turbulence acting on front strength

@inline function turbulence_h_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    b = fields.b

    u_dfm = dependency_fields.u_dfm
    v_dfm = dependency_fields.v_dfm
    w_dfm = dependency_fields.w_dfm

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))
    
    mean_velocities = (; u = SumOfArrays{2}(u_dfm, U),
                        v = SumOfArrays{2}(v_dfm, V),
                        w = SumOfArrays{2}(w_dfm, W))

    Fu_x = ∂Uu∂x_func(i, j, k, grid, centered, total_velocities.u, u) - ∂Uu∂x_func(i, j, k, grid, centered, total_velocities.u, u_dfm)

    Fu_y = ∂Vu∂y_func(i, j, k, grid, centered, total_velocities.v, u) - ∂Vu∂y_func(i, j, k, grid, centered, total_velocities.v, u_dfm)
    
    return @inbounds -(u_dfm[i, j, k] * (Fu_x + Fu_y))
end

@inline function turbulence_z_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    b = fields.b

    u_dfm = dependency_fields.u_dfm
    v_dfm = dependency_fields.v_dfm
    w_dfm = dependency_fields.w_dfm

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))
    
    mean_velocities = (; u = SumOfArrays{2}(u_dfm, U),
                        v = SumOfArrays{2}(v_dfm, V),
                        w = SumOfArrays{2}(w_dfm, W))

    Fu_z = ∂Wu∂z_func(i, j, k, grid, centered, total_velocities.w, u) - ∂Wu∂z_func(i, j, k, grid, centered, total_velocities.w, u_dfm)
    
    return @inbounds -(u_dfm[i, j, k] * Fu_z)
end

turbulence_dependencies = (:u_dfm, :v_dfm, :w_dfm,)
