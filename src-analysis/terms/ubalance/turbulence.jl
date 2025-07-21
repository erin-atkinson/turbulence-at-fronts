# Turbulence acting on front strength

@inline function turbulence_h_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w

    u_dfm = dependency_fields.u_dfm

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    Fu_x = ∂Uu∂x_func(i, j, k, grid, weno, total_velocities.u, u) - ∂Uu∂x_func(i, j, k, grid, weno, total_velocities.u, u_dfm)
    Fu_y = ∂Vu∂y_func(i, j, k, grid, weno, total_velocities.v, u) - ∂Vu∂y_func(i, j, k, grid, weno, total_velocities.v, u_dfm)
    
    return @inbounds -(Fu_x + Fu_y)
end

@inline function turbulence_z_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w

    u_dfm = dependency_fields.u_dfm

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    Fu_z = ∂Wu∂z_func(i, j, k, grid, weno, total_velocities.w, u) - ∂Wu∂z_func(i, j, k, grid, weno, total_velocities.w, u_dfm)
    
    return @inbounds -Fu_z
end

turbulence_dependencies = (:u_dfm, )
