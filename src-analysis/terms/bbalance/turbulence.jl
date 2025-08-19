# Turbulence acting on front

@inline function turbulence_h_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    b = fields.b

    b_dfm = dependency_fields.b_dfm

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))
    advection = centered

    Fb_x = ∂Uc∂x_func(i, j, k, grid, weno, total_velocities.u, b) - ∂Uc∂x_func(i, j, k, grid, advection, total_velocities.u, b_dfm)
    Fb_y = ∂Vc∂y_func(i, j, k, grid, weno, total_velocities.v, b) - ∂Vc∂y_func(i, j, k, grid, advection, total_velocities.v, b_dfm)
    
    return @inbounds -(Fb_x + Fb_y)
end

@inline function turbulence_z_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    b = fields.b

    b_dfm = dependency_fields.b_dfm

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))
    advection = centered

    Fb_z = ∂Wc∂z_func(i, j, k, grid, weno, total_velocities.w, b) - ∂Wc∂z_func(i, j, k, grid, advection, total_velocities.w, b_dfm)
    
    return @inbounds -Fb_z
end

turbulence_dependencies = (:b_dfm, )
