@inline function DvDt_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    u = fields.u
    v = fields.v
    w = fields.w

    v_dfm = dependency_fields.v_dfm
    v_next = fields.v_next
    
    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    ∂v∂t = @inbounds (v_next[i, j, k] - v[i, j, k]) / clock.last_Δt
    
    return ∂v∂t + div_𝐯v(i, j, k, grid, weno, total_velocities, v_dfm)
end

@inline function turbulence_h_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    u = fields.u
    v = fields.v
    w = fields.w

    v_dfm = dependency_fields.v_dfm

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    Fv_x = ∂Uv∂x_func(i, j, k, grid, weno, total_velocities.u, v) - ∂Uv∂x_func(i, j, k, grid, weno, total_velocities.u, v_dfm)
    Fv_y = ∂Vv∂y_func(i, j, k, grid, weno, total_velocities.v, v) - ∂Vv∂y_func(i, j, k, grid, weno, total_velocities.v, v_dfm)
    
    return -(Fv_x + Fv_y)
end

@inline function turbulence_z_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    u = fields.u
    v = fields.v
    w = fields.w

    v_dfm = dependency_fields.v_dfm

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    Fv_z = ∂Wv∂z_func(i, j, k, grid, weno, total_velocities.w, v) - ∂Wv∂z_func(i, j, k, grid, weno, total_velocities.w, v_dfm)
    
    return -Fv_z
end
