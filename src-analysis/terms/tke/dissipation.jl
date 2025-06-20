# Dissipation of turbulent kinetic energy

@inline function ε3D_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w

    u_dfm = dependency_fields.u_dfm
    v_dfm = dependency_fields.v_dfm
    w_dfm = dependency_fields.w_dfm

    U = dependency_fields.U
    V = dependency_fields.V
    W = dependency_fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    Fu = div_effective_Uu(i, j, k, grid, total_velocities, u)
    Fv = div_effective_Uv(i, j, k, grid, total_velocities, v)
    Fw = div_effective_Uw(i, j, k, grid, total_velocities, w)
    
    return -(
          (u - u_dfm) * Fu
        + (v - v_dfm) * Fv
        + (w - w_dfm) * Fw
    )
end

ε_dependencies = (
    :u_dfm,
    :v_dfm,
    :w_dfm
)
