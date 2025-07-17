@inline function uq_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    U = fields.U

    q = dependency_fields.q

    return advective_tracer_flux_density_x(i, j, k, grid, weno, SumOfArrays{2}(u, U), q)
end

uq_dependencies = (:q, )

@inline function wq_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    w = fields.w
    W = fields.W

    q = dependency_fields.q

    return advective_tracer_flux_density_z(i, j, k, grid, weno, SumOfArrays{2}(w, W), q)
end

wq_dependencies = (:q, )

@inline function div_Uq_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    
    u = fields.u
    v = fields.v
    w = fields.w

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    q = dependency_fields.q

    return div_Uc(i, j, k, grid, weno, total_velocities, q)
end

div_Uq_dependencies = (:q, )