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