@inline function Vq_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    return @inbounds dependency_fields.q[i, j, k] < -1e-10 ? 1.0 : 0.0
end

Vq_dependencies = (:q, )
