# Sub-grid diffusion acting on the front strength

@inline function subgrid_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u_dfm = dependency_fields.u_dfm
    
    u = fields.u
    v = fields.v
    w = fields.w

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    Fu = div_effective_Uu(i, j, k, grid, total_velocities, u)
    
    return @inbounds -Fu * u_dfm[i, j, k]
end

subgrid_dependencies = (:u_dfm, )
