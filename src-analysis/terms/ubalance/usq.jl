@inline function u²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    u_dfm = dependency_fields.u_dfm
    
    return @inbounds u_dfm[i, j, k] * u_dfm[i, j, k]
end

u²_dependencies = (:u_dfm, )
