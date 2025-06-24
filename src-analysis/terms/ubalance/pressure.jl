@inline function pressure_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    p_dfm = dependency_fields.p_dfm
    u_dfm = dependency_fields.u_dfm
    
    px = ∂xᶠᶜᶜ(i, j, k, grid, p_dfm)
    
    return @inbounds -px * u_dfm[i, j, k]
end

pressure_dependencies = (:u_dfm, :p_dfm)
