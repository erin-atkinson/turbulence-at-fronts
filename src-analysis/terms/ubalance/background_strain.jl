@inline function background_strain_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    u_dfm = dependency_fields.u_dfm
    α = variable_strain_rate(clock.time)
    
    return @inbounds α * u_dfm[i, j, k]
end

background_strain_dependencies = (:u_dfm, )
