@inline function pressure_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    p_dfm = dependency_fields.p_dfm
    
    px = ∂xᶠᶜᶜ(i, j, k, grid, p_dfm)
    
    return @inbounds -px
end

pressure_dependencies = (:p_dfm, )
