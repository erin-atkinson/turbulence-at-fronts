@inline function coriolis_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    v_dfm = ℑxyᶠᶜᵃ(i, j, k, grid, dependency_fields.v_dfm)
    
    return @inbounds sp.f * v_dfm
end

coriolis_dependencies = (:v_dfm, )
