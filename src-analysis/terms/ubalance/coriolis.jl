@inline function coriolis_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    u_dfm = dependency_fields.u_dfm
    v_dfm = ℑxyᶠᶜᵃ(i, j, k, grid, dependency_fields.v_dfm)
    
    return @inbounds u_dfm[i, j, k] * sp.f * v_dfm
end

coriolis_dependencies = (:u_dfm, :v_dfm, )
