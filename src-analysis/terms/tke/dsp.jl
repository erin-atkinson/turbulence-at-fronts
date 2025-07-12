@inline function DSP3D_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    α = variable_strain_rate(clock.time)

    u = fields.u
    v = fields.v
    
    u_next = fields.u_next
    v_next = fields.v_next

    u_dfm = dependency_fields.u_dfm
    v_dfm = dependency_fields.v_dfm
    
    u_next_dfm = dependency_fields.u_next_dfm
    v_next_dfm = dependency_fields.v_next_dfm

    uu = ℑxᶜᵃᵃ(i, j, k, grid, fGg, u, a_avg, u, u_next)
    vv = ℑyᵃᶜᵃ(i, j, k, grid, fGg, v, a_avg, v, v_next)

    uu_dfm = ℑxᶜᵃᵃ(i, j, k, grid, fGg, u_dfm, a_avg, u_dfm, u_next_dfm)
    vv_dfm = ℑyᵃᶜᵃ(i, j, k, grid, fGg, v_dfm, a_avg, v_dfm, v_next_dfm)
    
    return -α * ((uu - uu_dfm) - (vv - vv_dfm))
end

DSP_dependencies = (
    :u_dfm,
    :v_dfm,
    :u_next_dfm,
    :v_next_dfm,
)
