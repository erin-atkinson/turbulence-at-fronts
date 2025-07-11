@inline function DSP3D_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    α = variable_strain_rate(clock.time)

    u = fields.u
    v = fields.v
    
    u_dfm = dependency_fields.u_dfm
    v_dfm = dependency_fields.v_dfm

    U = fields.U
    V = fields.V

    total_u = SumOfArrays{2}(U, u)
    total_v = SumOfArrays{2}(V, v)
    
    uu = advective_momentum_flux_density_Uu(i, j, k, grid, centered, total_u, u)
    uu_dfm = advective_momentum_flux_density_Uu(i, j, k, grid, centered, total_u, u_dfm)

    vv = advective_momentum_flux_density_Vv(i, j, k, grid, centered, total_v, v)
    vv_dfm = advective_momentum_flux_density_Vv(i, j, k, grid, centered, total_v, v_dfm)
    
    return -α * ((uu - uu_dfm) - (vv - vv_dfm))
end

DSP_dependencies = (
    :u_dfm,
    :v_dfm,
)
