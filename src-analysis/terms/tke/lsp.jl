@inline function LSP3D_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w

    u_dfm = dependency_fields.u_dfm
    v_dfm = dependency_fields.v_dfm
    w_dfm = dependency_fields.w_dfm

    U = fields.U

    total_u = SumOfArrays{2}(u, U)
    
    uu = advective_momentum_flux_density_Uu(i, j, k, grid, weno, total_u, u)
    uv = advective_momentum_flux_density_Uv(i, j, k, grid, weno, total_u, v)
    uw = advective_momentum_flux_density_Uw(i, j, k, grid, weno, total_u, w)

    uu_dfm = advective_momentum_flux_density_Uu(i, j, k, grid, weno, total_u, u_dfm)
    uv_dfm = advective_momentum_flux_density_Uv(i, j, k, grid, weno, total_u, v_dfm)
    uw_dfm = advective_momentum_flux_density_Uw(i, j, k, grid, weno, total_u, w_dfm)
    
    ux = ∂xᶜᶜᶜ(i, j, k, grid, u_dfm)
    vx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, v_dfm)
    wx = ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, w_dfm)
    
    return -(
          (uu - uu_dfm) * ux
        + (uv - uv_dfm) * vx
        + (uw - uw_dfm) * wx
    )
end

LSP_dependencies = (
    :u_dfm,
    :v_dfm,
    :w_dfm,
)
