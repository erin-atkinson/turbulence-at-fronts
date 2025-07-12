@inline function VSP3D_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w

    u_dfm = dependency_fields.u_dfm
    v_dfm = dependency_fields.v_dfm
    w_dfm = dependency_fields.w_dfm

    u_next_dfm = dependency_fields.u_next_dfm
    v_next_dfm = dependency_fields.v_next_dfm
    w_next_dfm = dependency_fields.w_next_dfm

    W = fields.W
    
    total_w = SumOfArrays{2}(w, W)
    
    wu = advective_momentum_flux_density_Wu(i, j, k, grid, weno, total_w, u)
    wv = advective_momentum_flux_density_Wv(i, j, k, grid, weno, total_w, v)
    ww = advective_momentum_flux_density_Ww(i, j, k, grid, weno, total_w, w)

    wu_dfm = advective_momentum_flux_density_Wu(i, j, k, grid, weno, total_w, u_dfm)
    wv_dfm = advective_momentum_flux_density_Wv(i, j, k, grid, weno, total_w, v_dfm)
    ww_dfm = advective_momentum_flux_density_Ww(i, j, k, grid, weno, total_w, w_dfm)
    
    uz = ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, a_avg, u_dfm, u_next_dfm)
    vz = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, a_avg, v_dfm, v_next_dfm)
    wz = ∂zᶜᶜᶜ(i, j, k, grid, a_avg, w_dfm, w_next_dfm)
    
    return -(
          (wu - wu_dfm) * uz
        + (wv - wv_dfm) * vz
        + (ww - ww_dfm) * wz
    )
end

VSP_dependencies = (
    :u_dfm,
    :v_dfm,
    :w_dfm,
    :u_next_dfm,
    :v_next_dfm,
    :w_next_dfm,
)
