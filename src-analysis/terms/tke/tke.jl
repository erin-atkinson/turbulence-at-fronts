@inline function TKE3D_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    
    u_dfm = dependency_fields.u_dfm
    v_dfm = dependency_fields.v_dfm
    w_dfm = dependency_fields.w_dfm
    
    uu = ℑxᶜᵃᵃ(i, j, k, grid, fg, u, u) - ℑxᶜᵃᵃ(i, j, k, grid, fg, u_dfm, u_dfm)
    vv = ℑyᵃᶜᵃ(i, j, k, grid, fg, v, v) - ℑyᵃᶜᵃ(i, j, k, grid, fg, v_dfm, v_dfm)
    ww = ℑzᵃᵃᶜ(i, j, k, grid, fg, w, w) - ℑzᵃᵃᶜ(i, j, k, grid, fg, w_dfm, w_dfm)
    
    return 0.5 * (uu + vv + ww)
end

TKE_dependencies = (
    :u_dfm,
    :v_dfm,
    :w_dfm,
)
