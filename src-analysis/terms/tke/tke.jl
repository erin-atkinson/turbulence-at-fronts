@inline function TKE3D_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    
    u_dfm = dependency_fields.u_dfm
    v_dfm = dependency_fields.v_dfm
    w_dfm = dependency_fields.w_dfm
    
    uu = ℑxᶜᵃᵃ(i, j, k, grid, f′g′, u, u_dfm, u, u_dfm)
    vv = ℑyᵃᶜᵃ(i, j, k, grid, f′g′, v, v_dfm, v, v_dfm)
    ww = ℑzᵃᵃᶜ(i, j, k, grid, f′g′, w, w_dfm, w, w_dfm)
    
    return (uu + vv + ww) / 2
end

TKE_dependencies = (
    :u_dfm,
    :v_dfm,
    :w_dfm,
)
