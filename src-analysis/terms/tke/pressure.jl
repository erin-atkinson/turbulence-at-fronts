@inline function pressure_work(i, j, k, grid, clock, fields, dependency_fields, sp)
    u = fields.u
    w = fields.w
    p = fields.p
    
    u_dfm = dependency_fields.u_dfm
    w_dfm = dependency_fields.w_dfm
    p_dfm = dependency_fields.p_dfm

    upx = ℑxᶜᵃᵃ(i, j, k, grid, f′_avg_Gg, u, u_next, u_dfm, u_next_dfm, ∂xᶠᶜᶜ, f′, p, p_dfm)
    wpz = ℑzᵃᵃᶜ(i, j, k, grid, f′_avg_Gg, w, w_next, w_dfm, w_next_dfm, ∂zᶜᶜᶠ, f′, p, p_dfm)

    return -(upx + wpz) 
end
