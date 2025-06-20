# Sub-grid diffusion acting on the front strength

@inline function subgrid_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    b = fields.b

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    Fbx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, div_effective_Uc, total_velocities, b)
    Fby = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, div_effective_Uc, total_velocities, b)
    
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, fields.b)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, fields.b)
    
    return -(bx * Fbx + by * Fby)
end

subgrid_dependencies = ()
