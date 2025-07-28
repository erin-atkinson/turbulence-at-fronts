# Turbulence acting on front strength

@inline function turbulence_h_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    b = fields.b

    u_dfm = dependency_fields.u_dfm
    v_dfm = dependency_fields.v_dfm
    w_dfm = dependency_fields.w_dfm
    b_dfm = dependency_fields.b_dfm

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))
    
    mean_velocities = (; u = SumOfArrays{2}(u_dfm, U),
                        v = SumOfArrays{2}(v_dfm, V),
                        w = SumOfArrays{2}(w_dfm, W))

    Fbx_x = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂Uc∂x_func, weno, total_velocities.u, b) - ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂Uc∂x_func, weno, total_velocities.u, b_dfm)
    Fbx_y = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂Vc∂y_func, weno, total_velocities.v, b) - ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂Vc∂y_func, weno, total_velocities.v, b_dfm)
    
    return -(Fbx_x + Fbx_y)
end

@inline function turbulence_z_func(i, j, k, grid, clock, fields, dependency_fields, sp)

    u = fields.u
    v = fields.v
    w = fields.w
    b = fields.b

    b_dfm = dependency_fields.b_dfm

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; u = SumOfArrays{2}(u, U),
                        v = SumOfArrays{2}(v, V),
                        w = SumOfArrays{2}(w, W))

    Fbx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂Wc∂z_func, weno, total_velocities.w, b) - ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂Wc∂z_func, weno, total_velocities.w, b_dfm)
    
    return -Fbx
end

turbulence_dependencies = (:u_dfm, :v_dfm, :w_dfm, :b_dfm, )
