# Turbulence acting on front strength

@inline function turbulence_func(i, j, k, grid, clock, fields, dependency_fields, sp)

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

    Fbx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, div_Uc, centered, total_velocities, b) - ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, div_Uc, centered, total_velocities, b_dfm)
    Fby = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, div_Uc, centered, total_velocities, b) - ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, div_Uc, centered, total_velocities, b_dfm)
    
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, b_dfm)
    
    return -(bx * Fbx + by * Fby)
end

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

    Fbx_x = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂Uc∂x_func, centered, total_velocities.u, b) - ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂Uc∂x_func, centered, total_velocities.u, b_dfm)
    Fby_x = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, ∂Uc∂x_func, centered, total_velocities.u, b) - ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, ∂Uc∂x_func, centered, total_velocities.u, b_dfm)

    Fbx_y = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂Vc∂y_func, centered, total_velocities.v, b) - ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂Vc∂y_func, centered, total_velocities.v, b_dfm)
    Fby_y = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, ∂Vc∂y_func, centered, total_velocities.v, b) - ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, ∂Vc∂y_func, centered, total_velocities.v, b_dfm)
    
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, b_dfm)
    
    return -(bx * (Fbx_x + Fbx_y) + by * (Fby_x + Fby_y))
end

@inline function turbulence_z_func(i, j, k, grid, clock, fields, dependency_fields, sp)

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

    Fbx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂Wc∂z_func, centered, total_velocities.w, b) - ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, ∂Wc∂z_func, centered, total_velocities.w, b_dfm)
    Fby = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, ∂Wc∂z_func, centered, total_velocities.w, b) - ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, ∂Wc∂z_func, centered, total_velocities.w, b_dfm)
    
    bx = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
    by = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, b_dfm)
    
    return -(bx * Fbx + by * Fby)
end

turbulence_dependencies = (:u_dfm, :v_dfm, :w_dfm, :b_dfm, )
