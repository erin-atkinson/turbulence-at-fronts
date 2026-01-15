# Lagrangian tendency
function DN²Dt_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    div_UN² = @inbounds dependency_fields.div_UN²[i, j, k]

    b_next = dependency_fields.b_next_dfm
    b = dependency_fields.b_dfm

    ∂N²∂t = (∂zᶜᶜᶠ(i, j, k, grid, b_next) - ∂zᶜᶜᶠ(i, j, k, grid, b)) / clock.last_Δt

    return ∂N²∂t + div_UN²
end

# Advection
function div_UN²_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    N² = dependency_fields.N²

    u = fields.u
    v = fields.v
    w = fields.w

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (;
        u = SumOfArrays{2}(u, U),
        v = SumOfArrays{2}(v, V),
        w = SumOfArrays{2}(w, W)
    )
    
    return @inbounds div_𝐯w(i, j, k, grid, weno, total_velocities, N²) + sp.α * N²[i, j, k]
end
