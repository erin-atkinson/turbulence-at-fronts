# Estimation of dissipation of turbulent kinetic energy

@inline function DTKEDt3D_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    Δt = clock.last_Δt

    u = fields.u
    v = fields.v
    w = fields.w

    u_dfm = dependency_fields.u_dfm
    v_dfm = dependency_fields.v_dfm
    w_dfm = dependency_fields.w_dfm

    u_next = fields.u_next
    v_next = fields.v_next
    w_next = fields.w_next

    u_next_dfm = dependency_fields.u_next_dfm
    v_next_dfm = dependency_fields.v_next_dfm
    w_next_dfm = dependency_fields.w_next_dfm

    U = fields.U
    V = fields.V
    W = fields.W

    total_velocities = (; 
        u = SumOfArrays{2}(u, U),
        v = SumOfArrays{2}(v, V),
        w = SumOfArrays{2}(w, W)
    )

    u′div_𝐯u′ = ℑxᶜᵃᵃ(i, j, k, grid, f′_avg_Gg, u, u_next, u_dfm, u_next_dfm, div_𝐯u′, weno, total_velocities, u, u_dfm)
    v′div_𝐯v′ = ℑyᵃᶜᵃ(i, j, k, grid, f′_avg_Gg, v, v_next, v_dfm, v_next_dfm, div_𝐯v′, weno, total_velocities, v, v_dfm)
    w′div_𝐯w′ = ℑzᵃᵃᶜ(i, j, k, grid, f′_avg_Gg, w, w_next, w_dfm, w_next_dfm, div_𝐯w′, weno, total_velocities, w, w_dfm)

    u′du′dt = ℑxᶜᵃᵃ(i, j, k, grid, f′_avg_Gg, u, u_next, u_dfm, u_next_dfm, df′dt, u, u_next, u_dfm, u_next_dfm, Δt)
    v′dv′dt = ℑyᵃᶜᵃ(i, j, k, grid, f′_avg_Gg, v, v_next, v_dfm, v_next_dfm, df′dt, v, v_next, v_dfm, v_next_dfm, Δt)
    w′dw′dt = ℑzᵃᵃᶜ(i, j, k, grid, f′_avg_Gg, w, w_next, w_dfm, w_next_dfm, df′dt, w, w_next, w_dfm, w_next_dfm, Δt)
    
    return (u′du′dt + v′dv′dt + w′dw′dt) + (u′div_𝐯u′ + v′div_𝐯v′ + w′div_𝐯w′)
end

DTKEDt3D_dependencies = (
    :u_dfm,
    :v_dfm,
    :w_dfm,
    :u_next_dfm,
    :v_next_dfm,
    :w_next_dfm,
)

# Residual of the TKE balance
@inline function ε_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    DTKEDt = dependency_fields.DTKEDt

    LSP = dependency_fields.LSP
    VSP = dependency_fields.VSP

    BFLUX = dependency_fields.BFLUX

    DSP = dependency_fields.DSP
    
    return @inbounds DTKEDt[i, j, k] - (LSP[i, j, k] + VSP[i, j, k] + BFLUX[i, j, k] + DSP[i, j, k])
end

ε_dependencies = (
    :DTKEDt,
    :LSP,
    :VSP,
    :BFLUX,
    :DSP
)