using Oceanostics
using Oceananigans.Operators
@inline dfm(a) = Field(Average(a; dims=2))

u_dfm = dfm(u)
v_dfm = dfm(v)
w_dfm = dfm(w)
b_dfm = dfm(b)

#q = Field(@at (Center, Center, Center) KernelFunctionOperation{Face, Face, Face}(Oceanostics.FlowDiagnostics.ertel_potential_vorticity_fff, grid, u, v, w, b, 0, 0, sp.f))
#q_dfm = dfm(q)
mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)



#q_vol_func(i, j, k, grid, q) = @inbounds q[i, j, k] <= 0 ? grid.Δxᶜᵃᵃ * grid.Δyᵃᶜᵃ * grid.Δzᵃᵃᶜ[k] : 0
#q_vol_op = KernelFunctionOperation{Center, Center, Center}(q_vol_func, grid, q)

@inline function q_func(i, j, k, grid, v_dfm, b_dfm, f)
    a = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, v_dfm)
    b = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
    
    c = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, v_dfm)
    d = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b_dfm)
    
    return d * (a + f) - b * c
end

@inline fg(i, j, k, grid, f, g) = @inbounds f[i, j, k] * g[i, j, k]
@inline fGg(i, j, k, grid, f, G, g) = @inbounds f[i, j, k] * G(i, j, k, grid, g)

@inline function q_adv_func(i, j, k, grid, u_dfm, w_dfm, q)
    a = ∂xᶜᶜᶜ(i, j, k, grid,  fGg, u_dfm, ℑxᶠᵃᵃ, q)
    b = ∂zᶜᶜᶜ(i, j, k, grid,  fGg, w_dfm, ℑzᵃᵃᶠ, q)
    
    return -(a + b)
end

q_op = KernelFunctionOperation{Center, Nothing, Center}(q_func, grid, v_dfm, b_dfm, sp.f)
q = Field(q_op)

q_adv_op = KernelFunctionOperation{Center, Nothing, Center}(q_adv_func, grid, u_dfm, w_dfm, q)
q_adv = Field(q_adv_op)

outputs = (; q, q_adv)

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, outputs)
end
