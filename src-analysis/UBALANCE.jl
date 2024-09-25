using Oceananigans.Operators

@inline dfm(a) = Field(Average(a; dims=2))
@inline strain_function(t, x, α, f, ℓ) = -α * (1 - exp(-f * t / 15)) * (tanh((x-4ℓ) / ℓ) + tanh((-x-4ℓ) / ℓ)) / 2

# Terms that appear in the equation that aren't calculated elsewhere
# Face, Center, Center
u_dfm = dfm(u)
v_dfm = dfm(v)
w_dfm = dfm(w)
φ_dfm = dfm(φ)

mean_fields = (; u_dfm, v_dfm, w_dfm, φ_dfm)

@inline fg(i, j, k, grid, f, g) = @inbounds f[i, j, k] * g[i, j, k]
@inline fGg(i, j, k, grid, f, G, g) = @inbounds f[i, j, k] * G(i, j, k, grid, g)

@inline function u_dfm_u_dfm_x_func(i, j, k, grid, u_dfm)
    return -ℑxᶠᵃᵃ(i, j, k, grid, ∂xᶜᶜᶜ, fg, u_dfm, u_dfm)
end

@inline function u_dfm_w_dfm_z_func(i, j, k, grid, u_dfm, w_dfm)
    return -ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, fGg, u_dfm, ℑxzᶠᵃᶜ, w_dfm)
end

@inline function φ_dfm_x_func(i, j, k, grid, φ_dfm)
    return -∂xᶜᶜᶜ(i, j, k, grid, φ_dfm)
end

@inline function fv_dfm_func(i, j, k, grid, v_dfm, sp)
    return sp.f * ℑxyᶠᶜᵃ(i, j, k, grid, v_dfm)
end

@inline function strain_func(i, j, k, grid, u_dfm, sp, t, strain_function)
    α = @inbounds strain_function(t, grid.xᶠᵃᵃ[i, j, k], sp.α, sp.f, sp.ℓ)
    return @inbounds α * grid.xᶠᵃᵃ[i, j, k] * ℑxᶠᵃᵃ(i, j, k, grid, ∂xᶜᶜᶜ, u_dfm) + 2α * u_dfm[i, j, k]
end

@inline f′(i, j, k, grid, f, f_dfm) = @inbounds f[i, j, k] - f_dfm[i, 1, k]
@inline f′g′(i, j, k, grid, f, f_dfm, g, g_dfm) = f′(i, j, k, grid, f, f_dfm) * f′(i, j, k, grid, g, g_dfm)
@inline f′Gg′(i, j, k, grid, f, f_dfm, G, g, g_dfm) = f′(i, j, k, grid, f, f_dfm) * G(i, j, k, grid, f′, g, g_dfm)

@inline function u′u′_x_func(i, j, k, grid, u, u_dfm)
    -ℑxᶠᵃᵃ(i, j, k, grid, ∂xᶜᶜᶜ, f′g′, u, u_dfm, u, u_dfm)
end

@inline function u′w′_z_func(i, j, k, grid, u, u_dfm, w, w_dfm)
    -ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, f′Gg′, u, u_dfm, ℑxzᶠᵃᶜ, w, w_dfm)
end

u_dfm_u_dfm_x_op = KernelFunctionOperation{Face, Nothing, Center}(u_dfm_u_dfm_x_func, grid, u_dfm)
u_dfm_w_dfm_z_op = KernelFunctionOperation{Face, Nothing, Center}(u_dfm_w_dfm_z_func, grid, u_dfm, w_dfm)
φ_dfm_x_op = KernelFunctionOperation{Face, Nothing, Center}(φ_dfm_x_func, grid, φ_dfm)
fv_dfm_op = KernelFunctionOperation{Face, Nothing, Center}(fv_dfm_func, grid, v_dfm, sp)
strain_op = KernelFunctionOperation{Face, Nothing, Center}(strain_func, grid, u_dfm, sp, t, strain_function)
u′u′_x_op = KernelFunctionOperation{Face, Center, Center}(u′u′_x_func, grid, u, u_dfm)
u′w′_z_op = KernelFunctionOperation{Face, Center, Center}(u′w′_z_func, grid, u, u_dfm, w, w_dfm)

u_dfm_u_dfm_x = Field(u_dfm_u_dfm_x_op)
u_dfm_w_dfm_z = Field(u_dfm_w_dfm_z_op)
φ_dfm_x = Field(φ_dfm_x_op)
fv_dfm = Field(fv_dfm_op)
strain = Field(strain_op)
u′u′_x = dfm(u′u′_x_op)
u′w′_z = dfm(u′w′_z_op)

outputs = (; u_dfm_u_dfm_x, u_dfm_w_dfm_z, φ_dfm_x, fv_dfm, strain, u′u′_x, u′w′_z)

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, outputs)
end