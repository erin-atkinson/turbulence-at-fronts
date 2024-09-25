using Oceananigans.Operators

@inline dfm(a) = Field(Average(a; dims=2))
@inline strain_function(t, x, α, f, ℓ) = -α * (1 - exp(-f * t / 15)) * (tanh((x-4ℓ) / ℓ) + tanh((-x-4ℓ) / ℓ)) / 2

# Terms that appear in the equation that aren't calculated elsewhere
# Center, Center, Center
u_dfm = dfm(u)
v_dfm = dfm(v)
w_dfm = dfm(w)
b_dfm = dfm(b)
ϕ_dfm = dfm(ϕ)

mean_fields = (; u_dfm, v_dfm, w_dfm, ϕ_dfm)

@inline fg(i, j, k, grid, f, g) = @inbounds f[i, j, k] * g[i, j, k]
@inline fGg(i, j, k, grid, f, G, g) = @inbounds f[i, j, k] * G(i, j, k, grid, g)

@inline f′(i, j, k, grid, f, f_dfm) = @inbounds f[i, j, k] - f_dfm[i, 1, k]
@inline f′g′(i, j, k, grid, f, f_dfm, g, g_dfm) = f′(i, j, k, grid, f, f_dfm) * f′(i, j, k, grid, g, g_dfm)
@inline f′Gg′(i, j, k, grid, f, f_dfm, G, g, g_dfm) = f′(i, j, k, grid, f, f_dfm) * G(i, j, k, grid, f′, g, g_dfm)

@inline function DSP_func(i, j, k, grid, u_dfm, sp, t, strain_function)
    α = @inbounds strain_function(t[1], grid.xᶠᵃᵃ[i, j, k], sp.α, sp.f, sp.ℓ)
    uᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u_dfm)
    return α * uᶜᶜᶜ * uᶜᶜᶜ
end

@inline function BFLUX_func(i, j, k, grid, w_dfm, b_dfm)
    wᶜᶜᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w_dfm)
    return @inbounds wᶜᶜᶜ * b_dfm[i, j, k]
end

@inline function PWORK_x_func(i, j, k, grid, u_dfm, ϕ_dfm)
    return -ℑxᶜᵃᵃ(i, j, k, grid, fGg, u_dfm, ∂xᶠᶜᶜ, ϕ_dfm)
end

@inline function PWORK_z_func(i, j, k, grid, w_dfm, ϕ_dfm)
    return -ℑzᵃᵃᶜ(i, j, k, grid, fGg, w_dfm, ∂zᶜᶜᶠ, ϕ_dfm)
end

@inline function COR_func(i, j, k, grid, u_dfm, v_dfm, f)
    return @inbounds f * ℑxᶜᵃᵃ(i, j, k, grid, u_dfm) * v_dfm[i, j, k]
end

DSP_op = KernelFunctionOperation{Center, Nothing, Center}(DSP_func, grid, u_dfm, sp, t, strain_function)
BFLUX_op = KernelFunctionOperation{Center, Nothing, Center}(BFLUX_func, grid, w_dfm, b_dfm)
PWORK_x_op = KernelFunctionOperation{Center, Nothing, Center}(PWORK_x_func, grid, u_dfm, ϕ_dfm)
PWORK_z_op = KernelFunctionOperation{Center, Nothing, Center}(PWORK_z_func, grid, w_dfm, ϕ_dfm)
COR_op = KernelFunctionOperation{Center, Nothing, Center}(COR_func, grid, u_dfm, v_dfm, sp.f)

DSP = Field(DSP_op)
BFLUX = Field(BFLUX_op)
PWORK_x = Field(PWORK_x_op)
PWORK_z = Field(PWORK_z_op)
COR = Field(COR_op)

@inline function u′u′_func(i, j, k, grid, u, u_dfm)
    -ℑxᶜᵃᵃ(i, j, k, grid, f′g′, u, u_dfm, u, u_dfm)
end

@inline function w′w′_func(i, j, k, grid, w, w_dfm)
    -ℑzᵃᵃᶜ(i, j, k, grid, f′g′, w, w_dfm, w, w_dfm)
end

@inline function u′w′_func(i, j, k, grid, u, u_dfm, w, w_dfm)
    -ℑxᶜᵃᵃ(i, j, k, grid, f′Gg′, u, u_dfm, ℑxzᶠᵃᶜ, w, w_dfm)
end

u′u′_op = KernelFunctionOperation{Center, Center, Center}(u′u′_func, grid, u, u_dfm)
u′w′_op = KernelFunctionOperation{Center, Center, Center}(u′w′_func, grid, u, u_dfm, w, w_dfm)
w′w′_op = KernelFunctionOperation{Center, Center, Center}(w′w′_func, grid, w, w_dfm)

u′u′ = Field(u′u′_op)
u′w′ = Field(u′w′_op)
w′w′ = Field(w′w′_op)

u′u′_dfm = dfm(u′u′)
u′w′_dfm = dfm(u′w′)
w′w′_dfm = dfm(w′w′)

turb_fields = (; u′u′, u′w′, w′w′, u′u′_dfm, u′w′_dfm, w′w′_dfm)

@inline function XX_SP_func(i, j, k, grid, u_dfm, u′u′_dfm)
    -∂xᶜᶜᶜ(i, j, k, grid, u_dfm) * u′u′_dfm[i, j, k]
end

@inline function XZ_SP_func(i, j, k, grid, u_dfm, u′w′_dfm)
    -ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, u_dfm) * u′w′_dfm[i, j, k]
end

@inline function ZX_SP_func(i, j, k, grid, w_dfm, u′w′_dfm)
    -ℑxzᶜᵃᶜ(i, j, k, grid, ∂xᶠᶜᶠ, w_dfm) * u′w′_dfm[i, j, k]
end

@inline function ZZ_SP_func(i, j, k, grid, w_dfm, w′w′_dfm)
    -∂zᶜᶜᶜ(i, j, k, grid, w_dfm) * w′w′_dfm[i, j, k]
end

XX_SP_op = KernelFunctionOperation{Center, Nothing, Center}(XX_SP_func, grid, u_dfm, u′u′_dfm)
XZ_SP_op = KernelFunctionOperation{Center, Nothing, Center}(XZ_SP_func, grid, u_dfm, u′w′_dfm)
ZX_SP_op = KernelFunctionOperation{Center, Nothing, Center}(ZX_SP_func, grid, w_dfm, u′w′_dfm)
ZZ_SP_op = KernelFunctionOperation{Center, Nothing, Center}(ZZ_SP_func, grid, w_dfm, w′w′_dfm)

XX_SP = Field(XX_SP_op)
XZ_SP = Field(XZ_SP_op)
ZX_SP = Field(ZX_SP_op)
ZZ_SP = Field(ZZ_SP_op)

outputs = (; DSP, BFLUX, PWORK_x, PWORK_z, COR, XX_SP, XZ_SP, ZX_SP, ZZ_SP)

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, turb_fields)
    map(compute!, outputs)
end