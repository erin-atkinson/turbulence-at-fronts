include("terms/terms.jl")
using Oceananigans
using Oceananigans.Operators

u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
b_dfm = dfm(input_fields.b)

@inline M²_func(i, j, k, grid, b) = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b)
@inline N²_func(i, j, k, grid, b) = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b)
@inline S²_func(i, j, k, grid, u, v) = ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, u)^2 + ℑyzᵃᶜᶜ(i, j, k, grid, ∂zᶜᶠᶠ, v)^2
@inline ζ_func(i, j, k, grid, v) = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶠᶜ, v)

M² = Field(KernelFunctionOperation{Center, Nothing, Center}(M²_func, grid, b_dfm))
N² = Field(KernelFunctionOperation{Center, Nothing, Center}(N²_func, grid, b_dfm))
S² = Field(KernelFunctionOperation{Center, Nothing, Center}(S²_func, grid, u_dfm, v_dfm))
ζ = Field(KernelFunctionOperation{Center, Nothing, Center}(ζ_func, grid, v_dfm))

dependency_fields = (; u_dfm, v_dfm, b_dfm, N², M², S², ζ)

output_fields = (; N², M², S², ζ)

skip_update = (:w, :b_next, :u_next, :v_next, :w_next, :pNHS, :pNHS_next)
