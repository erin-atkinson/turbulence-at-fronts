using Oceananigans.Operators

@inline dfm(a) = Field(Average(a; dims=2))

δ = dfm(∂x(u))
ζ = dfm(∂x(v))

ω = dfm(∂z(u) - ∂x(w))


@inline function ω_slice_func(i, j, k, grid, u, w)
    
    du = ℑxzᶜᵃᶜ(i, j, 54, grid, ∂zᶠᶜᶠ, u)
    dw = ℑxzᶜᵃᶜ(i, j, 54, grid, ∂xᶠᶜᶠ, w)
    
    return du .- dw
end

ω_slice_op = KernelFunctionOperation{Center, Center, Nothing}(ω_slice_func, grid, u, w)
ω_slice = Field(ω_slice_op)

outputs = (; δ, ζ, ω, ω_slice)

function update_outputs!(outputs)
    map(compute!, outputs)
end