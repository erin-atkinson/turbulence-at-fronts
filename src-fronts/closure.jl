closure = AnisotropicMinimumDissipation()
#closure = (VerticalScalarDiffusivity(; ν=sp.νᵥ, κ=sp.νᵥ), HorizontalScalarDiffusivity(; ν=sp.νₕ, κ=sp.νₕ))
#closure = nothing