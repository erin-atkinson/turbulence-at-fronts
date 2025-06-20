include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")

include("terms/frontogenesis/dbdx.jl")
include("terms/frontogenesis/background_strain.jl")
include("terms/frontogenesis/strain.jl")
include("terms/frontogenesis/divergence.jl")
include("terms/frontogenesis/shear.jl")
include("terms/frontogenesis/subgrid.jl")

# Need to avoid computation of arguments of ∇bD∇bDt
const ∇b² = Field(KernelFunctionOperation{Center, Center, Center}(∇b²_func, grid, clock, fields, (; ), sp))
∇bD∇bDt = Field(KernelFunctionOperation{Center, Center, Center}(∇bD∇bDt_func, grid, clock, fields, (; ), sp))

# Resolved terms in frontogenesis equation
background_strain = Field(KernelFunctionOperation{Center, Center, Center}(background_strain_func, grid, clock, fields, (; ), sp))
strain = Field(KernelFunctionOperation{Center, Center, Center}(strain_func, grid, clock, fields, (; ), sp))
shear = Field(KernelFunctionOperation{Center, Center, Center}(shear_func, grid, clock, fields, (; ), sp))
divergence = Field(KernelFunctionOperation{Center, Center, Center}(divergence_func, grid, clock, fields, (; ), sp))

# Subgrid terms
subgrid = Field(KernelFunctionOperation{Center, Center, Center}(subgrid_func, grid, clock, fields, (; ), sp))

# Note that ∇bD∇bDt is computed before ∇b²
frontogenesis_fields = (; 
    ∇bD∇bDt, 
    ∇b², 
    background_strain,
    strain,
    shear,
    divergence,
    subgrid
)

(
    ∇bD∇bDt_dfm, 
    ∇b²_dfm, 
    background_strain_dfm,
    strain_dfm,
    shear_dfm,
    divergence_dfm,
    subgrid_dfm
) = dfm.(frontogenesis_fields)

frontogenesis_fields_dfm = (;
    ∇bD∇bDt_dfm, 
    ∇b²_dfm, 
    background_strain_dfm,
    strain_dfm,
    shear_dfm,
    divergence_dfm,
    subgrid_dfm
)

dependency_fields = merge(frontogenesis_fields, frontogenesis_fields_dfm)
output_fields = frontogenesis_fields_dfm
