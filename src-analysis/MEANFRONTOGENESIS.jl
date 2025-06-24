include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")

include("terms/meanfrontogenesis/dbdx.jl")
include("terms/meanfrontogenesis/background_strain.jl")
include("terms/meanfrontogenesis/divergence.jl")
include("terms/meanfrontogenesis/tilting.jl")

include("terms/meanfrontogenesis/turbulence.jl")
include("terms/meanfrontogenesis/subgrid.jl")


u_dfm = dfm(fields.u)
v_dfm = dfm(fields.v)
w_dfm = dfm(fields.w)
b_dfm = dfm(fields.b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)

∇b² = Field(KernelFunctionOperation{Center, Center, Center}(∇b²_func, grid, clock, fields, mean_fields, sp))

# Resolved terms in frontogenesis equation
background_strain = Field(KernelFunctionOperation{Center, Center, Center}(background_strain_func, grid, clock, fields, mean_fields, sp))
divergence = Field(KernelFunctionOperation{Center, Center, Center}(divergence_func, grid, clock, fields, mean_fields, sp))
tilting = Field(KernelFunctionOperation{Center, Center, Center}(tilting_func, grid, clock, fields, mean_fields, sp))

# Turbulence
turbulence_h = Field(KernelFunctionOperation{Center, Center, Center}(turbulence_h_func, grid, clock, fields, mean_fields, sp))
turbulence_z = Field(KernelFunctionOperation{Center, Center, Center}(turbulence_z_func, grid, clock, fields, mean_fields, sp))

# Subgrid terms
subgrid = Field(KernelFunctionOperation{Center, Center, Center}(subgrid_func, grid, clock, fields, mean_fields, sp))

frontogenesis_fields = (;
    ∇b², 
    background_strain,
    divergence,
    tilting,
    turbulence_h,
    turbulence_z,
    subgrid
)

(
    ∇b²_dfm, 
    background_strain_dfm,
    divergence_dfm,
    tilting_dfm,
    turbulence_h_dfm,
    turbulence_z_dfm,
    subgrid_dfm
) = map(dfm, frontogenesis_fields)

frontogenesis_fields_dfm = (;
    ∇b²_dfm, 
    background_strain_dfm,
    divergence_dfm,
    tilting_dfm,
    turbulence_h_dfm,
    turbulence_z_dfm,
    subgrid_dfm
)

dependency_fields = merge(mean_fields, frontogenesis_fields, frontogenesis_fields_dfm)
output_fields = frontogenesis_fields_dfm
