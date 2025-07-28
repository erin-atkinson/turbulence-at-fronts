include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")

include("terms/meanfrontogenesis/dbdx.jl")
include("terms/meanfrontogenesis/advection.jl")
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

dbdx = Field(KernelFunctionOperation{Center, Nothing, Center}(dbdx_func, grid, clock, fields, mean_fields, sp))

# Resolved terms in frontogenesis equation
background_strain = Field(KernelFunctionOperation{Center, Nothing, Center}(background_strain_func, grid, clock, fields, mean_fields, sp))
divergence = Field(KernelFunctionOperation{Center, Nothing, Center}(divergence_func, grid, clock, fields, mean_fields, sp))
tilting = Field(KernelFunctionOperation{Center, Nothing, Center}(tilting_func, grid, clock, fields, mean_fields, sp))

# Turbulence
turbulence_h = Field(KernelFunctionOperation{Center, Center, Center}(turbulence_h_func, grid, clock, fields, mean_fields, sp))
turbulence_z = Field(KernelFunctionOperation{Center, Center, Center}(turbulence_z_func, grid, clock, fields, mean_fields, sp))

# Advection
adv_x = Field(KernelFunctionOperation{Center, Center, Center}(adv_x_func, grid, clock, fields, (; dbdx), sp))
adv_y = Field(KernelFunctionOperation{Center, Center, Center}(adv_y_func, grid, clock, fields, (; dbdx), sp))
adv_z = Field(KernelFunctionOperation{Center, Center, Center}(adv_z_func, grid, clock, fields, (; dbdx), sp))

# Subgrid terms
# subgrid = Field(KernelFunctionOperation{Center, Center, Center}(subgrid_func, grid, clock, fields, mean_fields, sp))

frontogenesis_fields = (;
    dbdx, 
    background_strain,
    divergence,
    tilting
)

frontogenesis_fields_3D = (;
    turbulence_h,
    turbulence_z,
    adv_x,
    adv_y,
    adv_z
)

turbulence_h_dfm = dfm(turbulence_h)
turbulence_z_dfm = dfm(turbulence_z)

adv_x_dfm = dfm(adv_x)
adv_y_dfm = dfm(adv_y)
adv_z_dfm = dfm(adv_z)

frontogenesis_fields_dfm = (;
    turbulence_h_dfm,
    turbulence_z_dfm,
    adv_x_dfm,
    adv_y_dfm,
    adv_z_dfm
)

dependency_fields = merge(mean_fields, frontogenesis_fields, frontogenesis_fields_3D, frontogenesis_fields_dfm)
output_fields = merge(frontogenesis_fields, frontogenesis_fields_dfm)

frames = 891:901
