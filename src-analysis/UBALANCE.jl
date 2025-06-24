include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")

include("terms/ubalance/usq.jl")
include("terms/ubalance/coriolis.jl")
include("terms/ubalance/background_strain.jl")
include("terms/ubalance/pressure.jl")

include("terms/ubalance/turbulence.jl")
include("terms/ubalance/subgrid.jl")


u_dfm = dfm(fields.u)
v_dfm = dfm(fields.v)
w_dfm = dfm(fields.w)
p_dfm = dfm(fields.p)

mean_fields = (; u_dfm, v_dfm, w_dfm, p_dfm)

u² = Field(KernelFunctionOperation{Face, Center, Center}(u²_func, grid, clock, fields, mean_fields, sp))

# Resolved terms in u²/2 equation
background_strain = Field(KernelFunctionOperation{Face, Center, Center}(background_strain_func, grid, clock, fields, mean_fields, sp))
pressure = Field(KernelFunctionOperation{Face, Center, Center}(pressure_func, grid, clock, fields, mean_fields, sp))
coriolis = Field(KernelFunctionOperation{Face, Center, Center}(coriolis_func, grid, clock, fields, mean_fields, sp))

# Turbulence
turbulence_h = Field(KernelFunctionOperation{Face, Center, Center}(turbulence_h_func, grid, clock, fields, mean_fields, sp))
turbulence_z = Field(KernelFunctionOperation{Face, Center, Center}(turbulence_z_func, grid, clock, fields, mean_fields, sp))

# Subgrid terms
subgrid = Field(KernelFunctionOperation{Face, Center, Center}(subgrid_func, grid, clock, fields, mean_fields, sp))

u²_fields = (;
    u², 
    background_strain,
    pressure,
    coriolis,
    turbulence_h,
    turbulence_z,
    subgrid
)

(
    u²_dfm, 
    background_strain_dfm,
    pressure_dfm,
    coriolis_dfm,
    turbulence_h_dfm,
    turbulence_z_dfm,
    subgrid_dfm
) = map(dfm, u²_fields)

u²_fields_dfm = (; 
    u²_dfm, 
    background_strain_dfm,
    pressure_dfm,
    coriolis_dfm,
    turbulence_h_dfm,
    turbulence_z_dfm,
    subgrid_dfm
)

dependency_fields = merge(mean_fields, u²_fields, u²_fields_dfm)
output_fields = u²_fields_dfm
