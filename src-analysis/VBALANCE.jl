include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")

include("terms/vbalance/turbulence.jl")

v_dfm = dfm(input_fields.v)

mean_fields = (; v_dfm)

DvDt3D = Field(KernelFunctionOperation{Center, Face, Center}(DvDt_func, grid, clock, input_fields, mean_fields, sp))
turbulence_h3D = Field(KernelFunctionOperation{Center, Face, Center}(turbulence_h_func, grid, clock, input_fields, mean_fields, sp))
turbulence_z3D = Field(KernelFunctionOperation{Center, Face, Center}(turbulence_z_func, grid, clock, input_fields, mean_fields, sp))

turbulence3D_fields = (; DvDt3D, turbulence_h3D, turbulence_z3D)

(DvDt, turbulence_h, turbulence_z) = map(dfm, turbulence3D_fields)

turbulence_fields = (; DvDt, turbulence_h, turbulence_z)

dependency_fields = merge(mean_fields, turbulence3D_fields, turbulence_fields)
output_fields = turbulence_fields

skip_update = (:u_next, :w_next, :b_next, :pNHS, :pNHS_next)