include("terms/terms.jl")
include("terms/vorticity/vorticity.jl")


ηx = Field(KernelFunctionOperation{Center, Face, Face}(ηx_func, grid, clock, fields, (; ), sp))
ηy = Field(KernelFunctionOperation{Face, Center, Face}(ηy_func, grid, clock, fields, (; ), sp))
ηz = Field(KernelFunctionOperation{Face, Face, Center}(ηz_func, grid, clock, fields, (; ), sp))

frames = [401, 581, 761, 761 + 180*2]
dependency_fields = (; ηx, ηy, ηz)
output_fields = (; ηx, ηy, ηz)
