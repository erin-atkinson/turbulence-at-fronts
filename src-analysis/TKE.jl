include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")

include("terms/tke/tke.jl")
include("terms/tke/bflux.jl")
include("terms/tke/dissipation.jl")
include("terms/tke/dsp.jl")
include("terms/tke/lsp.jl")
include("terms/tke/vsp.jl")

u_dfm = dfm(fields.u)
v_dfm = dfm(fields.v)
w_dfm = dfm(fields.w)
b_dfm = dfm(fields.b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)

LSP3D = Field(KernelFunctionOperation{Center, Center, Center}(LSP3D_func, grid, clock, fields, mean_fields, sp))
VSP3D = Field(KernelFunctionOperation{Center, Center, Center}(VSP3D_func, grid, clock, fields, mean_fields, sp))
BFLUX3D = Field(KernelFunctionOperation{Center, Center, Center}(BFLUX3D_func, grid, clock, fields, mean_fields, sp))
DSP3D = Field(KernelFunctionOperation{Center, Center, Center}(DSP3D_func, grid, clock, fields, mean_fields, sp))
ε3D = Field(KernelFunctionOperation{Center, Center, Center}(ε3D_func, grid, clock, fields, mean_fields, sp))
TKE3D = Field(KernelFunctionOperation{Center, Center, Center}(TKE3D_func, grid, clock, fields, mean_fields, sp))

TKE3D_fields = (; LSP3D, VSP3D, BFLUX3D, DSP3D, ε3D, TKE3D)

(LSP, VSP, BFLUX, DSP, ε, TKE) = dfm.(TKE3D_fields)

TKE_fields = (; LSP, VSP, BFLUX, DSP, ε, TKE)

dependency_fields = merge(mean_fields, TKE3D_fields, TKE_fields)
output_fields = TKE_fields
