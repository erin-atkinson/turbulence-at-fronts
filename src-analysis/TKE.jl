include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")

include("terms/tke/tke.jl")
include("terms/tke/bflux.jl")
include("terms/tke/dissipation.jl")
include("terms/tke/dsp.jl")
include("terms/tke/lsp.jl")
include("terms/tke/vsp.jl")

u_dfm = dfm(rawfields.u)
v_dfm = dfm(rawfields.v)
w_dfm = dfm(rawfields.w)
u_next_dfm = dfm(rawfields.u_next)
v_next_dfm = dfm(rawfields.v_next)
w_next_dfm = dfm(rawfields.w_next)
b_dfm = dfm(rawfields.b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm, u_next_dfm, v_next_dfm, w_next_dfm)

LSP3D = Field(KernelFunctionOperation{Center, Center, Center}(LSP3D_func, grid, clock, rawfields, mean_fields, sp))
VSP3D = Field(KernelFunctionOperation{Center, Center, Center}(VSP3D_func, grid, clock, rawfields, mean_fields, sp))
BFLUX3D = Field(KernelFunctionOperation{Center, Center, Center}(BFLUX3D_func, grid, clock, rawfields, mean_fields, sp))
DSP3D = Field(KernelFunctionOperation{Center, Center, Center}(DSP3D_func, grid, clock, rawfields, mean_fields, sp))
TKE3D = Field(KernelFunctionOperation{Center, Center, Center}(TKE3D_func, grid, clock, rawfields, mean_fields, sp))
DTKEDt3D = Field(KernelFunctionOperation{Center, Center, Center}(DTKEDt3D_func, grid, clock, rawfields, mean_fields, sp))

TKE3D_fields = (; LSP3D, VSP3D, BFLUX3D, DSP3D, TKE3D, DTKEDt3D)

(LSP, VSP, BFLUX, DSP, TKE, DTKEDt) = map(dfm, TKE3D_fields)

TKE_fields = (; LSP, VSP, BFLUX, DSP, TKE, DTKEDt)

# ε = Field(KernelFunctionOperation{Center, Nothing, Center}(ε_func, grid, clock, rawfields, TKE_fields, sp))

# TKE_fields = merge(TKE_fields, (; ε))

dependency_fields = merge(mean_fields, TKE3D_fields, TKE_fields)
output_fields = TKE_fields
