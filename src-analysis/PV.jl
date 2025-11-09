include("terms/terms.jl")
include("terms/vorticity/vorticity.jl")
include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")


include("terms/pv/q.jl")
include("terms/pv/uq.jl")
include("terms/pv/Vq.jl")

q = Field(KernelFunctionOperation{Center, Center, Center}(q_func, grid, clock, rawfields, (; ), sp))

uq = Field(KernelFunctionOperation{Face, Center, Center}(uq_func, grid, clock, rawfields, (; q), sp))
wq = Field(KernelFunctionOperation{Center, Center, Face}(wq_func, grid, clock, rawfields, (; q), sp))
div_Uq = Field(KernelFunctionOperation{Center, Center, Center}(div_Uq_func, grid, clock, rawfields, (; q), sp))

Vq = Field(KernelFunctionOperation{Center, Center, Center}(Vq_func, grid, clock, rawfields, (; q), sp))

q_dfm = dfm(q)

uq_dfm = dfm(uq)
wq_dfm = dfm(wq)
div_Uq_dfm = dfm(div_Uq)

Vq_dfm = dfm(Vq)

dependency_fields = (; q, uq, wq, div_Uq, Vq, q_dfm, uq_dfm, wq_dfm, div_Uq_dfm, Vq_dfm)
output_fields = (; q_dfm, uq_dfm, wq_dfm, div_Uq_dfm, Vq_dfm)
