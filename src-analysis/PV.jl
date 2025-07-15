include("terms/terms.jl")
include("terms/vorticity/vorticity.jl")
include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")


include("terms/pv/q.jl")
include("terms/pv/uq.jl")
include("terms/pv/Vq.jl")

q = Field(KernelFunctionOperation{Center, Center, Center}(q_func, grid, clock, fields, (; ), sp))

uq = Field(KernelFunctionOperation{Face, Center, Center}(uq_func, grid, clock, fields, (; q), sp))
wq = Field(KernelFunctionOperation{Center, Center, Face}(wq_func, grid, clock, fields, (; q), sp))
Vq = Field(KernelFunctionOperation{Center, Center, Center}(Vq_func, grid, clock, fields, (; q), sp))

q_dfm = dfm(q)
uq_dfm = dfm(uq)
wq_dfm = dfm(wq)
Vq_dfm = dfm(Vq)

dependency_fields = (; q, uq, wq, Vq, q_dfm, uq_dfm, wq_dfm, Vq_dfm)
output_fields = (; q_dfm, uq_dfm, wq_dfm, Vq_dfm)
