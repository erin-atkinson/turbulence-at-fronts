include("terms/terms.jl")
include("terms/vorticity/vorticity.jl")
include("terms/advection/advection.jl")


include("terms/pv/q.jl")
include("terms/pv/uq.jl")
include("terms/pv/Vq.jl")

q = Field(KernelFunctionOperation{Center, Center, Center}(q_func, grid, clock, fields, (; ), sp))
uq = Field(KernelFunctionOperation{Face, Center, Center}(uq_func, grid, clock, fields, (; q), sp))
Vq = Field(KernelFunctionOperation{Center, Center, Center}(Vq_func, grid, clock, fields, (; q), sp))

q_dfm = dfm(q)
uq_dfm = dfm(uq)
Vq_dfm = dfm(Vq)

dependency_fields = (; q, uq, Vq, q_dfm, uq_dfm, Vq_dfm)
output_fields = (; q_dfm, Vq_dfm)
