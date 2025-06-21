include("terms/terms.jl")

include("terms/pv/q.jl")
include("terms/pv/Vq.jl")

q = Field(KernelFunctionOperation{Center, Center, Center}(q_func, grid, clock, fields, (; ), sp))
Vq = Field(KernelFunctionOperation{Center, Center, Center}(Vq_func, grid, clock, fields, (; q), sp))

q_dfm = dfm(q)
Vq_dfm = dfm(Vq)

dependency_fields = (; q, Vq, q_dfm, Vq_dfm)
output_fields = (; q_dfm, Vq_dfm)
