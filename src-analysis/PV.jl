include("terms/terms.jl")
include("terms/vorticity/vorticity.jl")
include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")

include("terms/pv/q.jl")
include("terms/pv/uq.jl")
include("terms/pv/Vq.jl")

u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)

mean_fields = (; u = u_dfm, v = v_dfm, w = w_dfm, b = b_dfm)

q = Field(KernelFunctionOperation{Center, Nothing, Center}(q_func, grid, clock, mean_fields, (; ), sp))

dependency_fields = (; q)
output_fields = (; q)

skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)
