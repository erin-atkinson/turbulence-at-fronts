include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")

include("terms/Mbalance/M.jl")
include("terms/Mbalance/DMDt.jl")
include("terms/Mbalance/Q.jl")
include("terms/Mbalance/F.jl")

u_dfm = dfm(input_fields.u)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)
b_next_dfm = dfm(input_fields.b_next)

M² = Field(KernelFunctionOperation{Face, Nothing, Center}(M²_func, grid, clock, input_fields, (; b_dfm), sp))
mean_fields = (; u_dfm, w_dfm, b_dfm, b_next_dfm, M²)

u∂M²∂x3D = Field(KernelFunctionOperation{Face, Center, Center}(u∂M²∂x_func, grid, clock, input_fields, mean_fields, sp))
w∂M²∂z3D = Field(KernelFunctionOperation{Face, Center, Center}(w∂M²∂z_func, grid, clock, input_fields, mean_fields, sp))
u∂M²∂x = dfm(u∂M²∂x3D)
w∂M²∂z = dfm(w∂M²∂z3D)

∂u∂xM² = Field(KernelFunctionOperation{Face, Nothing, Center}(∂u∂xM²_func, grid, clock, input_fields, mean_fields, sp))
∂w∂xN² = Field(KernelFunctionOperation{Face, Nothing, Center}(∂w∂xN²_func, grid, clock, input_fields, mean_fields, sp))
M²_aux_fields = (; mean_fields..., u∂M²∂x3D, w∂M²∂z3D, u∂M²∂x, w∂M²∂z, ∂u∂xM², ∂w∂xN²)

DM²Dt = Field(KernelFunctionOperation{Face, Nothing, Center}(DM²Dt_func, grid, clock, input_fields, M²_aux_fields, sp))
Fx3D = Field(KernelFunctionOperation{Face, Center, Center}(Fx_func, grid, clock, input_fields, M²_aux_fields, sp))
Fz3D = Field(KernelFunctionOperation{Face, Center, Center}(Fz_func, grid, clock, input_fields, M²_aux_fields, sp))
Fx = dfm(Fx3D)
Fz = dfm(Fz3D)

dependency_fields = (; M²_aux_fields..., DM²Dt, Fx3D, Fz3D, Fx, Fz)

output_fields = (; M², DM²Dt, u∂M²∂x, w∂M²∂z, ∂u∂xM², ∂w∂xN², Fx, Fz)

skip_update = (:u_next, :v_next, :w_next, :pNHS, :pNHS_next)
