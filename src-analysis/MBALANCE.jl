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

div_UM²3D = Field(KernelFunctionOperation{Face, Center, Center}(div_UM²_func, grid, clock, input_fields, mean_fields, sp))
div_UM² = dfm(div_UM²3D)

∂u∂xM² = Field(KernelFunctionOperation{Face, Nothing, Center}(∂u∂xM²_func, grid, clock, input_fields, mean_fields, sp))
∂w∂xN² = Field(KernelFunctionOperation{Face, Nothing, Center}(∂w∂xN²_func, grid, clock, input_fields, mean_fields, sp))
M²_aux_fields = (; mean_fields..., div_UM²3D, div_UM², ∂u∂xM², ∂w∂xN²)

DM²Dt = Field(KernelFunctionOperation{Face, Nothing, Center}(DM²Dt_func, grid, clock, input_fields, M²_aux_fields, sp))
F3D = Field(KernelFunctionOperation{Face, Center, Center}(F_func, grid, clock, input_fields, M²_aux_fields, sp))
F = dfm(F3D)

dependency_fields = (; M²_aux_fields..., DM²Dt, F3D, F)

output_fields = (; M², DM²Dt, ∂u∂xM², ∂w∂xN², F)

skip_update = (:u_next, :v_next, :w_next, :pNHS, :pNHS_next)
