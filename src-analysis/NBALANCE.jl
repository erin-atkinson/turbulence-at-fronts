include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")

include("terms/Nbalance/N.jl")
include("terms/Nbalance/DNDt.jl")
include("terms/Nbalance/Q.jl")
include("terms/Nbalance/F.jl")

u_dfm = dfm(input_fields.u)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)
b_next_dfm = dfm(input_fields.b_next)

N² = Field(KernelFunctionOperation{Center, Nothing, Face}(N²_func, grid, clock, input_fields, (; b_dfm), sp))
mean_fields = (; u_dfm, w_dfm, b_dfm, b_next_dfm, N²)

div_UN²3D = Field(KernelFunctionOperation{Center, Center, Face}(div_UN²_func, grid, clock, input_fields, mean_fields, sp))
div_UN² = dfm(div_UN²3D)

∂w∂zN² = Field(KernelFunctionOperation{Center, Nothing, Face}(∂w∂zN²_func, grid, clock, input_fields, mean_fields, sp))
∂u∂zM² = Field(KernelFunctionOperation{Center, Nothing, Face}(∂u∂zM²_func, grid, clock, input_fields, mean_fields, sp))
N²_aux_fields = (; mean_fields..., div_UN²3D, div_UN², ∂w∂zN², ∂u∂zM²)

DN²Dt = Field(KernelFunctionOperation{Center, Nothing, Face}(DN²Dt_func, grid, clock, input_fields, N²_aux_fields, sp))
F3D = Field(KernelFunctionOperation{Center, Center, Face}(F_func, grid, clock, input_fields, N²_aux_fields, sp))
F = dfm(F3D)

dependency_fields = (; N²_aux_fields..., DN²Dt, F3D, F)

output_fields = (; N², DN²Dt, ∂w∂zN², ∂u∂zM², F, div_UN²)

skip_update = (:u_next, :v_next, :w_next, :pNHS, :pNHS_next)
