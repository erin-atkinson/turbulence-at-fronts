using Oceanostics
using Oceananigans.Operators

# Kinetic energy balance in the down-front mean fields

include("terms/terms.jl")
include("terms/fluctuation_terms.jl")
include("terms/sgs_terms.jl")
include("terms/end_terms.jl")
include("terms/energy_terms.jl")

u_dfm = dfm(u)
v_dfm = dfm(v)
w_dfm = dfm(w)
p_dfm = dfm(p)
b_dfm = dfm(b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm, p_dfm)

# Kinetic energy
T_op = KernelFunctionOperation{Center, Nothing, Center}(T_func, grid, mean_fields)
T = Field(T_op)

# Advection
u_adv_op = KernelFunctionOperation{Center, Nothing, Center}(u_adv_func, grid, mean_fields, sp, t, T)
UADV = Field(u_adv_op)

w_adv_op = KernelFunctionOperation{Center, Nothing, Center}(w_adv_func, grid, mean_fields, T)
WADV = Field(w_adv_op)

# Shear production
LSP_op = KernelFunctionOperation{Center, Nothing, Center}(LSP_func, grid, u, v, w, mean_fields)
LSP = Field(LSP_op)

VSP_op = KernelFunctionOperation{Center, Nothing, Center}(VSP_func, grid, u, v, w, mean_fields)
VSP = Field(VSP_op)

DSP_op = KernelFunctionOperation{Center, Nothing, Center}(DSP_func, grid, mean_fields, sp, t)
DSP = Field(DSP_op)

# Pressure work
PWORK_op = KernelFunctionOperation{Center, Nothing, Center}(PWORK_func, grid, mean_fields)
PWORK = Field(PWORK_op)

# Buoyancy flux
BFLUX_op = KernelFunctionOperation{Center, Nothing, Center}(BFLUX_func, grid, mean_fields)
BFLUX = Field(BFLUX_op)

# Sub-grid scale
SGS_op = KernelFunctionOperation{Center, Nothing, Center}(SGS_func, grid, u, v, w, ν, mean_fields)
SGS = Field(SGS_op)

outputs = (; T, UADV, WADV, LSP, VSP, PWORK, BFLUX, DSP, SGS)

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, outputs)
end
