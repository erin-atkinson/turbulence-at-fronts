using Oceanostics
using Oceananigans.Operators
using ImageFiltering: imfilter, Kernel.gaussian

# Steady state potential vorticity and its transport
# The average u, v, w, b fields are first calculated over the whole interval
# Then the potential vorticity q is calculated from these
# And transport is given by uq and wq, which also gives non-conservative changes

include("terms/terms.jl")
include("terms/potential_vorticity.jl")

u_dfm = dfm(u)
v_dfm = dfm(v)
w_dfm = dfm(w)
b_dfm = dfm(b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)

N = length(ts)

function update_steady!(steady_field, field, N)
    steady_field .+= (field ./ N)
    return nothing
end

u_steady = Field{Face, Nothing, Center}(grid)
v_steady = Field{Center, Nothing, Center}(grid)
w_steady = Field{Center, Nothing, Face}(grid)
b_steady = Field{Center, Nothing, Center}(grid)

steady_fields = (; u_steady, v_steady, w_steady, b_steady)

ψ = Field(@at (Center, Nothing, Center) CumulativeIntegral(-u_steady; dims=3))

# Potential vorticity operations

q_op = KernelFunctionOperation{Center, Nothing, Center}(q_func, grid, v_steady, b_steady, sp.f)
q = Field(q_op)

UADV_op = KernelFunctionOperation{Center, Nothing, Center}(UADV_func, grid, u_steady, q)
UADV = Field(UADV_op)

WADV_op = KernelFunctionOperation{Center, Nothing, Center}(WADV_func, grid, w_steady, q)
WADV = Field(WADV_op)

# Horizontal transport
@inline function uq_func(i, j, k, grid, u_dfm, q)
    return @inbounds u_dfm[i, j, k] * ℑxᶠᵃᵃ(i, j, k, grid, q)
end

uq_op = KernelFunctionOperation{Face, Nothing, Center}(uq_func, grid, u_steady, q)
uq = Field(uq_op)

# Vertical transport
@inline function wq_func(i, j, k, grid, w_dfm, q)
    return @inbounds w_dfm[i, j, k] * ℑzᵃᵃᶠ(i, j, k, grid, q)
end

wq_op = KernelFunctionOperation{Center, Nothing, Face}(wq_func, grid, w_steady, q)
wq = Field(wq_op)

strain_op = KernelFunctionOperation{Center, Nothing, Center}(strain_func, grid, q, sp, t)
αADV = Field(strain_op)

Q_op = KernelFunctionOperation{Center, Nothing, Center}(Q_func, grid, UADV, WADV, αADV)
Q = Field(Q_op)

outputs = (; )

function update_outputs!(outputs)
    map(compute!, mean_fields)
    update_steady!(u_steady, u_dfm, N)
    update_steady!(v_steady, v_dfm, N)
    update_steady!(w_steady, w_dfm, N)
    update_steady!(b_steady, b_dfm, N)
end

macro postpostprocess()
    outputs = (; q, UADV, WADV, αADV, Q, uq, wq)
    map(compute!, outputs)
    write_outputs(filename, outputs, 0)
    compute!(ψ)
    write_outputs(filename, (; u=u_steady, v=v_steady, w=w_steady, b=b_steady, ψ), 0)
end