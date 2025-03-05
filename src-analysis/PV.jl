using Oceanostics
using Oceananigans.Operators

# Potential vorticity q and its evolution equation
# q is calculated as the PV of the mean state:
#  q = (∂v/∂x + f) * ∂b/∂z - ∂v/∂z * ∂b/∂x.
# The evolution equation is as follows:
#  ∂q/∂t = qA + qV + qB,
# where
#  qA = (-u_dfm + αx) * ∂q/∂x - w_dfm * ∂q/∂z
#  qV = ∂V/∂x * ∂b/∂z - ∂V/∂z * ∂b/∂x
#  qB = (∂v/∂x + f) * ∂B/∂z - ∂v/∂z * ∂B/∂x
# where V and B are the non-conservative terms in the velocity 
# and buoyancy equations

include("terms/terms.jl")
include("terms/end_terms.jl")
#include("terms/fluctuation_terms.jl")
#include("terms/sgs_terms.jl")

u_dfm = dfm(u)
v_dfm = dfm(v)
w_dfm = dfm(w)
b_dfm = dfm(b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)

#=
    Values of the fields at the end of the domain
=#

η_x_end_op = KernelFunctionOperation{Center, Nothing, Center}(η_x_end_func, grid, (; u, v, w), (; u=u_dfm, v=v_dfm, w=w_dfm))
η_z_end_op = KernelFunctionOperation{Center, Nothing, Center}(η_z_end_func, grid, (; u, v, w), (; u=u_dfm, v=v_dfm, w=w_dfm))
b_end_op = KernelFunctionOperation{Center, Nothing, Center}(center_end_func, grid, b, b_dfm)

η_x_end = Field(η_x_end_op)
η_z_end = Field(η_z_end_op)
b_end = Field(b_end_op)

end_fields = (; η_x_end, η_z_end, b_end)

#=
    Potential vorticity
=#

@inline function q_func(i, j, k, grid, v_dfm, b_dfm, f)
    
    a = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, v_dfm)
    b = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_dfm)
    
    c = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, v_dfm)
    d = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b_dfm)
    
    return d * (a + f) - b * c
end

q_op = KernelFunctionOperation{Center, Nothing, Center}(q_func, grid, v_dfm, b_dfm, sp.f)
q = Field(q_op)

# Horizontal transport
@inline function uq_func(i, j, k, grid, u_dfm, q)
    return @inbounds u_dfm[i, j, k] * ℑxᶠᵃᵃ(i, j, k, grid, q)
end

uq_op = KernelFunctionOperation{Face, Nothing, Center}(uq_func, grid, u_dfm, q)
uq = Field(uq_op)

# Vertical transport
@inline function wq_func(i, j, k, grid, w_dfm, q)
    return @inbounds w_dfm[i, j, k] * ℑzᵃᵃᶠ(i, j, k, grid, q)
end

wq_op = KernelFunctionOperation{Center, Nothing, Face}(wq_func, grid, w_dfm, q)
wq = Field(wq_op)

# Change in q due to strain
@inline function αq_func(i, j, k, grid, v_dfm, b_dfm, q, sp, t, η_x_end, η_z_end, b_end)
    
    α = @inbounds strain(grid.xᶜᵃᵃ[i], t[1], sp)
    
    a = α * ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_end) * ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, v_dfm)
    b = -α * ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b_end) * ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, v_dfm)
    
    c = @inbounds -α * η_x_end[i, j, k] * ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b_end)
    d = @inbounds -α * η_z_end[i, j ,k] * ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, b_end)
    
    
    αx = @inbounds α * grid.xᶜᵃᵃ[i]
    
    e = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, q)
    
    return a + b + c + d + αx * e
end

αq_op = KernelFunctionOperation{Center, Nothing, Center}(αq_func, grid, v_dfm, b_dfm, q, sp, t, η_x_end, η_z_end, b_end)
αq = Field(αq_op)

# Total conservative change in q
@inline function qA_func(i, j, k, grid, uq, wq, αq)
    return @inbounds -(∂xᶜᶜᶜ(i, j, k, grid, uq) + ∂zᶜᶜᶜ(i, j, k, grid, wq)) + αq[i, j, k]
end

# Non-conservative change to be calculated in post

qA_op = KernelFunctionOperation{Center, Nothing, Center}(qA_func, grid, uq, wq, αq)
qA = Field(qA_op)

outputs = (; q, αq, uq, wq, qA)

function update_outputs!(outputs)
    map(compute!, mean_fields)
    map(compute!, end_fields)
    map(compute!, outputs)
end
