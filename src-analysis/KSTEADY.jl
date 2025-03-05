using Oceanostics
using Oceananigans.Operators
using ImageFiltering: imfilter, Kernel.gaussian

# Steady state diffusivity

include("terms/terms.jl")

u_dfm = dfm(u)
v_dfm = dfm(v)
w_dfm = dfm(w)
b_dfm = dfm(b)

uu_dfm = dfm(u * u)
uv_dfm = dfm(u * v)
uw_dfm = dfm(u * w)
vw_dfm = dfm(v * w)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm, uu_dfm, uv_dfm, uw_dfm, vw_dfm)

N = length(ts)

function update_steady!(steady_field, field, N)
    steady_field .+= (field ./ N)
    return nothing
end

u_steady = Field{Face, Nothing, Center}(grid)
v_steady = Field{Center, Nothing, Center}(grid)
w_steady = Field{Center, Nothing, Face}(grid)
b_steady = Field{Center, Nothing, Center}(grid)

uu_steady = Field{Face, Nothing, Center}(grid)
uv_steady = Field{Center, Nothing, Center}(grid)
uw_steady = Field{Center, Nothing, Face}(grid)
vw_steady = Field{Center, Nothing, Center}(grid)

steady_fields = (; 
    u=u_steady, v=v_steady, w=w_steady, b=b_steady,
    uu=uu_steady, uv=uv_steady, uw=uw_steady, vw=vw_steady
)

# Horizontal viscosity
@inline function Kh_func(i, j, k, grid, steady_fields)
    u = steady_fields.u
    v = steady_fields.v
    uuᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, steady_fields.uu)
    uvᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, steady_fields.uv)
    
    uᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u)
    
    u_x = ∂xᶜᶜᶜ(i, j, k, grid, u)
    v_x = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, v)
    
    SP1 = (uuᶜᶜᶜ - uᶜᶜᶜ * uᶜᶜᶜ) * u_x
    SP2 = @inbounds (uvᶜᶜᶜ - uᶜᶜᶜ * v[i, j, k]) * v_x
    
    return -(SP1 + SP2) / (u_x^2 + v_x^2)
end

# Vertical viscosity
@inline function Kv_func(i, j, k, grid, steady_fields)
    u = steady_fields.u
    v = steady_fields.v
    w = steady_fields.w
    uwᶜᶜᶜ = ℑxzᶜᵃᶜ(i, j, k, grid, steady_fields.uw)
    vwᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, steady_fields.vw)
    
    uᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u)
    wᶜᶜᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w)
    
    u_z = ℑxzᶜᵃᶜ(i, j, k, grid, ∂zᶠᶜᶠ, u)
    v_z = ℑzᵃᵃᶜ(i, j, k, grid, ∂zᶜᶜᶠ, v)
    
    SP1 = (uwᶜᶜᶜ - uᶜᶜᶜ * wᶜᶜᶜ) * u_z
    SP2 = @inbounds (vwᶜᶜᶜ - wᶜᶜᶜ * v[i, j, k]) * v_z
    
    return -(SP1 + SP2) / (u_z^2 + v_z^2)
end


Kh_op = KernelFunctionOperation{Center, Nothing, Center}(Kh_func, grid, steady_fields)
Kh = Field(Kh_op)

Kv_op = KernelFunctionOperation{Center, Nothing, Center}(Kv_func, grid, steady_fields)
Kv = Field(Kv_op)

outputs = (; )

function update_outputs!(outputs)
    map(compute!, mean_fields)
    update_steady!(u_steady, u_dfm, N)
    update_steady!(v_steady, v_dfm, N)
    update_steady!(w_steady, w_dfm, N)
    update_steady!(b_steady, b_dfm, N)
    
    update_steady!(uu_steady, uu_dfm, N)
    update_steady!(uv_steady, uv_dfm, N)
    update_steady!(uw_steady, uw_dfm, N)
    update_steady!(vw_steady, vw_dfm, N)
end

macro postpostprocess()
    outputs = (; Kh, Kv)
    map(compute!, outputs)
    write_outputs(filename, outputs, 0)
    write_outputs(filename, (; u=u_steady, v=v_steady, w=w_steady, b=b_steady), 0)
end