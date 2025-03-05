using Oceanostics
using Oceananigans.Operators
using ImageFiltering: imfilter, Kernel.gaussian

# Steady state potential vorticity and its transport
# The average u, v, w, b fields are first calculated over the whole interval
# Then the potential vorticity q is calculated from these
# And transport is given by uq and wq, which also gives non-conservative changes

include("terms/terms.jl")

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


@inline function DcDt(i, j, k, grid, u, w, args...)
    u = ℑxᶜᵃᵃ(i, j, k, grid, fGg, u, ∂xᶠᶜᶜ, args...)
    w = ℑzᵃᵃᶜ(i, j, k, grid, fGg, w, ∂zᶜᶜᶠ, args...)
    return u + w
end

# Divergence operations
@inline function δ_func(i, j, k, grid, u_steady)
    ∂xᶜᶜᶜ(i, j, k, grid, u_steady)
end

δ_op = KernelFunctionOperation{Center, Nothing, Center}(δ_func, grid, u_steady)
δ = Field(δ_op)

DδDt_op = KernelFunctionOperation{Center, Nothing, Center}(DcDt, grid, u_steady, w_steady, δ)
DδDt = Field(DδDt_op)

# Divergence operations with strain?

@inline function αδ_func(i, j, k, grid, u_steady, sp)
    @inbounds ∂xᶜᶜᶜ(i, j, k, grid, u_steady) + strain(grid.xᶜᵃᵃ[i], 100 / sp.f, sp)
end

@inline function αDcDt(i, j, k, grid, u, w, δ, sp)
    u = ℑxᶜᵃᵃ(i, j, k, grid, fGg, u, ∂xᶠᶜᶜ, δ)
    w = ℑzᵃᵃᶜ(i, j, k, grid, fGg, w, ∂zᶜᶜᶠ, δ)
    
    αx = @inbounds strain(grid.xᶜᵃᵃ[i], 100 / sp.f, sp) * grid.xᶜᵃᵃ[i]
    c = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, δ)
    
    return u + w - αx * c
end

αδ_op = KernelFunctionOperation{Center, Nothing, Center}(αδ_func, grid, u_steady, sp)
αδ = Field(αδ_op)

αDδDt_op = KernelFunctionOperation{Center, Nothing, Center}(αDcDt, grid, u_steady, w_steady, αδ, sp)
αDδDt = Field(αDδDt_op)

outputs = (; )

function update_outputs!(outputs)
    map(compute!, mean_fields)
    update_steady!(u_steady, u_dfm, N)
    update_steady!(v_steady, v_dfm, N)
    update_steady!(w_steady, w_dfm, N)
    update_steady!(b_steady, b_dfm, N)
end

macro postpostprocess()
    outputs = (; δ, DδDt, αδ, αDδDt)
    map(compute!, outputs)
    write_outputs(filename, outputs, 0)
    write_outputs(filename, (; u=u_steady, v=v_steady, w=w_steady, b=b_steady), 0)
end