using Oceananigans
using Oceananigans.Operators
using Oceananigans.BoundaryConditions: fill_halo_regions!
# Defines some functions for use in dfm.jl

# Down-front mean
@inline dfm(a) = Field(Average(a; dims=2))

# Helper function for computing the correlations
@inline a′b′(i, j, k, grid, a, a_dfm, b, b_dfm) = @inbounds (a[i, j, k]-a_dfm[i, j, k]) * (b[i, j, k]-b_dfm[i, j, k])

# Variance of a single field (at same location)
function variance_kernel_func(i, j, k, grid, a, a_dfm)
    return a′b′(i, j, k, grid, a, a_dfm, a, a_dfm)
end

# All velocity cross-correlations will be interpolated onto ᶜᶜᶜ
# Names are cyclic
function u′v′ᶜᶜᶜ_kernel_func(i, j, k, grid, u, u_dfm, v, v_dfm)
    uᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u)
    u_dfmᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u_dfm)
    vᶜᶜᶜ = ℑyᵃᶜᵃ(i, j, k, grid, v)
    return a′b′(i, j, k, grid, uᶜᶜᶜ, u_dfmᶜᶜᶜ, v, v_dfm)
end

function v′w′ᶜᶜᶜ_kernel_func(i, j, k, grid, v, v_dfm, w, w_dfm)
    vᶜᶜᶜ = ℑyᵃᶜᵃ(i, j, k, grid, v)
    wᶜᶜᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w)
    w_dfmᶜᶜᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w_dfm)
    return a′b′(i, j, k, grid, vᶜᶜᶜ, v_dfm, wᶜᶜᶜ, w_dfmᶜᶜᶜ)
end

function w′u′ᶜᶜᶜ_kernel_func(i, j, k, grid, u, u_dfm, w, w_dfm)
    wᶜᶜᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w)
    w_dfmᶜᶜᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w_dfm)
    uᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u)
    u_dfmᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u_dfm)
    return a′b′(i, j, k, grid, wᶜᶜᶜ, w_dfmᶜᶜᶜ, uᶜᶜᶜ, u_dfmᶜᶜᶜ)
end

# Buoyancy correlations on center
function u′b′ᶜᶜᶜ_kernel_func(i, j, k, grid, u, u_dfm, b, b_dfm)
    uᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, b)
    u_dfmᶜⁿᶜ = ℑxᶜᵃᵃ(i, j, k, grid, b_dfm)
    return a′b′(i, j, k, grid, uᶜᶜᶜ, u_dfmᶜⁿᶜ, b, b_dfm)
end

function v′b′ᶜᶜᶜ_kernel_func(i, j, k, grid, v, v_dfm, b, b_dfm)
    vᶜᶜᶜ = ℑyᵃᶜᵃ(i, j, k, grid, v)
    return a′b′(i, j, k, grid, vᶜᶜᶜ, v_dfm, b, b_dfm)
end

function w′b′ᶜᶜᶜ_kernel_func(i, j, k, grid, w, w_dfm, b, b_dfm)
    wᶜᶜᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w)
    w_dfmᶜⁿᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w_dfm)
    return a′b′(i, j, k, grid, wᶜᶜᶜ, w_dfmᶜⁿᶜ, b, b_dfm)
end

function write_output(mean_fields, correlation_fields, frame, path)
    jldopen(path, "a") do file
        for (k, v) in pairs(mean_fields)
            file["$k/$frame"] = Array(v.data[-2:end, 1:1, -2:end])
        end
        for (k, v) in pairs(correlation_fields)
            file["$k/$frame"] = Array(v.data[-2:end, 1:1, -2:end])
        end
    end
    return nothing
end

function write_grid_times(grid, frames, ts, path)
    jldopen(path, "a") do file
        for (i, frame) in enumerate(frames)
            file["t/$frame"] = ts[i]
        end
        # Now copy over grid things so Oceananigans isn't needed
        for k in fieldnames(typeof(grid))
            file["grid/$k"] = getproperty(grid, k)
        end
    end
    return nothing
end