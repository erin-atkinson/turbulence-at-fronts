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

function w′u′ᶜᶜᶜ_kernel_func(i, j, k, grid, v, v_dfm, w, w_dfm)
    wᶜᶜᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w)
    w_dfmᶜᶜᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w_dfm)
    uᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u)
    u_dfmᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u_dfm)
    return a′b′(i, j, k, grid, wᶜᶜᶜ, w_dfmᶜᶜᶜ, uᶜᶜᶜ, u_dfmᶜᶜᶜ)
end

# Buoyancy correlations on respective velocity location
function u′b′ᶜᶜᶠ_kernel_func(i, j, k, grid, u, u_dfm, b, b_dfm)
    bᶠᶜᶜ = ℑxᶠᵃᵃ(i, j, k, grid, b)
    b_dfmᶠᶜᶜ = ℑxᶠᵃᵃ(i, j, k, grid, b_dfm)
    return a′b′(i, j, k, grid, wᶜᶜᶜ, w_dfmᶜᶜᶜ, bᶠᶜᶜ, b_dfmᶠᶜᶜ)
end

function v′b′ᶜᶜᶠ_kernel_func(i, j, k, grid, v, v_dfm, b, b_dfm)
    bᶜᶠᶜ = ℑyᵃᶠᵃ(i, j, k, grid, b)
    return a′b′(i, j, k, grid, v, v_dfm, bᶜᶠᶜ, b_dfm)
end

function w′b′ᶜᶜᶠ_kernel_func(i, j, k, grid, w, w_dfm, b, b_dfm)
    bᶜᶜᶠ = ℑzᵃᵃᶠ(i, j, k, grid, b)
    b_dfmᶜᶜᶠ = ℑzᵃᵃᶠ(i, j, k, grid, b_dfm)
    return a′b′(i, j, k, grid, w, w_dfm, bᶜᶜᶠ, b_dfmᶜᶜᶠ)
end

# Function to update input fields
function update_fields!(input_fields, frame, path)
    jldopen(path) do file
        input_fields.u.data .= file["timeseries/u/$frame"]
        input_fields.v.data .= file["timeseries/v/$frame"]
        input_fields.w.data .= file["timeseries/w/$frame"]
        input_fields.b.data .= file["timeseries/b/$frame"]
        input_fields.φ.data .= file["timeseries/φ/$frame"]
        map(fill_halo_regions!, input_fields)
    end
    return nothing
end

function write_output(mean_fields, correlation_fields, frame, path)
    jldopen(path, "w") do file
        file["uᶠⁿᶜ/$frame"] = mean_fields.uᶠⁿᶜ
        file["vᶜⁿᶜ/$frame"] = mean_fields.vᶜⁿᶜ
        file["wᶜⁿᶠ/$frame"] = mean_fields.wᶜⁿᶠ
        file["bᶜⁿᶜ/$frame"] = mean_fields.bᶜⁿᶜ
        file["φᶜⁿᶜ/$frame"] = mean_fields.φᶜⁿᶜ
        
    end
    return nothing
end