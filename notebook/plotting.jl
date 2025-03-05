using CairoMakie
using JLD2
using ImageFiltering: imfilter, Kernel.gaussian
using OffsetArrays: no_offset_view

# Some things for making plots

@inline function bin_counts(data, bins)
    map(bins[1:end-1], bins[2:end]) do l, r
        count(x->l<=x<r, data)
    end
end

@inline filt(a, σ::Tuple{T, T} where T<:Number) = imfilter(a, gaussian(σ))
@inline filt(a, σ::Number) = filt(a, (σ, σ)) 
@inline filt(a, σ1::Number, σ2::Number) = filt(a, (σ1, σ2))

@inline function grid_cells(foldername; halo=false)
   jldopen("$foldername/output.jld2") do file
        
        xsᶜ = no_offset_view(file["grid/xᶜᵃᵃ"])
        xsᶠ = no_offset_view(file["grid/xᶠᵃᵃ"])

        ysᶜ = no_offset_view(file["grid/yᵃᶜᵃ"])
        ysᶠ = no_offset_view(file["grid/yᵃᶠᵃ"])

        zsᶜ = no_offset_view(file["grid/zᵃᵃᶜ"])
        zsᶠ = no_offset_view(file["grid/zᵃᵃᶠ"])
        
        if !halo
            xsᶜ = xsᶜ[6:end-5]
            xsᶠ = xsᶠ[6:end-5]
            ysᶜ = ysᶜ[6:end-5]
            ysᶠ = ysᶠ[6:end-5]
            zsᶜ = zsᶜ[6:end-5]
            zsᶠ = zsᶠ[6:end-5]
        end
        
        xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ
    end 
end

@inline function iterations_times(foldername)
   jldopen("$foldername/output.jld2") do file

        iterations = keys(file["timeseries/t"])
        ts = map(iteration->file["timeseries/t/$iteration"], iterations)

        iterations, ts
    end 
end

@inline function simulation_parameters(foldername)
   jldopen("$foldername/parameters.jld2") do file
        file["simulation"]
    end
end

@inline function get_field(file, field, iteration)
    data = no_offset_view(file["timeseries/$field/$iteration"])
    # Remove singleton dimensions
    indices = map(i-> i==1 ? 1 : 6:i-5, size(data))
    data[indices...]
end

@inline function get_field(filename::String, field, iteration)
    jldopen(file->get_field(file, field, iteration), filename)
end

@inline function get_field(f::Function, args...)
    f(get_field(args...))
end

@inline function time_average_field(filename, field, iterations, ts, l, r)
    mapreduce(.+, iterations[l:r]) do iteration
        get_field(filename, field, iteration)
    end ./ (ts[r] - ts[l])
end

@inline function time_average_field(f::Function, args...)
    f(time_average_field(args...))
end

@inline function mean_field(filename, field, iteration, sp)
    get_field(filename, field, iteration) do a
        mean(a[(sp.Nx-sp.Nₕ)÷2:(sp.Nx+sp.Nₕ)÷2, :])
    end
end

@inline strain_function(t, x, α, f, ℓ) = -α * (1 - exp(-max(f * t, 0.0) / 15)) * (tanh((x-4ℓ) / ℓ) + tanh((-x-4ℓ) / ℓ)) / 2

@inline strain_function(t, x, sp) = strain_function(t, x, sp.α, sp.f, sp.ℓ)

@inline function field_timeseries(filename, field, iterations, sp)
    map(iterations) do iteration
        mean_field(filename, field, iteration, sp)
    end
end

@inline function hovmöller_field(func, filename, field, iterations)
    # Reduce field to 1D using func, then get a timeseries
    data = jldopen(filename) do file
        map(iteration->get_field(func, file, field, iteration), iterations)
    end
    return [data[i][j] for i in 1:length(data), j in 1:length(data[1])]
end

@inline function box_limits(foldername, limits)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_cells(foldername)

    l = argmin(abs.(xsᶜ .- limits[1]))
    r = argmin(abs.(xsᶜ .- limits[2]))
    b = argmin(abs.(zsᶜ .- limits[3]))
    t = argmin(abs.(zsᶜ .- limits[4]))
    
    return (l, r, b, t)
end

@inline function box_limits(foldername, l::Number, r::Number, b::Number, t::Number)
    box_limits(foldername, (l, r, b, t))
end

function x_surface_sum(iteration, file, fieldname, z_limits, x, Δz)
    get_field(file, fieldname, iteration) do field
        sum(field[x, z_limits[1]:z_limits[2]]) * Δz
    end
end

function z_surface_sum(iteration, file, fieldname, x_limits, z, Δx)
    get_field(file, fieldname, iteration) do field
        sum(field[x_limits[1]:x_limits[2], z]) * Δx
    end
end

function box_sum(iteration, file, fieldname, limits, ΔV)
    get_field(file, fieldname, iteration) do field
        sum(filt(field, 1, 1)[limits[1]:limits[2], limits[3]:limits[4]]) * ΔV
    end
end

function box_timeseries(iterations, filename, fieldname, limits, ΔV)
    jldopen(filename) do file
        map(iteration->box_sum(iteration, file, fieldname, limits, ΔV), iterations)
    end
end

function save_3d(foldername, field, n)
    
    iterations, ts = iterations_times(foldername)
    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_cells(foldername)
    sp = simulation_parameters(foldername)
    
    data = jldopen(file->file["timeseries/$field/$(iterations[n])"], "$foldername/output.jld2")
    
    filename = "$(rsplit(foldername, '/'; limit=2)[2])-$field-$n.jld2"
    
    jldopen(filename, "w") do file
        file["xsᶜ"] = xsᶜ
        file["ysᶜ"] = ysᶜ
        file["zsᶜ"] = zsᶜ
        file["xsᶠ"] = xsᶠ
        file["ysᶠ"] = ysᶠ
        file["zsᶠ"] = zsᶠ
        file["t"] = ts[n]
        file["sp"] = sp
        file["data"] = data
    end
    return nothing
end