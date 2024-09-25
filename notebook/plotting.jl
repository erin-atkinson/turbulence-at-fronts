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

@inline function grid_cells(foldername)
   jldopen("$foldername/output.jld2") do file

        xsᶜ = no_offset_view(file["grid/xᶜᵃᵃ"])[4:end-3]
        xsᶠ = no_offset_view(file["grid/xᶠᵃᵃ"])[4:end-3]

        ysᶜ = no_offset_view(file["grid/yᵃᶜᵃ"])[4:end-3]
        ysᶠ = no_offset_view(file["grid/yᵃᶠᵃ"])[4:end-3]

        zsᶜ = no_offset_view(file["grid/zᵃᵃᶜ"])[4:end-3]
        zsᶠ = no_offset_view(file["grid/zᵃᵃᶠ"])[4:end-3]
        
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

@inline function get_field(filename, field, iteration)
    jldopen(file->no_offset_view(file["timeseries/$field/$iteration"][:, 1, :])[4:end-3, 4:end-3], filename)
end

@inline function get_field(f::Function, filename, field, iteration)
    f(get_field(filename, field, iteration))
end

@inline function time_average_field(filename, field, iterations, ts, l, r)
    mapreduce(.+, iterations[l:r]) do iteration
        get_field(filename, field, iteration)
    end ./ (ts[r] - ts[l])
end

@inline function time_average_field(f::Function, filename, field, iterations, ts, l, r)
    f(time_average_field(filename, field, iterations, ts, l, r))
end

@inline function mean_field(filename, field, iteration, sp)
    get_field(filename, field, iteration) do a
        mean(a[(sp.Nx-sp.Nₕ)÷2:(sp.Nx+sp.Nₕ)÷2, :])
    end
end

@inline strain_function(t, x, α, f, ℓ) = -α * (1 - exp(-f * t / 15)) * (tanh((x-4ℓ) / ℓ) + tanh((-x-4ℓ) / ℓ)) / 2

@inline strain_function(t, x, sp) = strain_function(t, x, sp.α, sp.f, sp.ℓ)

@inline function field_timeseries(filename, field, iterations, sp)
    map(iterations) do iteration
        mean_field(filename, field, iteration, sp)
    end
end