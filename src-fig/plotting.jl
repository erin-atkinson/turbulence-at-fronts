# plotting.jl
# Various helper functions for making figures
# -------------------------------------------------------------


# -------------------------------------------------------------
using JLD2
using ImageFiltering: imfilter, Kernel.gaussian
using OffsetArrays: no_offset_view
using CairoMakie
using CairoMakie.Colors
# -------------------------------------------------------------

# -------------------------------------------------------------
const b_step = 0.002 / 15
const b_levels = range(-0.01, -0.004, 60)

const ψ_levels = range(-1, 1, 10)
# -------------------------------------------------------------

# -------------------------------------------------------------
nov = no_offset_view

# For changing length and time scales

# Histogram
# Is there not an existing function for this?
@inline function bin_counts(data, bins)
    map(bins[1:end-1], bins[2:end]) do l, r
        count(x->l<=x<r, data)
    end
end

# Shortcut for filtering...
@inline filt(a, σ::Tuple) = imfilter(a, gaussian(σ))
@inline filt(a, σ...) = filt(a, σ)

# If only one σ, repeat for length of a
@inline filt(a::A, σ::Integer) where {T, n, A<:Array{T, n}} = filt(a, (repeat([σ], n)..., ))

@inline halos(file) = file["grid/Hx"], file["grid/Hz"], file["grid/Hz"]

# Take a vector of nd arrays and turn it into a single n+1d array
function expand(v::V) where {T, n, A<:Array{T, n}, V<:Vector{A}}
    l = length(v)
    s = size(v[1])
    new_size = (l, s...)
    a = zeros(T, new_size)
    
    inds = [Colon() for x in s] # There must be some syntax I'm missing
    for i in 1:l
        a[i, inds...] .= v[i]
    end
    a
end

expand(v::V) where {T, V<:Vector{T}} = v
# -------------------------------------------------------------


# -------------------------------------------------------------
function subfig_label!(scene, label; fontsize=14)
    gl = GridLayout(scene, 
        tellwidth = false, 
        tellheight = false, 
        halign = :left, 
        valign = :top,
    )
    Box(gl[1, 1], color = :white, strokecolor = :black, strokewidth = 0)
    Label(gl[1, 1], label;
        fontsize,
        padding = (1, 1, 2, 2),
    )
end

# Pass an integer to get that letter 
function subfig_label!(scene, label::Integer; fontsize=14)
    char = 'a' + (label - 1)
    subfig_label!(scene, "$char)"; fontsize)
end

function prettytime(t; l=3, digits=3)
    a = split(string(round(t; digits)), '.')
    s = t < 0 ? "-" : " "
    return string(s, lpad(a[1], l, '0'), ".", rpad(a[2], digits, '0'))
end
# -------------------------------------------------------------


# -------------------------------------------------------------
# Grid, time and parameter retrieval
@inline function grid_nodes(file; halo=false)
    xsᶜ = nov(file["grid/xᶜᵃᵃ"])
    xsᶠ = nov(file["grid/xᶠᵃᵃ"])

    ysᶜ = nov(file["grid/yᵃᶜᵃ"])
    ysᶠ = nov(file["grid/yᵃᶠᵃ"])

    zsᶜ = nov(file["grid/zᵃᵃᶜ"])
    zsᶠ = nov(file["grid/zᵃᵃᶠ"])

    if !halo
        Hx, Hy, Hz = halos(file)

        xsᶜ = xsᶜ[(Hx+1):(end-Hx)]
        xsᶠ = xsᶠ[(Hx+1):(end-Hx)]
        ysᶜ = ysᶜ[(Hy+1):(end-Hy)]
        ysᶠ = ysᶠ[(Hy+1):(end-Hy)]
        zsᶜ = zsᶜ[(Hz+1):(end-Hz)]
        zsᶠ = zsᶠ[(Hz+1):(end-Hz)]
    end

    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ
end

@inline function grid_nodes(foldername::String; halo=false)
    
    filename = joinpath(foldername, "output.jld2")
    !isfile(filename) && (filename = joinpath(foldername, "DFM.jld2"))
    
    jldopen(filename) do file
        grid_nodes(file; halo)
    end 
end

@inline function iterations_times(file)
    
    iterations = keys(file["timeseries/t"])
    ts = map(iteration->file["timeseries/t/$iteration"], iterations)
    
    iterations, ts
end

@inline function iterations_times(foldername::String)
    
    filename = joinpath(foldername, "output.jld2")
    !isfile(filename) && (filename = joinpath(foldername, "DFM.jld2"))
    
    jldopen(iterations_times, filename)
end

@inline function simulation_parameters(foldername)
    
    filename = joinpath(foldername, "parameters.jld2")
    
    jldopen(filename) do file
        file["simulation"]
    end
end

@inline function centre_indices(foldername)
    sp = simulation_parameters(foldername)
    return ((sp.Nx - sp.Nh) ÷ 2):((sp.Nx + sp.Nh) ÷ 2)
end
# -------------------------------------------------------------


# -------------------------------------------------------------
# Get field from the JLD2 file
@inline function get_field(file, field, iteration; halo=false)
    
    jldpath = "timeseries/$field/$iteration"
    data = nov(file[jldpath])
    
    # Remove singleton dimensions and halo points
    Hs = halo ? zeros(Int, 3) : halos(file)
    indices = map(Hs, size(data)) do H, s
        s==1 ? 1 : (1+H):(s-H)
    end
    
    data[indices...]
end

@inline function get_field(files::AbstractVector, fields::AbstractVector, iteration; halo=false)
    map(files, fields) do file, field
        get_field(file, field, iteration; halo)
    end
end

# Get from a string
@inline function get_field(filename::String, args...; kwargs...)
    jldopen(file->get_field(file, args...; kwargs...), filename)
end

# Apply a function to the field after getting
@inline function get_field(f::Function, args...; kwargs...)
    f(get_field(args...; kwargs...))
end
# -------------------------------------------------------------


# -------------------------------------------------------------
# Reduce Nd field using func (Nd->Md), then find it for all times

@inline function timeseries_of(args...)
    timeseries_of(identity, args...)
end

@inline function timeseries_of(func::Function, file, field, iterations)
    data = map(iteration->get_field(func, file, field, iteration), iterations)
    expand(data)
end

# From a string
@inline function timeseries_of(func::Function, filename::String, args...)
    jldopen(filename) do file
        timeseries_of(func, file, args...)
    end
end
# -------------------------------------------------------------


# -------------------------------------------------------------
# Time average of a function
@inline function time_average_of(func, file, field, iterations)
    mapreduce(.+, iterations) do iteration
        get_field(func, file, field, iteration)
    end ./ length(iterations)
end

# From a string
@inline function time_average_of(func, filename::String, args...)
    jldopen(filename) do file
        time_average_of(func, file, args...)
    end
end
# -------------------------------------------------------------


# -------------------------------------------------------------
for i in [:x, :y, :z], j in [:ᶜ, :ᶠ]
    fname = Symbol(i, j, :bounds)
    nodes = Symbol(i, :s, j)
    @eval function $fname(foldername, args::Tuple)
        xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(foldername)
        map(arg->argmin(abs.($nodes .- arg)), args)
    end
    @eval $fname(foldername, args...) = $fname(foldername, args)
    @eval $fname(foldername, i) = $fname(foldername, (i, ))[1]
end

function tbounds(foldername, args::Tuple)
    iterations, times = iterations_times(foldername)
    map(arg->argmin(abs.(times .- arg)), args)
end
tbounds(foldername, args...) = tbounds(foldername, args)
tbounds(foldername, i) = tbounds(foldername, (i, ))[1]
# -------------------------------------------------------------

# -------------------------------------------------------------
function combine_output(target::AbstractString, source::AbstractString)
    jldopen(target, "a") do targetfile
        jldopen(source) do sourcefile
            combine_output(targetfile, sourcefile)
        end
    end
end

function combine_output(targetfile, sourcefile)
    timeseries = keys(sourcefile["timeseries"])
    for x in timeseries
        combine_timeseries(targetfile, sourcefile, x)
    end
    return nothing
end

function combine_timeseries(targetfile, sourcefile, timeseries)
    iterations = keys(sourcefile["timeseries/$timeseries"])
    for iteration in iterations
        targetfile["timeseries/$timeseries/$iteration"] = targetfile["timeseries/$timeseries/$iteration"]
    end
    return nothing
end
# -------------------------------------------------------------

# -------------------------------------------------------------
function b_weight(b::T, b_start, b_end) where T <: Number
    if b_start <= b < b_end
        return one(T)
    else
        return zero(T)
    end
end

function avg_field_between(field::A, b::A, b₁, b₂) where {T, n, A<:Array{T, n}}
    weights = b_weight.(b, b₁, b₂)
    result = sum(field .* weights; dims=1:(n-1)) ./ sum(weights; dims=1:(n-1))
    return result[:]
end

function field_by_bins(field, b, b_bins)
    map(b_bins[1:end-1], b_bins[2:end]) do b₁, b₂
        field_between(field, b, b₁, b₂)
    end
end
# -------------------------------------------------------------