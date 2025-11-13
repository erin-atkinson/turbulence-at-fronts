# plotting.jl
# Various helper functions for making figures
# -------------------------------------------------------------

# -------------------------------------------------------------
using JLD2
using ImageFiltering: imfilter, Kernel.gaussian
using OffsetArrays: no_offset_view
using Printf
# -------------------------------------------------------------

# -------------------------------------------------------------
const b_step = 0.002 / 15
const b_levels = range(-0.01, -0.004, 60)
# -------------------------------------------------------------

# -------------------------------------------------------------
function iterations_times(file)
    iterations = keys(file["timeseries/t"])
    times = [file["timeseries/t/$i"] for i in iterations]

    return iterations, times
end

iterations_times(filename::String) = jldopen(iterations_times, filename)

@inline prettytime(t) = @sprintf "%06.1f" t

simulation_parameters(file) = file["metadata/parameters"]
simulation_parameters(filename::String) = jldopen(simulation_parameters, filename)

@inline function variable_strain_rate(t, sp)
    turnon = max(1-exp(-sp.f * t / 15), 0)
    return sp.α * turnon
end
# -------------------------------------------------------------

# -------------------------------------------------------------
include("filt.jl")
include("subfig_label.jl")
include("grid_nodes.jl")
include("bounds.jl")
include("get_field.jl")
include("timeseries_of.jl")
include("time_average_of.jl")
# -------------------------------------------------------------
