# postprocess.jl
# Boiler-plate code for offline processing

using Oceananigans
using JLD2

using Oceananigans.Fields: AbstractField, compute_at!
using Oceananigans.Utils: SumOfArrays
using Oceananigans: fill_halo_regions!

using Oceananigans.OutputWriters: saveproperty!, jld2output!

initial_time = Int(time_ns())
prev_time = Int(time_ns())

function timestring(t_ns; digits=3)
    return string(round(1e-9t_ns; digits), " s")
end

function eltimestring()
    t_ns = Int(time_ns())
    str = timestring(t_ns - prev_time)
    global prev_time = Int(t_ns)
    return str
end

function update_field!(field, fieldtimeseries, frame)
    parent(field) .= parent(fieldtimeseries[frame])
    return nothing
end

function update_fields!(fields, fieldstimeseries, clock, frame)
    update_field!(fields.u, fieldstimeseries.u, frame)
    update_field!(fields.v, fieldstimeseries.v, frame)
    update_field!(fields.w, fieldstimeseries.w, frame)
    update_field!(fields.b, fieldstimeseries.b, frame)
    update_field!(fields.pNHS, fieldstimeseries.pNHS, frame)

    l = length(fieldstimeseries.u)
    update_field!(fields.u_next, fieldstimeseries.u, min(frame + 1, l))
    update_field!(fields.v_next, fieldstimeseries.v, min(frame + 1, l))
    update_field!(fields.w_next, fieldstimeseries.w, min(frame + 1, l))
    update_field!(fields.b_next, fieldstimeseries.b, min(frame + 1, l))
    update_field!(fields.pNHS_next, fieldstimeseries.pNHS, min(frame + 1, l))

    compute_background!(fields.U, fields.V, fields.W, clock)
    return nothing
end

function update_clock!(clock, iterations, times, frame)
    clock.time = times[frame]
    clock.iteration = parse(Int, iterations[frame])
    clock.last_Δt = frame > 1 ? times[frame] - times[frame-1] : Inf
    return nothing
end

function compute_fields_at!(dependency_fields, frame)
    compute_at!(dependency_fields, frame)
    return nothing
end

function write_outputs(filename, iteration, time, outputs)
    data = map(parent, outputs)
    jld2output!(filename, iteration, time, data, (; ))
    return nothing
end
cleanup() = nothing
temp_fields = (; )

read_grid(file) = file["serialized/grid"]
read_iterations(file) = keys(file["timeseries/t"])
read_times(file) = map(iteration -> file["timeseries/t/$iteration"], read_iterations(file))
read_parameters(file) = file["simulation"]

foldername = ARGS[1]
scriptname = ARGS[2]

# Possible third argument is a temporary location
buffer = length(ARGS) > 2 ? ARGS[3] : ARGS[1]

# Filenames
inputfilename = joinpath(foldername, "output.jld2")
outputfilename = joinpath(buffer, "$scriptname.jld2")
tempfilename = joinpath(buffer, "temp_$scriptname.jld2")
parameterfilename = joinpath(foldername, "parameters.jld2")

finalfilename = joinpath(foldername, "$scriptname.jld2")

@info "Reading timeseries from file"

fieldsymbols = (; u=:u, v=:v, w=:w, b=:b, pNHS=:pNHS)

fieldstimeseries = map(fieldsymbols) do ξ
    ξ_string = String(ξ)
    FieldTimeSeries(inputfilename, ξ_string; backend=OnDisk())
end

fields = map(x->x[1], fieldstimeseries)

next_fields = (; 
    u_next=fieldstimeseries.u[2],
    v_next=fieldstimeseries.v[2],
    w_next=fieldstimeseries.w[2],
    b_next=fieldstimeseries.b[2],
    pNHS_next=fieldstimeseries.pNHS[2],
)

# Initialise a clock
clock = Clock(; time=0.0)

grid = jldopen(read_grid, inputfilename)
iterations = jldopen(read_iterations, inputfilename)
times = jldopen(read_times, inputfilename)
sp = jldopen(read_parameters, parameterfilename)

include("terms/strainflow.jl")
fields = merge(fields, next_fields, (; U, V, W))

frames = 1:length(iterations)

#= 
Input Julia file should define some things:
    `dependency_fields`:
        Ordered collection of fields to call compute! on
    `output_fields`:
        NamedTuple of fields that will get saved. Note that these should also be in
        `calculated_fields` if they need to be computed
    `temp_fields`:
        List of fields that will get saved temporarily. Note that these should also be in
        `calculated_fields` if they need to be computed.
    `cleanup`:
        Function to be called before temp_fields is deleted
In addition, it can redefine `frames` to process only a subset of frames
=#

@info "Including $scriptname.jl"
include("$scriptname.jl")

# Write grid to file
jldopen(file->saveproperty!(file, "grid", grid), outputfilename, "a")
jldopen(file->saveproperty!(file, "grid", grid), tempfilename, "a")

@info "Finished setup! Elapsed: $(eltimestring())"

@info "Performing first computation..."
print("Updating clock...\r")
update_clock!(clock, iterations, times, frames[1])
println("Updated clock! Elapsed: $(eltimestring())")

print("Updating fields...\r")
update_fields!(fields, fieldstimeseries, clock, frames[1])
println("Updated fields! Elapsed: $(eltimestring())")

map(keys(dependency_fields), dependency_fields) do k, dependency_field
    print("Calculating $k...\r")
    compute_at!(dependency_field, frames[1])
    println("Calculated $(k)! Elapsed: $(eltimestring())")
end

start_time = Int(time_ns())
eltimestring()
for (i, frame, iteration, time) in zip(1:length(frames), frames, iterations[frames], times[frames])
    # Reset counter after giving kernel functions chance to compile for better estimate
    i == 11 && global setup_time = Int(time_ns())
    update_clock!(clock, iterations, times, frame)
    update_fields!(fields, fieldstimeseries, clock, frame)

    compute_fields_at!(dependency_fields, frame)
    
    write_outputs(outputfilename, iteration, time, output_fields)
    write_outputs(tempfilename, iteration, time, temp_fields)
    progstring = if i > 10
        setupstr = timestring(setup_time - start_time)
        
        avg_time = (Int(time_ns()) - setup_time)/(i-10)
        avgstr = timestring(avg_time)
        elapsed_time = timestring(Int(time_ns()) - initial_time)
        total_time = timestring(setup_time - initial_time + avg_time * (length(frames) - 10))
        
        string("$frame of $(frames[end]), ", "$(eltimestring()), ", "setup: $(setupstr), ", "avg: $(avgstr), ", "est: $(elapsed_time) / $(total_time)")
    else
        setupstr = timestring(Int(time_ns()) - start_time)
        string("$frame of $(frames[end]), ", "$(eltimestring()), ", "setup: $(setupstr)")
    end
    print(rpad(progstring, 80, ' '), "\r")
end
println()
cleanup()

if !isequal(outputfilename, finalfilename)
    @info "Moving from $buffer to $foldername"
    mv(outputfilename, finalfilename; force=true)
end
rm(tempfilename)

@info "Finished!"
