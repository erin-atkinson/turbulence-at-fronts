using Oceananigans
using JLD2

using Oceananigans.Fields: AbstractField, compute_at!
using Oceananigans.Utils: SumOfArrays
using Oceananigans: fill_halo_regions!

using Oceananigans.OutputWriters: saveproperty!, jld2output!

initial_time = Int(time_ns())
prev_time = Int(time_ns())

function update_field!(field, fieldtimeseries, frame)
    parent(field) .= parent(fieldtimeseries[frame])
    return nothing
end

function update_fields!(fields, fieldstimeseries, clock, frame)
    update_field!(fields.u, fieldstimeseries.u, frame)
    update_field!(fields.v, fieldstimeseries.v, frame)
    update_field!(fields.w, fieldstimeseries.w, frame)
    update_field!(fields.b, fieldstimeseries.b, frame)
    update_field!(fields.p, fieldstimeseries.p, frame)

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

fieldsymbols = (; u=:u, v=:v, w=:w, b=:b, p=:p)

fieldstimeseries = map(fieldsymbols) do ξ
    ξ_string = String(ξ)
    FieldTimeSeries(inputfilename, ξ_string; backend=OnDisk())
end

fields = map(x->x[1], fieldstimeseries)

# Initialise a clock
clock = Clock(; time=0.0)

grid = jldopen(read_grid, inputfilename)
iterations = jldopen(read_iterations, inputfilename)
times = jldopen(read_times, inputfilename)
sp = jldopen(read_parameters, parameterfilename)

include("terms/strainflow.jl")
fields = merge(fields, (; U, V, W))

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
=#
@info "Including $scriptname.jl"
include("$scriptname.jl")

# Write grid to file
jldopen(file->saveproperty!(file, "grid", grid), outputfilename, "a")
jldopen(file->saveproperty!(file, "grid", grid), tempfilename, "a")

function eltimestring()
    str = string(round((Int(time_ns()) - prev_time)/1e9; digits=3), " s")
    global prev_time = Int(time_ns())
    str
end

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
prev_time = Int(time_ns())
for (frame, iteration, time) in zip(frames, iterations, times)
    # Reset counter after giving kernel functions chance to compile
    frame == 11 && global setup_time = Int(time_ns())
    update_clock!(clock, iterations, times, frame)
    update_fields!(fields, fieldstimeseries, clock, frame)

    compute_fields_at!(dependency_fields, frame)
    
    write_outputs(outputfilename, iteration, time, output_fields)
    write_outputs(tempfilename, iteration, time, temp_fields)
    progstring = if frame > 10
        setupstr = round((setup_time - start_time)/1e9; digits=3)
        
        avg_time = (Int(time_ns()) - setup_time)/(1e9*(frame-10))
        avgstr = round(avg_time; digits=3)
        elapsed_time = round((Int(time_ns()) - initial_time)/1e9; digits=3)
        total_time = round((setup_time - initial_time)/1e9 + avg_time * (length(frames) - 10); digits=3)
        
        string("$frame of $(frames[end]), ", "$(eltimestring()), ", "setup: $(setupstr)s, ", "avg: $(avgstr)s, ", "est: $(elapsed_time)s / $(total_time)s")
    else
        setupstr = round((Int(time_ns()) - start_time)/1e9; digits=3)
        string("$frame of $(frames[end]), ", "$(eltimestring()), ", "setup: $(setupstr)s")
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
