using Oceananigans
using JLD2
using Printf

using Oceananigans.Fields: compute_at!
using Oceananigans.OutputWriters: saveproperty!

include("update_clock.jl")
include("update_fields.jl")

t0 = time()

# Folder containing simulation output
foldername = ARGS[1]

# Defining output terms
scriptname = ARGS[2]

# Possible third argument is a temporary location
buffer = length(ARGS) > 2 ? ARGS[3] : ARGS[1]

# Path to input data and output
RAW = joinpath(foldername, "output.jld2")
BUFFER = joinpath(buffer, "$scriptname.jld2")

# Path to temp output
TEMP = joinpath(buffer, "temp_$scriptname.jld2")

# Final output folder
PROCESSED = joinpath(foldername, "$scriptname.jld2")

# Read simulation state
fds = FieldDataset(RAW; backend=OnDisk())

if "parameters" ∉ keys(fds.metadata)
    parameterfilename = joinpath(foldername, "parameters.jld2")
    fds.metadata["parameters"] = jldopen(file->file["simulation"], parameterfilename)
end

fieldnames = Symbol.(keys(fds.fields))
const sp = fds.metadata["parameters"]

# Setup grid, times and iterations
grid = fds.u.grid
times = fds.u.times
iterations = jldopen(file->keys(file["timeseries/t"]), RAW)
iterations = parse.(Int, iterations)
frames = 1:length(iterations)

# Named tuple of current simulation state fields
rawfields = NamedTuple(k => fds[k][1] for k in fieldnames)
nextrawfields = NamedTuple(Symbol(k, :_next) => fds[k][2] for k in fieldnames)

# Setup background strain
include("../terms/strainflow.jl")
input_fields = merge(rawfields, nextrawfields, (; U, V, W))

# Initialise a clock
clock = Clock(; time=times[1])

#= 
Input Julia file should define some things:
    `dependency_fields`:
        NamedTuple of fields to call compute! on
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
dependency_fields = NamedTuple()
temp_fields = NamedTuple()
skip_update = ()
cleanup() = nothing
@info "Including $scriptname.jl"
include("../$scriptname.jl")

output_fds = FieldDataset(times, output_fields; 
    backend = OnDisk(), 
    path = BUFFER,
    metadata = fds.metadata
)

temp_fds = if length(temp_fields) > 0
    FieldDataset(times, temp_fields; 
        backend = OnDisk(), 
        path = TEMP,
        metadata = fds.metadata
    )
else
    nothing
end

# Helpful for debugging to print times for all
@info "Performing first computation..."
print("Updating clock...\r")
dt = @elapsed update_clock!(clock, iterations, times, frames[1])
@printf "Updated clock! Elapsed: %.2f\n" dt

print("Updating fields...\r")
dt = @elapsed update_fields!(input_fields, fds, clock, frames[1]; skip_update)
@printf "Updated fields! Elapsed: %.2f\n" dt

for (k, dependency_field) in pairs(dependency_fields)
    print("Calculating $k...\r")
    local dt = @elapsed compute_at!(dependency_field, frames[2])
    @printf "Calculated %s! Elapsed: %.2f\n" k dt
end

@info "Computing..."
t1 = time()
t2 = t1
for (i, frame) in enumerate(frames)
    iteration = iterations[frame]
    t = times[frame]

    # Update clock
    update_clock!(clock, iterations, times, frame)

    # Update inputs
    update_fields!(input_fields, fds, clock, frame; skip_update)

    for field in dependency_fields
        compute_at!(field, t)
    end
    
    set!(output_fds, iteration, t; output_fields...)
    temp_fds != nothing && set!(temp_fds, iteration, t; temp_fields...)

    # Little bit of timekeeping for convenience
    tstr = if i < 11
        global t2 = time()
        setup_time = t2 - t0
        tstr = @sprintf "Setup: %.2f s" setup_time
    else
        t3 = time()
        setup_time = t2 - t0
        avg_time = (t3 - t2) / (i - 10)
        total_time = setup_time + avg_time * (length(frames) - 10)
        @sprintf "Setup: %.2f s, Avg: %.2f s, Total: %.2f s" setup_time avg_time total_time
    end
    print("$(frames[1]) -> $frame -> $(frames[end]) | $tstr\r")
end
println()
@info "Cleaning up..."
cleanup()

# Write grid to file
jldopen(file->saveproperty!(file, "grid", grid), BUFFER, "a")

if !isequal(BUFFER, PROCESSED)
    @info "Moving from $BUFFER to $PROCESSED"
    mv(BUFFER, PROCESSED; force=true)
end
rm(TEMP; force=true)



@info "Done!"