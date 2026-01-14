using Oceananigans
using CUDA
using JLD2
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Printf

include("base_state.jl")
include("grid_faces.jl")
include("parameters.jl")

const sp = create_simulation_parameters(simulation_parameters)

function init_jld2!(file, comment, parameters)
    file["metadata/author"] = "Erin Atkinson"
    file["metadata/comment"] = comment
    file["metadata/parameters"] = parameters
    return nothing
end

if !isdir("$output_folder")
    mkdir(output_folder)
end

# Get the grid
grid_faces = get_grid_faces(sp)
@info "Created grid faces"
(xs, ys, zs) = grid_faces

grid = RectilinearGrid(GPU();
    x=xs,
    y=ys,
    z=zs,
    size=(sp.Nx, sp.Ny, sp.Nz),
    topology=(Bounded, Periodic, Bounded),
    halo=(5, 5, 5)
)
@info grid

init_state = front_initial_conditions(grid, sp)

# Forcing functions
include("forcing.jl")

# Boundary conditions
include("boundary_conditions.jl")

# Closure
include("closure.jl")

model = NonhydrostaticModel(; grid,
    clock=Clock(time=sp.start_time),
    advection=WENO(; order=9),
    coriolis = FPlane(; sp.f),
    tracers = (:b, ),
    closure,
    buoyancy = BuoyancyTracer(),
    forcing,
    boundary_conditions
)

@info model
set!(model; init_state...)

# Some initial timestep...
Δt = 1e-3 * sp.save_time

checkpoint_files = filter(readdir(output_folder)) do x
    occursin(r"^checkpoint", x)
end

# Take the latest checkpoint file
prev_time = mapreduce(max, checkpoint_files; init=sp.start_time * 1.0) do checkpoint_file
    jldopen(file->file["NonhydrostaticModel/clock"].time, joinpath(output_folder, checkpoint_file))
end

simulation = Simulation(model; Δt, stop_time=prev_time+sp.run_time)

# Save pressure anomaly and velocities and tracers
u, v, w = model.velocities
b = model.tracers.b
pNHS = model.pressures.pNHS

# This is a workaround for a quirk with loading a checkpoint
writing_times = prev_time:sp.save_time:(prev_time+sp.run_time)

simulation.output_writers[:fields] = JLD2Writer(model, (; u, v, w, b, pNHS); 
    filename="$output_folder/output.jld2", 
    schedule=SpecifiedTimes(writing_times),
    overwrite_existing=false,
    with_halos=true,
    init=(file, model)->init_jld2!(file, comment, sp)
)

# Add a checkpointer
simulation.output_writers[:checkpointer] = Checkpointer(model;
    schedule=SpecifiedTimes(writing_times[end]),
    dir=output_folder,
    cleanup=true,
    overwrite_existing=false,
    verbose=true
)

# Variable time step
wizard = TimeStepWizard(; cfl=0.5, max_Δt=1e-3/sp.f)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Compute the advection
simulation.callbacks[:advection] = Callback(calculate_U_callback, IterationInterval(1); parameters=sp)

# Output progress
const prev_t = [0.0]
const prev_wall_time = [0.0]
function progress(sim)
    i = iteration(sim)

    t = time(sim)
    Δ_time = t - prev_t[1]
    prev_t[1] = t
    t_str = @sprintf " -- Time: %.3e" t

    wall_time = sim.run_wall_time
    Δ_wall_time = wall_time - prev_wall_time[1]
    prev_wall_time[1] = wall_time

    t_per_hour = Δ_time / (Δ_wall_time / 3600)
    tph_str = @sprintf " -- Time / wall hour: %.3e" t_per_hour

    remaining_time = wall_time / 3600 + (sim.stop_time - t) / t_per_hour
    remaining_str = @sprintf "%.1f hr / %.1f hr" (wall_time / 3600) remaining_time

    str = string("Iteration: ", i, t_str, tph_str, " -- Progress: ", remaining_str)
    
    print(rpad("\r$str", 100))

    return nothing
end
simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

@info simulation
run!(simulation; pickup=true)
