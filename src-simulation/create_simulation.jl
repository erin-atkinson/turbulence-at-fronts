using Oceananigans
using JLD2
using Oceananigans.BoundaryConditions: fill_halo_regions!

include("base_state.jl")
include("parameters.jl")

const sp = create_simulation_parameters(simulation_parameters)

function write_comment!(file, comment)
    file["author"] = "Erin Atkinson"
    file["comment"] = comment
    return nothing
end

# Save the parameters to a file
if !isdir("$output_folder")
    mkdir(output_folder)
end
jldopen("$output_folder/parameters.jld2", "w") do file
    file["simulation"] = sp
    write_comment!(file, comment)
end

# Get the grid
include("grid_faces.jl")
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
    advection=WENO(grid; order=9),
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
Δt = 1e-1 * sp.save_time

checkpoint_files = filter(readdir(output_folder)) do x
    occursin(r"^checkpoint", x)
end

# Take the latest checkpoint file (there should only be one though...)
prev_time = mapreduce(max, checkpoint_files; init=sp.start_time * 1.0) do checkpoint_file
    jldopen(file->file["clock"].time, joinpath(output_folder, checkpoint_file))
end

simulation = Simulation(model; Δt, stop_time=prev_time + sp.run_time)

# Save pressure anomaly and velocities and tracers
u, v, w = model.velocities
b = model.tracers.b
pNHS = model.pressures.pNHS

# This is a workaround for a quirk with loading a checkpoint
writing_times = prev_time:sp.save_time:(prev_time+sp.run_time)

simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w, b, pNHS); 
    filename="$output_folder/output.jld2", 
    schedule=SpecifiedTimes(writing_times),
    overwrite_existing=false,
    with_halos=true,
    init=(file, model)->write_comment!(file, comment)
)

# Add a checkpointer that saves at every 10%
simulation.output_writers[:checkpointer] = Checkpointer(model;
    schedule=SpecifiedTimes(range(prev_time, prev_time + sp.run_time, 11)),
    dir=output_folder,
    cleanup=true,
    overwrite_existing=false,
    verbose=true
)

# Variable time step
wizard = TimeStepWizard(; cfl=0.5, max_Δt=1e-3/sp.f)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Compute the advection
simulation.callbacks[:advection] = Callback(calculate_UV_callback, IterationInterval(1); parameters=sp)

# Output progress
progress(sim) = print("Running simulation t=$(round(time(sim); digits=2)) iter=$(iteration(sim))\r")
simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))

@info simulation
run!(simulation; pickup=true)
