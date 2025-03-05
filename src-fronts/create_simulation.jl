using Oceananigans
using JLD2
using Oceananigans.BoundaryConditions: fill_halo_regions!

include("base_state.jl")
include("parameters.jl")
include("grid_cells.jl")

const sp = let a = create_simulation_parameters(simulation_parameters)
    (; a..., δ=0)
end

# Save the parameters to a file
if !isdir("$output_folder")
    mkdir(output_folder)
end
jldopen("$output_folder/parameters.jld2", "w") do file
    file["simulation"] = sp
end

# Get the grid
(xs, ys, zs) = get_grid_faces(sp)

grid = RectilinearGrid(GPU();
    x=xs,
    y=ys,
    z=zs,
    size=(sp.Nx, sp.Ny, sp.Nz),
    topology=(Bounded, Periodic, Bounded),
    halo=(5, 5, 5)
)
@info grid

# Avoid indexing in GPU by generating a CPU grid
base_state = let CPU_grid = RectilinearGrid(CPU();
        x=xs,
        y=ys,
        z=zs,
        size=(sp.Nx, sp.Ny, sp.Nz),
        topology=(Bounded, Periodic, Bounded)
    )
    get_base_state(CPU_grid, sp)
end

# Initialise with some random secondary circulation
init_state = (;
    u=base_state.u .+ 1e-6 * (rand(size(base_state.u)...) .- 0.5),
    w=base_state.w .+ 1e-7 * (rand(size(base_state.w)...) .- 0.5),
    base_state.v,
    base_state.b
)

#side_b_bcs = FieldBoundaryConditions(grid, (Nothing, Nothing, Center),
#    bottom=GradientBoundaryCondition(sp.N₀²),
#    top=FluxBoundaryCondition(sp.B)
#)

# These store the buoyancy profiles to relax to at the edges of the domain
b_east = Field{Nothing, Nothing, Center}(grid) #; boundary_conditions=side_b_bcs)
b_west = Field{Nothing, Nothing, Center}(grid) #; boundary_conditions=side_b_bcs)

# Forcing functions
include("forcing_functions.jl")

# Boundary conditions
include("boundary_conditions.jl")

# Closure
include("closure.jl")

model = NonhydrostaticModel(; grid,
    clock=Clock(time=sp.start_time),
    coriolis = FPlane(; sp.f),
    tracers = (:b, :ξ),
    closure,
    buoyancy = BuoyancyTracer(),
    forcing,
    boundary_conditions,
    hydrostatic_pressure_anomaly=CenterField(grid)
)
@info model
@inline ξ_init(x, y, z) = exp(-(x^2+y^2) / (sp.Ly / 10.0)^2)#exp(-(z+sp.Lz)/(sp.Lz - sp.H))
set!(model; init_state..., ξ=ξ_init)

b = model.tracers.b

b_west_op = KernelFunctionOperation{Nothing, Nothing, Center}(b_west_func, grid, b, sp)
b_east_op = KernelFunctionOperation{Nothing, Nothing, Center}(b_east_func, grid, b, sp)

b_west_st = Field(b_west_op)
b_east_st = Field(b_east_op)

compute!(b_west_st)
compute!(b_east_st)

b_west .= b_west_st
b_east .= b_east_st

# Callback that computes the b damping fields
function compute_side_b!(sim)
    
    b_west .= b_west_op
    b_east .= b_east_op
    #fill_halo_regions!(b_west)
    #fill_halo_regions!(b_east)
    return nothing
end

Δt = 1e-3 * sp.save_time

checkpoint_files = filter(readdir(output_folder)) do x
    occursin(r"^checkpoint", x)
end

# Take the latest checkpoint file (there should only be one though...)
prev_time = mapreduce(max, checkpoint_files; init=sp.start_time * 1.0) do checkpoint_file
    jldopen(file->file["clock"].time, "$output_folder/$checkpoint_file")
end

simulation = Simulation(model; Δt, stop_time=prev_time + sp.run_time)

wizard = TimeStepWizard(; cfl=0.5, diffusive_cfl=0.5, max_Δt=0.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Save pressures and velocities and tracers
u, v, w = model.velocities
b = model.tracers.b
p = model.pressures.pHY′ + model.pressures.pNHS
ξ = model.tracers.ξ
ν = model.diffusivity_fields.νₑ

# This is a workaround for a quirk with loading a checkpoint
writing_times = prev_time:sp.save_time:(prev_time+sp.run_time)

simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w, b, ξ, p, ν); 
    filename="$output_folder/output.jld2", 
    schedule=SpecifiedTimes(writing_times),
    overwrite_existing=false,
    with_halos=true
)

# Add a checkpointer that saves at every 10%
simulation.output_writers[:checkpointer] = Checkpointer(model;
    schedule=SpecifiedTimes(range(prev_time, prev_time + sp.run_time, 11)),
    dir=output_folder,
    cleanup=true,
    overwrite_existing=false,
    verbose=true
)

# Compute the forcing buoyancy

simulation.callbacks[:side_b] = Callback(compute_side_b!, IterationInterval(10))
progress(sim, p) = print("Running simulation $(round(time(sim); digits=2)), iterations per simulated second $(round(iteration(sim) / (time(sim) - p.prev_time); digits=4))\r")
# Output progress
simulation.callbacks[:progress] = Callback(progress, IterationInterval(50), parameters=(; prev_time, sp.run_time))

@info simulation
run!(simulation; pickup=true)
