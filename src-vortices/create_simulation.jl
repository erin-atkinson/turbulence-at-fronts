using Oceananigans
using JLD2

include("parameters.jl")
include("forcing_functions.jl")
include("boundary_conditions.jl")
include("grid_cells.jl")
include("base_state.jl")
include("closure.jl")

@inline function create_simulation(run_time, output_folder, simulation_parameters)
    
    sp = create_simulation_parameters(simulation_parameters)
    
    base_state = get_base_state(sp)
    (xs, ys, zs) = get_grid_faces(sp)
    forcing = create_forcings(base_state, sp)
    closure = create_closure(xs, ys, zs, base_state, sp)
    boundary_conditions = create_boundary_conditions(base_state, sp)
    
    grid = RectilinearGrid(GPU();
        x=xs,
        y=ys,
        z=zs,
        size=(sp.Nx, sp.Ny, sp.Nz),
        topology=(Bounded, Bounded, Bounded)
    )
    @info grid
    model = NonhydrostaticModel(; grid,
        coriolis = FPlane(; sp.f),
        closure,
        tracers = (:b, ),
        buoyancy = BuoyancyTracer(),
        forcing,
        boundary_conditions
    )
    @info model
    set!(model; base_state...)
    Δt = 1e-3/sp.f
    simulation = Simulation(model; Δt, stop_time=run_time)
    
    wizard = TimeStepWizard(cfl=0.2, diffusive_cfl=0.2)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
    u, v, w = model.velocities
    b = model.tracers.b
    pHY′ = model.pressures.pHY′
    pNHS = model.pressures.pNHS
    #ν = model.diffusivity_fields[2].νₑ
    if !isdir("$output_folder")
        mkdir(output_folder)
    end
    simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w, b, pNHS, pHY′); 
        filename="$output_folder/output.jld2", 
        schedule=TimeInterval(1e-1/sp.f),
        overwrite_existing=true
    )
    progress(sim) = print("$(round(100*time(sim) / run_time; digits=1))%\r")
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(500))
    jldopen("$output_folder/parameters.jld2", "w") do file
        for key in keys(sp)
            file[string(key)] = sp[key]
        end
    end
    @info simulation
    simulation
end


















