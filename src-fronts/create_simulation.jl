using Oceananigans
using JLD2

include("parameters.jl")
include("forcing_functions.jl")
include("boundary_conditions.jl")
include("grid_cells.jl")
include("base_state.jl")
include("closure.jl")
include("initialise.jl")

@inline function create_simulation(run_time, output_folder, simulation_parameters)
    
    sp = create_simulation_parameters(simulation_parameters)
    
    # Save the parameters to a file
    if !isdir("$output_folder")
        mkdir(output_folder)
    end
    jldopen("$output_folder/parameters.jld2", "w") do file
        for key in keys(sp)
            file[string(key)] = sp[key]
        end
    end
    
    # Run a simulation to create the initial state
    # Write
    @info "Generating an initialised state..."
    init_state = initialise_simulation(sp; init_write_freq=10*sp.f)
    
    # For debugging initialisation just return here
    if run_time == 0
        @warn "Returning no simulation, will LoadError later"
        return nothing
    end
    
    base_state = get_base_state(sp)
    (xs, ys, zs) = get_grid_faces(sp)
    
    grid = RectilinearGrid(GPU();
        x=xs,
        y=ys,
        z=zs,
        size=(sp.Nx, sp.Ny, sp.Nz),
        topology=(Bounded, Periodic, Bounded)
    )
    
    forcing = create_forcings(base_state, sp)
    closure = create_closure(xs, ys, zs, base_state, sp)
    boundary_conditions = create_boundary_conditions(base_state, sp)
    
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
    set!(model; init_state...)
    
    Δt = 1e-3/sp.f
    simulation = Simulation(model; Δt, stop_time=sp.turbulence_spinup)
    
    wizard = TimeStepWizard(cfl=0.2, diffusive_cfl=0.2)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
    
    # Save pressures and velocities and tracers
    u, v, w = model.velocities
    b = model.tracers.b
    #pHY′ = model.pressures.pHY′
    #pNHS = model.pressures.pNHS
    φ = model.pressures.pHY′ + model.pressures.pNHS
    νₑ = model.diffusivity_fields.νₑ
    κₑ = model.diffusivity_fields.κₑ
    simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w, b, φ, νₑ, κₑ); 
        filename="$output_folder/output.jld2", 
        schedule=TimeInterval(sp.write_freq),
        overwrite_existing=true,
        with_halos=true
    )
    
    # Output progress
    progress(sim) = print("Running simulation $(round(100*time(sim) / run_time; digits=1))%\r")
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(500))
    
    @info simulation
    simulation
end


















