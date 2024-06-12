using Oceananigans
using JLD2
using CUDA
using Statistics

include("parameters.jl")
include("forcing_functions.jl")
include("boundary_conditions.jl")
include("grid_cells.jl")
include("base_state.jl")
include("closure.jl")

@inline function initialise_simulation(simulation_parameters; init_write_freq=0)
    
    # To initialise first spin up the strain then the turbulence?
    
    # Spin up the turbulence
    sp = (simulation_parameters..., α=0)
    
    base_state = get_base_state(sp)
    (xs, ys, zs) = get_grid_faces(sp)
    
    # For the sake of the initialisation, the state should be with no front (at infinity)
    pre_init_state = map(base_state) do f
        (x, y, z) -> f(x, y, z) + 1e-8 * rand()
    end
    
    grid = RectilinearGrid(GPU();
        x=xs,
        y=ys,
        z=zs,
        size=(sp.Nx, sp.Ny, sp.Nz),
        topology=(Bounded, Periodic, Bounded)
    )
    
    forcing = create_forcings(base_state, sp)
    # Closure as usual
    closure = create_closure(xs, ys, zs, base_state, sp)
    # Boundary conditions as usual
    boundary_conditions = create_boundary_conditions(base_state, sp)
    
    @info grid
    model = NonhydrostaticModel(; grid,
        coriolis = FPlane(; sp.f),
        closure,
        tracers = (:b, ),
        buoyancy = BuoyancyTracer(),
        boundary_conditions,
        forcing
    )
    @info model
    
    # Set the model state to the pre-initialisation state
    set!(model; pre_init_state...)
    
    Δt = 1e-3/sp.f
    simulation = Simulation(model; Δt, stop_time=sp.turbulence_spinup)
    
    wizard = TimeStepWizard(cfl=0.2, diffusive_cfl=0.2)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
    
    u, v, w = model.velocities
    b = model.tracers.b
    φ = model.pressures.pHY′ + model.pressures.pNHS
    
    if !isdir("$output_folder")
        mkdir(output_folder)
    end
    # Write at least the final state to file
    simulation.output_writers[:initialisation] = if init_write_freq == 0
        JLD2OutputWriter(model, (; u, v, w, b); 
            filename="$output_folder/initialisation.jld2", 
            schedule=SpecifiedTimes([sp.turbulence_spinup]),
            overwrite_existing=true
        )
    else
        JLD2OutputWriter(model, (; u, v, w, b, φ); 
            filename="$output_folder/initialisation.jld2", 
            schedule=TimeInterval(1/init_write_freq),
            overwrite_existing=true
        )
    end
    progress(sim) = print("Running turbulent initialisation $(round(100*time(sim) / sp.turbulence_spinup; digits=1))%\r")
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(20))
    
    run!(simulation)
    
    # This adds the reference state but since TTW frontogenesis doesn't exist
    # it's better to start with the whole front and cool it and wait for
    # transients to die
    xsᶠᶜᶜ = Array(xnodes(grid, Face(), Center(), Center()))
    zsᶜᶜᶠ = Array(znodes(grid, Center(), Center(), Face()))
    xsᶜᶜᶜ = Array(xnodes(grid, Center(), Center(), Center()))
    zsᶜᶜᶜ = Array(znodes(grid, Center(), Center(), Center()))
    
    # Return the initialised fields, with the front added
    init_state = CUDA.@allowscalar (;
        u = Array(interior(model.velocities.u)), 
        v = Array(interior(model.velocities.v)),
        w = Array(interior(model.velocities.w)),
        b = Array(interior(model.tracers.b))
        )
    
    # Remove the mean
    init_state = map(init_state) do init_field
        init_field .-= mean(init_field; dims=2)
    end
    #@info  "init_state sizes = $(map(size, init_state))"
    #@info  "grid sizes = $(size(xsᶠᶜᶜ)), $(size(zsᶜᶜᶠ)), $(size(xsᶜᶜᶜ)), $(size(zsᶜᶜᶜ))"
    # Add the reference state
    init_state.u .+= [base_state.u(x, 0, z) for x in xsᶠᶜᶜ, y in [1], z in zsᶜᶜᶜ]
    init_state.v .+= [base_state.v(x, 0, z) for x in xsᶜᶜᶜ, y in [1], z in zsᶜᶜᶜ]
    init_state.w .+= [base_state.w(x, 0, z) for x in xsᶜᶜᶜ, y in [1], z in zsᶜᶜᶠ]
    init_state.b .+= [base_state.b(x, 0, z) for x in xsᶜᶜᶜ, y in [1], z in zsᶜᶜᶜ]
    return init_state
end


















