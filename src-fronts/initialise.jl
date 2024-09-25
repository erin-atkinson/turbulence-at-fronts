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
    
    # Spin up the turbulence
    sp = (simulation_parameters..., α=0)
    
    (xs, ys, zs) = get_grid_faces(sp)
    
    grid = RectilinearGrid(GPU();
        x=xs,
        y=ys,
        z=zs,
        size=(sp.Nx, sp.Ny, sp.Nz),
        topology=(Bounded, Periodic, Bounded)
    )
    base_state = get_base_state(grid, sp)
    # For the sake of the initialisation, the state should be with no front (at infinity)
    pre_init_state = (;
        u=base_state.u .+ 1e-5 * (rand(size(base_state.u)...) .- 0.5),
        w=base_state.w .+ 1e-6 * (rand(size(base_state.w)...) .- 0.5),
        base_state.v,
        base_state.b
    )
    
    if sp.turbulence_spinup == 0
        @info "No turbulence spinup period specified, skipping init"
        return pre_init_state
    end
    
    forcing = create_forcings(base_state, sp)
    # Closure as usual
    closure = create_closure(xs, ys, zs, base_state, sp)
    # Boundary conditions as usual
    boundary_conditions = create_boundary_conditions(base_state, sp)
    Δt = 1e-3/sp.f
    @info grid
    model = NonhydrostaticModel(; grid,
        coriolis = FPlane(; sp.f),
        closure,
        tracers = (:b, ),
        buoyancy = BuoyancyTracer(),
        boundary_conditions,
        forcing, 
        hydrostatic_pressure_anomaly=CenterField(grid)
    )
    @info model
    # Set the model state to the pre-initialisation state
    set!(model; pre_init_state...)
    
    simulation = Simulation(model; Δt, stop_time=sp.turbulence_spinup)
    wizard = TimeStepWizard(cfl=0.2, diffusive_cfl=0.2)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
    
    u, v, w = model.velocities
    b = model.tracers.b
    φ = model.pressures.pNHS + model.pressures.pHY′
    νₑ = sp.Ek == 0 ? model.diffusivity_fields.νₑ : model.diffusivity_fields[2].νₑ
    
    if !isdir("$output_folder")
        mkdir(output_folder)
    end
    # Write at least the final state to file
    simulation.output_writers[:initialisation] = if init_write_freq == 0
        JLD2OutputWriter(model, (; u, v, w, b); 
            filename="$output_folder/initialisation.jld2", 
            schedule=SpecifiedTimes([sp.turbulence_spinup]),
            overwrite_existing=true,
            with_halos=true
        )
    else
        JLD2OutputWriter(model, (; u, v, w, b, φ, νₑ); 
            filename="$output_folder/initialisation.jld2", 
            schedule=TimeInterval(1/init_write_freq),
            overwrite_existing=true,
            with_halos=true
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
    
    # Return the initialised fields
    init_state = CUDA.@allowscalar (;
        u = Array(interior(model.velocities.u)), 
        v = Array(interior(model.velocities.v)),
        w = Array(interior(model.velocities.w)),
        b = Array(interior(model.tracers.b))
        )
    #=
    # Remove the mean
    init_state = map(init_state) do init_field
        init_field .-= mean(init_field; dims=2)
    end
    #@info  "init_state sizes = $(map(size, init_state))"
    #@info  "grid sizes = $(size(xsᶠᶜᶜ)), $(size(zsᶜᶜᶠ)), $(size(xsᶜᶜᶜ)), $(size(zsᶜᶜᶜ))"
    # Add the reference state
    init_state.u .+= base_state.u
    init_state.v .+= base_state.v
    init_state.w .+= base_state.w
    init_state.b .+= base_state.b
    =#
    return init_state
end


















