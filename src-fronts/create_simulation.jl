using Oceananigans
using JLD2

include("parameters.jl")
include("forcing_functions.jl")
include("boundary_conditions.jl")
include("grid_cells.jl")
include("base_state.jl")
include("closure.jl")
include("initialise.jl")

@inline function create_simulation(output_folder, simulation_parameters)
    
    sp = create_simulation_parameters(simulation_parameters)
    
    # Save the parameters to a file
    if !isdir("$output_folder")
        mkdir(output_folder)
    end
    jldopen("$output_folder/parameters.jld2", "w") do file
        file["simulation"] = sp
    end
    
    (xs, ys, zs) = get_grid_faces(sp)
    
    grid = RectilinearGrid(GPU();
        x=xs,
        y=ys,
        z=zs,
        size=(sp.Nx, sp.Ny, sp.Nz),
        topology=(Bounded, Periodic, Bounded)
    )
    @info grid
    
    # Avoid indexing in GPU...
    base_state = let CPU_grid = RectilinearGrid(CPU();
            x=xs,
            y=ys,
            z=zs,
            size=(sp.Nx, sp.Ny, sp.Nz),
            topology=(Bounded, Periodic, Bounded)
        )
        get_base_state(CPU_grid, sp)
    end
    #base_state = get_base_state(grid, sp)
    
    init_state = (;
        u=base_state.u .+ 1e-5 * (rand(size(base_state.u)...) .- 0.5),
        w=base_state.w .+ 1e-6 * (rand(size(base_state.w)...) .- 0.5),
        base_state.v,
        base_state.b
    )
    
    forcing = create_forcings(base_state, sp)
    @info forcing
    
    closure = create_closure(xs, ys, zs, base_state, sp)
    @info closure
    
    boundary_conditions = create_boundary_conditions(base_state, sp)
    @info boundary_conditions
    
    
    model = NonhydrostaticModel(; grid,
        clock=Clock(time=sp.start_time),
        coriolis = FPlane(; sp.f),
        closure,
        tracers = (:b,),
        buoyancy = BuoyancyTracer(),
        forcing,
        boundary_conditions,
        hydrostatic_pressure_anomaly=CenterField(grid)
    )
    @info model
    set!(model; init_state...)
    
    Δt = 1e-3/sp.f
    
    # 
    checkpoint_files = filter(readdir(output_folder)) do x
        occursin(r"^checkpoint", x)
    end
    
    # Take the latest checkpoint file (there should only be one though...)
    prev_time = mapreduce(max, checkpoint_files; init=sp.start_time) do checkpoint_file
        jldopen(file->file["clock"].time, "$output_folder/$checkpoint_file")
    end
    
    simulation = Simulation(model; Δt, stop_time=prev_time + sp.run_time)
    
    wizard = TimeStepWizard(cfl=0.2, diffusive_cfl=0.2)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
    
    # Save pressures and velocities and tracers
    u, v, w = model.velocities
    b = model.tracers.b
    ϕ = model.pressures.pHY′ + model.pressures.pNHS
    ν = model.diffusivity_fields.νₑ
    
    # This is a workaround for a quirk with loading a checkpoint
    writing_times = prev_time:0.1/sp.f:(prev_time+sp.run_time)
    
    simulation.output_writers[:fields] = JLD2OutputWriter(model, (; u, v, w, b, ϕ, ν); 
        filename="$output_folder/output.jld2", 
        schedule=SpecifiedTimes(writing_times),
        overwrite_existing=false,
        with_halos=true
    )
    
    # Add a checkpointer that saves at every 10%
    simulation.output_writers[:checkpointer] = Checkpointer(model;
        schedule=SpecifiedTimes(range(prev_time, prev_time + sp.run_time, 10)),
        dir=output_folder,
        cleanup=true,
        overwrite_existing=false,
        verbose=true
    )
    
    # Output progress
    progress(sim) = print("Running simulation $(round(100*time(sim) / sim.stop_time; digits=2))($(round(100*(time(sim)-prev_time) / sp.run_time; digits=2)))%\r")
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))
    
    @info simulation
    run!(simulation; pickup=true)
end


















