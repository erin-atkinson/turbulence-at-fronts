ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0

using Oceananigans
include("create_simulation.jl")

#= 
    Create a simulation of a front forced by strain flow and surface cooling

    Call with 
        julia -julia_opts -- simulation.jl ARGS...

    ARGS is a set of input arguments:
        [01]: Path to output folder
        [02]: Simulation run time in seconds
        [03]: Coriolis parameter in 1/seconds
        [04]: Height of the mixed layer in metres
        [05]: Number of across-front (x) grid points
        [06]: Number of along-front (y) grid points
        [07]: Number of vertical (z) grid points
        [08]: Maximum Rossby number of front
        [09]: Minimum Richardson number of front
        [10]: Strain rate in 1/seconds
        [11]: Ekman number
        [12]: Surface buoyancy loss (positive for cooling) in watts / sq. metre
        [13]: Sponge layer damping multiplier
        [14]: Sponge layer damping size relative to front
        [15]: Cooling spin-up duration

    The simulation is run with only surface cooling to generate initial fluctuation fields,
    then these are combined with a mean reference front generated to satisfy the given non-dimensional parameters.
    
    A sponge layer is added around the front to prevent the strain flow from causing numerical instabilities, and
    to prevent internal wave reflections from the frontogenesis.
=#
output_folder = ARGS[1]
run_time = Meta.parse(ARGS[2])

f = Meta.parse(ARGS[3])
H = Meta.parse(ARGS[4])
Nx = Meta.parse(ARGS[5])
Ny = Meta.parse(ARGS[6])
Nz = Meta.parse(ARGS[7])
Ro = Meta.parse(ARGS[8])
Ri = Meta.parse(ARGS[9])
α = Meta.parse(ARGS[10])
Ek = Meta.parse(ARGS[11])
Q = Meta.parse(ARGS[12])
c = Meta.parse(ARGS[13])
damping_frac = Meta.parse(ARGS[14])
turbulence_spinup = Meta.parse(ARGS[15])

simulation_parameters = (; f, H, Nx, Ny, Nz, Ro, Ri, α, Q, c, damping_frac, turbulence_spinup, Ek)

simulation = create_simulation(run_time, output_folder, simulation_parameters)
run!(simulation)
