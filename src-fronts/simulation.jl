ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0

using Oceananigans


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
        [06]: Number of across-front grid points in high resolution area
        [07]: Number of along-front (y) grid points
        [08]: Number of vertical (z) grid points
        [09]: Maximum Rossby number of front
        [10]: Minimum Richardson number of front
        [11]: Strain rate in 1/seconds
        [13]: Surface buoyancy loss (positive for cooling) in watts / sq. metre
        [14]: Sponge layer damping multiplier
        [15]: Sponge layer damping size relative to front

    The simulation is run with only surface cooling to generate initial fluctuation fields,
    then these are combined with a mean reference front generated to satisfy the given non-dimensional parameters.
    
    A sponge layer is added around the front to prevent the strain flow from causing numerical instabilities, and
    to prevent internal wave reflections from the frontogenesis.
=#
output_folder = ARGS[1]

run_time = parse(Float64, ARGS[2])
f = parse(Float64, ARGS[3])
H = parse(Float64, ARGS[4])
Nx = parse(Int64, ARGS[5])
Nₕ = parse(Int64, ARGS[6])
Ny = parse(Int64, ARGS[7])
Nz = parse(Int64, ARGS[8])
Ro = parse(Float64, ARGS[9])
Ri = parse(Float64, ARGS[10])
α = parse(Float64, ARGS[11])
Q = parse(Float64, ARGS[12])
c = parse(Float64, ARGS[13])
damping_frac = parse(Float64, ARGS[14])
start_time = parse(Float64, ARGS[15])
save_time = parse(Float64, ARGS[16])

simulation_parameters = (; run_time, f, H, Nx, Nₕ, Ny, Nz, Ro, Ri, α, Q, c, damping_frac, start_time, save_time)

include("create_simulation.jl")
