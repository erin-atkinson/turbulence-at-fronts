#= simulation.jl
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
        [11]: Surface buoyancy loss (positive for cooling) in watts / sq. metre
        [12]: Sponge layer damping multiplier
        [13]: Start time of simulation (negative for a cooling pre-initialisation)
        [14]: Save interval
        [15]: Stretch factor
=#
ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0
using Oceananigans

output_folder = ARGS[1]

run_time = parse(Float64, ARGS[2])
f = parse(Float64, ARGS[3])
H = parse(Float64, ARGS[4])
Nx = parse(Int64, ARGS[5])
Ny = parse(Int64, ARGS[6])
Nz = parse(Int64, ARGS[7])
Ro = parse(Float64, ARGS[8])
Ri = parse(Float64, ARGS[9])
α = parse(Float64, ARGS[10])
Q = parse(Float64, ARGS[11])
c = parse(Float64, ARGS[12])
start_time = parse(Float64, ARGS[13])
save_time = parse(Float64, ARGS[14])
s = parse(Float64, ARGS[15])
comment = ARGS[16]

simulation_parameters = (; run_time, f, H, Nx, Ny, Nz, Ro, Ri, α, Q, c, start_time, save_time, s)

include("create_simulation.jl")
