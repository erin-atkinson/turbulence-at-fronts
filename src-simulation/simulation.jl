#= simulation.jl
    Create a simulation of a front forced by strain flow and surface cooling

    Call with 
        julia -julia_opts -- simulation.jl ARGS...

    ARGS is a set of input arguments:
        [01]: Path to output folder

        [02]: Simulation run time in seconds
        [03]: Start time of simulation (negative for a cooling pre-initialisation)
        [04]: Save interval

        [05]: Coriolis parameter in 1/seconds

        [06]: Across-front domain size (relative to front)
        [07]: Central region domain size (relative to front)
        [08]: Height of the mixed layer in metres

        [09]: Number of across-front (x) grid points
        [10]: Number of central region (x) grid points
        [11]: Number of along-front (y) grid points
        [12]: Number of vertical (z) grid points

        [13]: Maximum Rossby number of front
        [14]: Minimum Richardson number of front
        [15]: Front buoyancy change

        [16]: Strain rate in 1 / seconds
        [17]: Surface buoyancy loss (positive for cooling) in watts / sq. metre
        [18]: Background stratification in 1 / seconds
        
        [19]: Comment
=#
ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0
using Oceananigans

output_folder = ARGS[1]

# Simulation times
run_time = parse(Float64, ARGS[2])
start_time = parse(Float64, ARGS[3])
save_time = parse(Float64, ARGS[4])

# Coriolis frequency
f = parse(Float64, ARGS[5])

# Grid extent
Lx = parse(Float64, ARGS[6])
Lh = parse(Float64, ARGS[7])
H = parse(Float64, ARGS[8])

# Grid sizes
Nx = parse(Int64, ARGS[9])
Nh = parse(Int64, ARGS[10])
Ny = parse(Int64, ARGS[11])
Nz = parse(Int64, ARGS[12])

# Front
Ro = parse(Float64, ARGS[13])
Ri = parse(Float64, ARGS[14])
Δb = parse(Float64, ARGS[15])

# Background
α = parse(Float64, ARGS[16])
Q = parse(Float64, ARGS[17])
N₀² = parse(Float64, ARGS[18])

comment = ARGS[19]

simulation_parameters = (;
    run_time, start_time, save_time,
    f,
    Lx, Lh, H,
    Nx, Nh, Ny, Nz,
    Ro, Ri, Δb,
    α, Q, N₀²,
    comment
)

include("create_simulation.jl")
