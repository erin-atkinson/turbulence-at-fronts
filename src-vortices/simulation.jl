ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0

using Oceananigans


include("create_simulation.jl")

# (; run_time, a=25000, f=1e-4, H=100, Nx=200, Ny=200, Nz=50, Ro=0.6, Ri=1, α=1e-5, damping_rate=1e-4, damping_frac=1)

output_folder = ARGS[1]

run_time = Meta.parse(ARGS[2])
a = Meta.parse(ARGS[3])
f = Meta.parse(ARGS[4])
H = Meta.parse(ARGS[5])
Nx = Meta.parse(ARGS[6])
Ny = Meta.parse(ARGS[7])
Nz = Meta.parse(ARGS[8])
Ro = Meta.parse(ARGS[9])
Ri = Meta.parse(ARGS[10])
α = Meta.parse(ARGS[11])
damping_rate = Meta.parse(ARGS[12])
damping_frac = Meta.parse(ARGS[13])

simulation_parameters = (; a, f, H, Nx, Ny, Nz, Ro, Ri, α, damping_rate, damping_frac)

simulation = create_simulation(run_time, output_folder, simulation_parameters)
run!(simulation)
