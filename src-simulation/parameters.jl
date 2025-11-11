# parameters.jl

default_inputs = (;
    run_time = 1e6, start_time = -5e6, save_time = 1e3,
    f = 1e-4,
    Lx = 2, Lh = 0.5, H = 100,
    Nx = 1024, Nh = 800, Ny = 128, Nz = 128,
    Ro = 0.1, Ri = 2, Δb = 5.5e-4, 
    α = 1e-5, Q = 100, N₀² = 1e-4,
    comment = ""
)

@inline function create_simulation_parameters(input_parameters=(; ))
    ip = (; default_inputs..., input_parameters...)

    sp = create_front_parameters(ip)

    # Domain size
    Lh = sp.Lh * sp.ℓ # About 2 for Ro = 0.4 works
    Lx = sp.Lx * sp.ℓ
    Ly = sp.Lh * sp.Ny / sp.Nh
    Lz = 1.5sp.H

    # Physical constants needed to construct surface BC
    αV = 2.0678e-4 # K⁻¹
    cₚ = 4.1819e3 # J kg⁻¹ K⁻¹
    ρ = 1.027e3 # kg m⁻³
    g = 9.81 # m s⁻²

    # Surface buoyancy flux
    B = αV * g * sp.Q / (cₚ * ρ)

    # Sponge layer damping rate
    σ = 0.5 * sqrt(sp.N₀²) / (2π)
    
    return merge(sp, (; Lh, Lx, Ly, Lz, B, σ))
end

@inline function create_simulation_parameters(; input_parameters...)
    create_simulation_parameters(input_parameters)
end
