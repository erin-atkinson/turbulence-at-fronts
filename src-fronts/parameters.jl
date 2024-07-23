include("base_state.jl")

default_inputs = (; f=1e-4, H=100, Nx=200, Ny=200, Nz=50, Ro=0.1, Ri=2, α=1e-5, Q=0, c=0.5, damping_frac=0.5, turbulence_spinup=1e4, Ek=0)

@inline function create_simulation_parameters(input_parameters=(; ))
    ip = (; default_inputs..., input_parameters...)
    let f=ip.f, Ro=ip.Ro, Ri=ip.Ri, Nx=ip.Nx, Ny=ip.Ny, H=ip.H, Q=ip.Q, Ek=ip.Ek
        # Velocity scale
        U = 0.1
        # Mixed layer depth change
        δ = 0.1
        # Transition width
        λ = 0.01
        # Mixed layer stratification relative to deep water
        ε = 0.0001
        ℓ = U / (f * Ro)
        # Compute the appropriate parameters here
        
        # Choose a with an iterative process
        
        # First make a guess of N₀²
        
        N₀² = f^2 * Ro * ℓ^2 / H^2 / A(δ)
        
        # This gives a guess of a
        a = Ri * N₀² * δ * H / (sqrt(π) * f^2 * ℓ)
        
        # Refine the guess of a?
        
        ν = Ek * f * H
        # Thermal properties of water
        α = 2.0678e-4 # K⁻¹
        cₚ = 4.1819e3 # J kg⁻¹ K⁻¹
        
        ρ = 1.027e3 # kg m⁻³
        g = 9.81 # m s⁻²
        
        # Buoyancy flux
        B = α * g * Q / (cₚ * ρ)
        
        Pr = 1
        
        Lx = 5ℓ
        # Down-front length should be such that the the grid cells are isotropic in horizontal
        Ly = Lx * Ny / Nx
        (; ip..., Lx, Ly, ℓ, N₀², B, λ, δ, ε, a, ν, Pr)
    end
end

@inline function create_simulation_parameters(; input_parameters...)
    create_simulation_parameters(input_parameters)
end