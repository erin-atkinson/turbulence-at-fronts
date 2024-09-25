include("base_state.jl")

default_inputs = (; f=1e-4, H=100, Nx=200, Nₕ=100, Ny=200, Nz=50, Ro=0.1, Ri=2, α=1e-5, Q=0, c=0.5, damping_frac=0.5, turbulence_spinup=1e4, Ek=0, )

@inline function create_simulation_parameters(input_parameters=(; ))
    ip = (; default_inputs..., input_parameters...)
    let f=ip.f, Ro=ip.Ro, Ri=ip.Ri, Nx=ip.Nx, Ny=ip.Ny, Nz=ip.Nz, H=ip.H, Q=ip.Q, Ek=ip.Ek, Nₕ=ip.Nₕ, α=ip.α
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
        αV = 2.0678e-4 # K⁻¹
        cₚ = 4.1819e3 # J kg⁻¹ K⁻¹
        
        ρ = 1.027e3 # kg m⁻³
        g = 9.81 # m s⁻²
        
        # Buoyancy flux
        B = αV * g * Q / (cₚ * ρ)
        
        Pr = 1
        
        # width of high-resolution area
        Lₕ = 2.5ℓ
        # total width
        Lx = 16ℓ
        # Down-front length should be such that the the grid cells are isotropic in horizontal
        # within the high-resolution area
        Ly = Lₕ * Ny / Nₕ
        Lz = 1.2H
        
        #νₕ = 0.1 * max(α, f / 10) * Lₕ^2 / (4Nₕ^2)
        #νᵥ = νₕ * (H * Nₕ / (Lₕ * Nz))^2
        
        (; ip..., Lx, Lₕ, Ly, Lz, ℓ, N₀², B, λ, δ, ε, a, ν, Pr)
    end
end

@inline function create_simulation_parameters(; input_parameters...)
    create_simulation_parameters(input_parameters)
end