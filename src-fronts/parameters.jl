default_inputs = (; f=1e-4, H=100, Nx=200, Ny=200, Nz=50, Ro=0.1, Ri=2, α=1e-5, Q=0, c=0.5, damping_frac=0.5, turbulence_spinup=1e4)

@inline function create_simulation_parameters(input_parameters=(; ))
    ip = (; default_inputs..., input_parameters...)
    let f=ip.f, Ro=ip.Ro, Ri=ip.Ri, Nx=ip.Nx, Ny=ip.Ny, H=ip.H, Q=ip.Q
        # Velocity scale
        U = 0.1
        # Mixed layer depth change
        δ = 0.1
        # Transition width
        λ = 0.001
        
        ℓ = U / (f * Ro)
        
        ΔN² = -f^2 * Ro * ℓ^2 * sqrt(π * exp(1)) / (2 * δ * H^2)
        N² = 4Ri * ΔN²^2 * δ^2 * H^2 / (π * f^2 * ℓ^2)
        N₀² = N² - ΔN²
        
        # Thermal properties of water
        α = 2.0678e-4 # K⁻¹
        cₚ = 4.1819e3 # J kg⁻¹ K⁻¹
        
        ρ = 1.027e3 # kg m⁻³
        g = 9.81 # m s⁻²
        
        # Buoyancy flux
        B = α * g * Q / (cₚ * ρ)
        
        Lx = 6ℓ
        # Down-front length should be such that the the grid cells are isotropic in horizontal
        Ly = Lx * Ny / Nx
        (; ip..., Lx, Ly, ℓ, ΔN², N², N₀², B, λ, δ)
    end
end

@inline function create_simulation_parameters(; input_parameters...)
    create_simulation_parameters(input_parameters)
end