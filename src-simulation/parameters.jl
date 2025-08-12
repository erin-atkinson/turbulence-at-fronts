# parameters.jl

default_inputs = (; 
    f=1e-4, # Coriolis parameter
    H=100, # Mixed layer depth
    Nx=200, # Zonal grid cells
    Ny=200, # Meridional grid cells
    Nz=50, # Vertical grid cells
    Ro=0.1, # Rossby number
    Ri=2, # Richardson number
    α=1e-5, # Strain rate
    Q=0, # Cooling rate
    c=0.5, # Damping rate / IW frequency
    s=1.05
)

@inline function create_simulation_parameters(input_parameters=(; ))
    ip = (; default_inputs..., input_parameters...)
    let f=ip.f, Ro=ip.Ro, Ri=ip.Ri, Nx=ip.Nx, Ny=ip.Ny, H=ip.H, Q=ip.Q, c=ip.c, s=ip.s
        # Velocity scale is set
        U = 0.5

        # Mixed layer depth change
        δ = 0.1
        
        # Transition width
        λ = 0.03
        
        # Pick a length scale using Rossby number
        ℓ = U / (f * Ro)
        
        Lh = ℓ / 2
        Nh = 2*((3Nx) ÷ 8) # Even
        Lx = x_width(Nx, Nh, Lh, s)
        
        Ly = Lh * Ny / Nh
        Lz = 1.5H
        
        
        # Make a guess of N₀² and a
        N₀² = f^2 * Ro * ℓ^2 / H^2 / A(δ)
        a = Ri * N₀² * δ * H / (sqrt(π) * f^2 * ℓ)
        
        # Physical constants needed to construct surface BC
        αV = 2.0678e-4 # K⁻¹
        cₚ = 4.1819e3 # J kg⁻¹ K⁻¹
        ρ = 1.027e3 # kg m⁻³
        g = 9.81 # m s⁻²
        
        # Buoyancy flux from cooling rate
        B = αV * g * Q / (cₚ * ρ)
        
        # Maximum damping rate
        σ = c * sqrt(N₀²) / (2π)
        
        (; ip..., Lx, Ly, Lz, ℓ, N₀², B, λ, δ, a, σ, Lh, Nh)
    end
end

@inline function create_simulation_parameters(; input_parameters...)
    create_simulation_parameters(input_parameters)
end
