@inline function x_faces(Nx, Nₕ, Lₕ, Lx, ℓ)
    Δx = Lₕ / Nₕ
    Nl = (Nx - Nₕ) ÷ 2
    Nr = (Nx + Nₕ) ÷ 2
    Lr = Lₕ / 2
    
    # Scale factor
    s =  3 * ((Lx / 2 - Lr) / Δx - (Nx - Nr)) / (Nx - Nr)^3
    
    Δxs = map(1:(Nx ÷ 2)) do n
        n < Nₕ ÷ 2 ? Δx : Δx * (1 + s*(n - (Nₕ ÷ 2))^2) 
    end
    xs = vcat(-reverse(cumsum(Δxs)), [0], cumsum(Δxs))
    
    # Above is for continuous n, make consistent by just scaling result
    xs .* (Lx / (xs[end] - xs[1]))
end


@inline function get_grid_faces(simulation_parameters)
    sp = simulation_parameters
    let Lx=sp.Lx, Ly=sp.Ly, H=sp.H, Nx=sp.Nx, Ny=sp.Ny, Nz=sp.Nz, Nₕ=sp.Nₕ, Lₕ=sp.Lₕ, ℓ=sp.ℓ, Lz=sp.Lz
        
        # stretching scale factor
        
        
        # x spacing varies
        xs = x_faces(Nx, Nₕ, Lₕ, Lx, ℓ)
        
        # Other dimensions are uniform spacing
        ys = (-Ly/2, Ly/2)
        zs = (-Lz, 0)
        
        @info "Created grid faces"
        (; xs, ys, zs)
    end
end
