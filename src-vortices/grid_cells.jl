@inline function get_grid_faces(simulation_parameters)
    sp = simulation_parameters
    let L=sp.L, H=sp.H, Nx=sp.Nx, Ny=sp.Ny, Nz=sp.Nz
        zs = (-H, 0)
        xs = (-L/2, L/2)
        ys = xs .* (Ny / Nx)
        @info "Created grid faces"
        (; xs, ys, zs)
    end
end
