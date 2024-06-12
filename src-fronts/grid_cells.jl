@inline function get_grid_faces(simulation_parameters)
    sp = simulation_parameters
    let Lx=sp.Lx, Ly=sp.Ly, H=sp.H, Nx=sp.Nx, Ny=sp.Ny, Nz=sp.Nz
        xs = (-Lx/2, Lx/2)
        ys = (-Ly/2, Ly/2)
        zs = (-H, 0)
        @info "Created grid faces"
        (; xs, ys, zs)
    end
end
