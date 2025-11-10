# grid_cells.jl

# ---------------------------------------
# Calculating the variable spaced grid
@inline function x_faces(Lx, Lh, Nx, Nh)
    ks = (-Nx ÷ 2):(Nx ÷ 2)
    return linear_to_cubic_x_faces.(ks, Lx, Lh, Nx, Nh)
end

@inline function linear_to_cubic_x_faces(k, Lx, Lh, Nx, Nh)
    a = (Lx/2 - Lh/2) / (Nx÷2 - Nh÷2)^3
    abs(k) <= Nh÷2 && return k * Lh / Nh
    k > Nh÷2 && return Lh/2 + a * (k - Nh÷2)^3
    k < -Nh÷2 && return -Lh/2 + a * (k + Nh÷2)^3
end

@inline function x_width(Nx, Nh, Lh, s)
    return x_faces(Nx, Nh, Lh, s)[end] - x_faces(Nx, Nh, Lh, s)[1]
end
# ---------------------------------------

# ---------------------------------------
# All grid faces
@inline function get_grid_faces(simulation_parameters)
    sp = simulation_parameters
    
    # x spacing varies
    xs = x_faces(sp.Lx, sp.Lh, sp.Nx, sp.Nh)
    
    # Other dimensions are uniform spacing
    ys = (-sp.Ly/2, sp.Ly/2)
    zs = (-sp.Lz, 0)
    
    
    (; xs, ys, zs)
end
# ---------------------------------------
