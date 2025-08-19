# grid_cells.jl

# ---------------------------------------
# Calculating the variable spaced grid
@inline function x_faces(Nx, Nh, Lh, s)
    ks = (-Nx ÷ 2):(Nx ÷ 2)
    return uniform_to_exponential_x_faces.(ks, Nh, Lh, s)
end

@inline function uniform_to_exponential_x_faces(k, Nh, Lh, s)
    abs(k) <= Nh÷2 && return k * Lh / Nh
    k > Nh÷2 && return Lh/2 + Lh * (s^(k - Nh÷2) - 1) / (log(s)*Nh)
    k < -Nh÷2 && return -Lh/2 - Lh * (s^(-k - Nh÷2) - 1) / (log(s)*Nh)
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
    xs = x_faces(sp.Nx, sp.Nh, sp.Lh, sp.s)
    
    # Other dimensions are uniform spacing
    ys = (-sp.Ly/2, sp.Ly/2)
    zs = (-sp.Lz, 0)
    
    
    (; xs, ys, zs)
end
# ---------------------------------------

grid_faces = get_grid_faces(sp)
@info "Created grid faces"
