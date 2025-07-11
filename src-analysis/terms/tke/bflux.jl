@inline function BFLUX3D_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    
    w = fields.w
    b = fields.b

    b_dfm = dependency_fields.b_dfm

    W = fields.W

    total_w = SumOfArrays{2}(w, W)
    
    wb = advective_tracer_flux_density_z(i, j, k, grid, weno, total_w, b)
    wb_dfm = advective_tracer_flux_density_z(i, j, k, grid, weno, total_w, b_dfm)
    
    return wb - wb_dfm
end

BFLUX_dependencies = (
    :b_dfm,
)
