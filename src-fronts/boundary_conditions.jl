# Boundary condition on buoyancy should agree with the forcing

@inline b_bcs_west_func(j, k, grid, clock, model_fields, p) = @inbounds p.b_west[1, 1, k]
@inline b_bcs_east_func(j, k, grid, clock, model_fields, p) = @inbounds p.b_east[1, 1, k]
    
w_bcs=FieldBoundaryConditions(;
    east=ValueBoundaryCondition(0),
    west=ValueBoundaryCondition(0),
)

# Prevent reflections
#=
b_bcs = FieldBoundaryConditions(;
    bottom=GradientBoundaryCondition(sp.N₀²),
    top=FluxBoundaryCondition(sp.B),
    west=ValueBoundaryCondition(b_bcs_west_func; discrete_form=true, parameters=(; b_west)),
    east=ValueBoundaryCondition(b_bcs_east_func; discrete_form=true, parameters=(; b_east))
)


b_bcs = FieldBoundaryConditions(;
    bottom=GradientBoundaryCondition(sp.N₀²),
    top=FluxBoundaryCondition(sp.B),
    west=ValueBoundaryCondition((y, z, t, sp)->base_b_func(-sp.Lx/2, z, sp); parameters=sp),
    east=ValueBoundaryCondition((y, z, t, sp)->base_b_func(sp.Lx/2, z, sp); parameters=sp)
)
=#

b_flux_func(x, y, t, p) = p.B * (1 - exp(-p.f*(t - p.start_time) / 20))

b_bcs = FieldBoundaryConditions(;
    bottom=GradientBoundaryCondition(sp.N₀²),
    top=FluxBoundaryCondition(b_flux_func; parameters=(; sp.B, sp.f, sp.start_time))
)

v_bcs=FieldBoundaryConditions(;
    bottom=ValueBoundaryCondition(0)
)

ξ_bcs = FieldBoundaryConditions(;
    bottom=ValueBoundaryCondition(1),
)

@info "Created boundary conditions"
boundary_conditions = (; v=v_bcs, b=b_bcs, w=w_bcs)#, ξ=ξ_bcs)