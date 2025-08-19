# boundary_conditions.jl

# ---------------------------------------
# Cooling turns on slowly
@inline function b_flux_func(x, y, t, sp) 
    turnon = 1 - exp(-sp.f*(t - sp.start_time) / 20)
    return p.B * turnon
end
# ---------------------------------------

# ---------------------------------------
# Boundary conditions have no flow outside of the domain
b_bcs = FieldBoundaryConditions(;
    bottom=GradientBoundaryCondition(sp.N₀²),
    top=FluxBoundaryCondition(b_flux_func; parameters=(; sp.B, sp.f, sp.start_time))
)
v_bcs=FieldBoundaryConditions(;
    bottom=ValueBoundaryCondition(0),
    east=ValueBoundaryCondition(0),
    west=ValueBoundaryCondition(0)
)
u_bcs=FieldBoundaryConditions(;
    bottom=ValueBoundaryCondition(0),
)
w_bcs=FieldBoundaryConditions(;
    east=ValueBoundaryCondition(0),
    west=ValueBoundaryCondition(0)
)
# ---------------------------------------

@info "Created boundary conditions"
boundary_conditions = (; v=v_bcs, b=b_bcs, w=w_bcs)