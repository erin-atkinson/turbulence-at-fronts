using Oceananigans
using Oceananigans.Operators

# Boundary conditions should be closed?

@inline function create_boundary_conditions(base_state, simulation_parameters)
    Lx = simulation_parameters.L
    Ly = simulation_parameters.L * (simulation_parameters.Ny / simulation_parameters.Nx)
    
    u_bcs=FieldBoundaryConditions(;
        north=ValueBoundaryCondition((x, z, t, b)->u(x, Ly/2, z); parameters=base_state.u),
        south=ValueBoundaryCondition((x, z, t, b)->u(x, -Ly/2, z); parameters=base_state.u),
        east=ValueBoundaryCondition((y, z, t, b)->u(Lx/2, y, z); parameters=base_state.u),
        west=ValueBoundaryCondition((y, z, t, b)->u(-Lx/2, y, z); parameters=base_state.u)
        )
    v_bcs=FieldBoundaryConditions(;
        north=ValueBoundaryCondition((x, z, t, b)->v(x, Ly/2, z); parameters=base_state.v),
        south=ValueBoundaryCondition((x, z, t, b)->v(x, -Ly/2, z); parameters=base_state.v),
        east=ValueBoundaryCondition((y, z, t, b)->v(Lx/2, y, z); parameters=base_state.v),
        west=ValueBoundaryCondition((y, z, t, b)->v(-Lx/2, y, z); parameters=base_state.v)
        )
    w_bcs=FieldBoundaryConditions(;
        north=ValueBoundaryCondition((x, z, t, b)->w(x, Ly/2, z); parameters=base_state.w),
        south=ValueBoundaryCondition((x, z, t, b)->w(x, -Ly/2, z); parameters=base_state.w),
        east=ValueBoundaryCondition((y, z, t, b)->w(Lx/2, y, z); parameters=base_state.w),
        west=ValueBoundaryCondition((y, z, t, b)->w(-Lx/2, y, z); parameters=base_state.w)
        )
    b_bcs=FieldBoundaryConditions(;
        north=ValueBoundaryCondition((x, z, t, b)->b(x, Ly/2, z); parameters=base_state.b),
        south=ValueBoundaryCondition((x, z, t, b)->b(x, -Ly/2, z); parameters=base_state.b),
        east=ValueBoundaryCondition((y, z, t, b)->b(Lx/2, y, z); parameters=base_state.b),
        west=ValueBoundaryCondition((y, z, t, b)->b(-Lx/2, y, z); parameters=base_state.b)
        )
    @info "Created boundary conditions"
    (; b=b_bcs)
end