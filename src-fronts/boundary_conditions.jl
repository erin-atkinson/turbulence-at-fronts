using Oceananigans
using Oceananigans.Operators

# Boundary conditions should be closed?

@inline function create_boundary_conditions(base_state, simulation_parameters)
    sp = simulation_parameters
    Lx = sp.Lx
    Ly = sp.Ly
    # z is Bounded but cooled on the top
    # y is more complicated, I'm going to try a couple of things
    #=
    # I think that we should have that y gradients are 0
    u_bcs=FieldBoundaryConditions(;
        north=GradientBoundaryCondition(0),
        south=GradientBoundaryCondition(0)
        )
    # I don't think that this works, perhaps needs to be open
    
    
    =#
    b_func = get_base_state_function(simulation_parameters).b
    
    w_bcs=FieldBoundaryConditions(;
        east=ValueBoundaryCondition(0),
        west=ValueBoundaryCondition(0),
        )
    
    # Prevent reflections
    b_bcs = if sp.α == 0 
        FieldBoundaryConditions(;
            bottom=GradientBoundaryCondition(simulation_parameters.N₀²),
            top=FluxBoundaryCondition(simulation_parameters.B)
        )
    else
        FieldBoundaryConditions(;
            bottom=GradientBoundaryCondition(simulation_parameters.N₀²),
            top=FluxBoundaryCondition(simulation_parameters.B),
            west=ValueBoundaryCondition((y, z, t)->b_func(-Lx/2, z)),
            east=ValueBoundaryCondition((y, z, t)->b_func(Lx/2, z))
        )
    end
    
    v_bcs=FieldBoundaryConditions(;
        bottom=ValueBoundaryCondition(0)
        )
    @info "Created boundary conditions"
    (; v=v_bcs, w=w_bcs, b=b_bcs)
end