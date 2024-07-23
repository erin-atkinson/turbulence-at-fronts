using Oceananigans

@inline function create_closure(xs, ys, zs, base_state, simulation_parameters)
    sp = simulation_parameters
    
    sp.Ek == 0 && return AnisotropicMinimumDissipation()
    
    closure = (VerticalScalarDiffusivity(; ν=sp.ν, κ=sp.ν / sp.Pr), AnisotropicMinimumDissipation())
    @info "Created closure"
    return closure
end