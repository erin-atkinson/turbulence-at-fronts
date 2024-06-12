using SpecialFunctions

@inline g(s) = s < 0 ? 0 : sqrt(1 + s^2) - 1
@inline γ(s, δ) = -1 + δ * (erf(s) + 1) / 2
@inline γ′(s, δ) = 2δ * exp(-s^2 / 2) / sqrt(π)

#= 
The base/reference state is a front above a mixed layer, with the varying mixed
layer height providing the 
=#
#@inline function get_base_state_filament(simulation_parameters) end

@inline function get_base_state_front(simulation_parameters)
    sp = simulation_parameters
    # Create the front
    
    @inline h(x) = sp.H * γ(x/sp.ℓ, sp.δ)
    @inline h′(x) = sp.H * γ′(x/sp.ℓ, sp.δ) / sp.ℓ
    
    @inline function b(x, z)
        sp.N₀² * z + sp.ΔN² * sp.λ * sp.H * g((z - h(x))/(sp.λ * sp.H))
    end
    @inline function v(x, z)
        -sp.ΔN² * sp.λ * sp.H * g((z - h(x))/(sp.λ * sp.H)) * h′(x) / sp.f
    end
    
    @inline w₀(x, y, z) = 0
    @inline b₀(x, y, z) = b(x, z)
    @inline u₀(x, y, z) = 0
    @inline v₀(x, y, z) = v(x, z)
    
    (; b=b₀, u=u₀, v=v₀, w=w₀)
end

@inline get_base_state(simulation_parameters) = get_base_state_front(simulation_parameters)