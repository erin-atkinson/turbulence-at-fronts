using SpecialFunctions

@inline function get_base_state_vortex(simulation_parameters)
    sp = simulation_parameters
    let H=sp.H, f=sp.f, Ro=sp.Ro, a=sp.a, Ri=sp.Ri
        
        
        @inline g(t) = (t+1)^2 / 2
        @inline g′(t) = (t+1)

        @inline h(t) = exp(-t^2 / 2)
        @inline h′(t) = -t * exp(-t^2 / 2)
        
        Δb = -f^2*a^2*Ro/(H)
        S² = (Δb/(a*f))^2
        N² = S² * Ri - Δb/H

        @inline function b(r, z)
            N²*z + Δb*h(r/a)*g′(z/H)
        end
        @inline function vₜ(r, z)
            H*Δb*h′(r/a)*g(z/H)/(f*a)
        end
        
        @inline function v(r, z)
            -f*r*(0.5 - sqrt(0.25 + vₜ(r, z)/(f*r)))
        end
        w₀(x, y, z) = 0
        
        b₀(x, y, z) = b(sqrt(x^2 + y^2), z)
        u₀(x, y, z) = -v(sqrt(x^2 + y^2), z) * y / sqrt(x^2 + y^2)
        v₀(x, y, z) = v(sqrt(x^2 + y^2), z) * x / sqrt(x^2 + y^2)
        
        xs = range(0, sp.L/2, 500)
        zs = range(-H, 0, 100)
        Δx = xs[2] - xs[1]
        Δz = zs[2] - zs[1]
        Ri_min = minimum([2Δz * (b(r, z+Δz) - b(r, z-Δz))/(v(r, z+Δz) - v(r, z-Δz))^2 for r in xs, z in zs])
        Ro_max = maximum([(v(r+Δx, z) - v(r-Δx, z))/(2Δx * f) for r in xs, z in zs])
        Ro_min = minimum([(v(r+Δx, z) - v(r-Δx, z))/(2Δx * f) for r in xs, z in zs])
        
        @info "Created front state with $(Ri_min) < Ri, $(Ro_min) < Ro < $(Ro_max)"
        
        (; b=b₀, u=u₀, v=v₀, w=w₀)
    end
end

@inline get_base_state(simulation_parameters) = get_base_state_vortex(simulation_parameters)