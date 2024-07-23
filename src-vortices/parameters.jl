default_inputs = (; f=1e-4, H=100, Nx=200, Ny=200, Nz=50, Ro=0.6, Ri=1, α=1e-5, damping_rate=1e-4, damping_frac=1, a=25000)

@inline function create_simulation_parameters(input_parameters=(; ))
    ip = (; default_inputs..., input_parameters...)
    let a=ip.a
        L = 7a
        (; ip..., L)
    end
end

@inline function create_simulation_parameters(; input_parameters...)
    create_simulation_parameters(input_parameters)
end