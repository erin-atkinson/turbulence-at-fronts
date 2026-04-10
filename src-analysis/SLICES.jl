using Oceananigans: location
include("terms/slices.jl")

# x, y, z slices of the model fields
dependency_fields = NamedTuple()

# x = Lh/2, x = -Lh/2
for ξ in (:u, :v, :w, :b)
    left_ξ = Symbol(:left_, ξ)
    right_ξ = Symbol(:right_, ξ)
    @info "Slices $left_ξ and $right_ξ"

    @eval begin
        field = input_fields.$ξ
        (Lx, Ly, Lz) = location(field)

        $left_ξ = Field(KernelFunctionOperation{Nothing, Ly, Lz}(x_slice_func, grid, field, -sp.Lh/2, (Lx(), Ly(), Lz())))
        $right_ξ = Field(KernelFunctionOperation{Nothing, Ly, Lz}(x_slice_func, grid, field, sp.Lh/2, (Lx(), Ly(), Lz())))

        dependency_fields = merge(dependency_fields, (; $left_ξ, $right_ξ))
    end
end

# y = Ly/2
for ξ in (:u, :v, :w, :b)
    front_ξ = Symbol(:front_, ξ)
    @info "Slices $front_ξ"

    
    @eval begin
        field = input_fields.$ξ
        (Lx, Ly, Lz) = location(field)

        $front_ξ = Field(KernelFunctionOperation{Lx, Nothing, Lz}(y_slice_func, grid, field, sp.Ly/2, (Lx(), Ly(), Lz())))

        dependency_fields = merge(dependency_fields, (; $front_ξ))
    end
end

# z = -0.05H, -0.25H, -0.5H, -0.75H, -0.95H
for ξ in (:u, :v, :w, :b)
    z005_ξ = Symbol(:z005_, ξ)
    z025_ξ = Symbol(:z025_, ξ)
    z050_ξ = Symbol(:z050_, ξ)
    z075_ξ = Symbol(:z075_, ξ)
    z095_ξ = Symbol(:z095_, ξ)
    @info "Slices $z005_ξ, $z025_ξ, $z050_ξ, $z075_ξ and $z095_ξ"

    @eval begin
        field = input_fields.$ξ
        (Lx, Ly, Lz) = location(field)

        $z005_ξ = Field(KernelFunctionOperation{Lx, Ly, Nothing}(z_slice_func, grid, field, -0.05sp.H, (Lx(), Ly(), Lz())))
        $z025_ξ = Field(KernelFunctionOperation{Lx, Ly, Nothing}(z_slice_func, grid, field, -0.25sp.H, (Lx(), Ly(), Lz())))
        $z050_ξ = Field(KernelFunctionOperation{Lx, Ly, Nothing}(z_slice_func, grid, field, -0.50sp.H, (Lx(), Ly(), Lz())))
        $z075_ξ = Field(KernelFunctionOperation{Lx, Ly, Nothing}(z_slice_func, grid, field, -0.75sp.H, (Lx(), Ly(), Lz())))
        $z095_ξ = Field(KernelFunctionOperation{Lx, Ly, Nothing}(z_slice_func, grid, field, -0.95sp.H, (Lx(), Ly(), Lz())))

        dependency_fields = merge(dependency_fields, (; $z005_ξ, $z025_ξ, $z050_ξ, $z075_ξ, $z095_ξ, ))
    end
end

skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)

output_fields = dependency_fields
