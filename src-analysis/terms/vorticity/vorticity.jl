@inline function ηx_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    return ∂yᶜᶠᶠ(i, j, k, grid, fields.w) - ∂zᶜᶠᶠ(i, j, k, grid, fields.v)
end

@inline function ηy_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    return ∂zᶠᶜᶠ(i, j, k, grid, fields.u) - ∂xᶠᶜᶠ(i, j, k, grid, fields.w)
end

@inline function ηz_func(i, j, k, grid, clock, fields, dependency_fields, sp)
    return ∂xᶠᶠᶜ(i, j, k, grid, fields.v) - ∂yᶠᶠᶜ(i, j, k, grid, fields.u) + sp.f
end

η_dependencies = ()
