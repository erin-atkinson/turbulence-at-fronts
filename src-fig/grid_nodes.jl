get_grid(file) = file["serialized/grid"]

function grid_nodes(file; with_halos=false, reshape=false)
    grid = get_grid(file)

    xsᶜ, ysᶜ, zsᶜ = nodes(grid, Center(), Center(), Center(); reshape, with_halos)
    xsᶠ, ysᶠ, zsᶠ = nodes(grid, Face(), Face(), Face(); reshape, with_halos)

    return xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ
end

function center_condition(file)
    sp = simulation_parameters(file)
    return (x, y, z) -> -sp.Lh / 2 < x < sp.Lh / 2
end

function center_indices(file)
    sp = simulation_parameters(file)
    return ((sp.Nx - sp.Nh) ÷ 2):((sp.Nx + sp.Nh) ÷ 2)
end