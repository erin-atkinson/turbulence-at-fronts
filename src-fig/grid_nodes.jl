get_grid(file) = file["serialized/grid"]

function serialized_grid_nodes(file; with_halos=false, reshape=false)
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

nov = no_offset_view
@inline halos(file) = file["grid/Hx"], file["grid/Hz"], file["grid/Hz"]

function grid_nodes(file; with_halos=false, reshape=false)
    xsᶜ = nov(file["grid/xᶜᵃᵃ"])
    xsᶠ = nov(file["grid/xᶠᵃᵃ"])

    ysᶜ = nov(file["grid/yᵃᶜᵃ"])
    ysᶠ = nov(file["grid/yᵃᶠᵃ"])

    zsᶠ = nov(file["grid/z/cᵃᵃᶠ"])
    zsᶜ = nov(file["grid/z/cᵃᵃᶜ"])

    if !with_halos
        Hx, Hy, Hz = halos(file)

        xsᶜ = xsᶜ[(Hx+1):(end-Hx)]
        xsᶠ = xsᶠ[(Hx+1):(end-Hx)]
        ysᶜ = ysᶜ[(Hy+1):(end-Hy)]
        ysᶠ = ysᶠ[(Hy+1):(end-Hy)]
        zsᶜ = zsᶜ[(Hz+1):(end-Hz)]
        zsᶠ = zsᶠ[(Hz+1):(end-Hz)]
    end

    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ
end