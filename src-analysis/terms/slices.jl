using Oceananigans.Grids: node
using Oceananigans.Fields: interpolate

@inline function x_slice_func(i, j, k, grid, field, x, loc)
    (_, y, z) = node(i, j, k, grid, loc...)
    at_node = (x, y, z)

    return interpolate(at_node, field, loc, grid)
end

@inline function y_slice_func(i, j, k, grid, field, y, loc)
    (x, _, z) = node(i, j, k, grid, loc...)
    at_node = (x, y, z)
    
    return interpolate(at_node, field, loc, grid)
end

@inline function z_slice_func(i, j, k, grid, field, z, loc)
    (x, y, _) = node(i, j, k, grid, loc...)
    at_node = (x, y, z)
    
    return interpolate(at_node, field, loc, grid)
end
