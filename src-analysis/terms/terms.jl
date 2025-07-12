using Oceananigans.Operators

@inline down_front_mean(a) = Field(Average(a; dims=2))
@inline dfm(a) = down_front_mean(a)

# Calculate a down front-mean of a field
# Grid spacing is uniform in down-front direction
@inline function discrete_down_front_mean(i, j, k, grid, f, args...)
    j_indices = (grid.Hy+1):(grid.Hy+grid.Ny)
    total = mapreduce(+, j_indices) do _j
        f(i, _j, k, grid, args...)
    end
    total / grid.Ny
end

@inline function discrete_down_front_mean(i, j, k, grid, field)
    j_indices = (grid.Hy+1):(grid.Hy+grid.Ny)
    total = mapreduce(+, j_indices) do _j
        @inbounds field[i, _j, k]
    end
    total / grid.Ny
end

@inline ddfm(args...) = discrete_down_front_mean(args...)

@inline fg(i, j, k, grid, f, g) = @inbounds f[i, j, k] * g[i, j, k]
@inline fGg(i, j, k, grid, f, G, args...) = @inbounds f[i, j, k] * G(i, j, k, grid, args...)
@inline FfGg(i, j, k, grid, F, f, G, args...) = @inbounds F(i, j, k, grid, f) * G(i, j, k, grid, args...)

@inline a_avg(i, j, k, grid, a, a_next) = @inbounds (a[i, j, k] + a_next[i, j, k]) / 2

@inline f_avg_Gg(i, j, k, grid, f, f_next, G, args...) = a_avg(i, j, k, grid, f, f_next) * G(i, j, k, grid, args...)
@inline f′_avg_Gg(i, j, k, grid, f, f_next, f_dfm, f_next_dfm, G, args...) = (a_avg(i, j, k, grid, f, f_next) - a_avg(i, j, k, grid, f_dfm, f_next_dfm)) * G(i, j, k, grid, args...)

@inline df′dt(i, j, k, grid, f, f_next, f_dfm, f_next_dfm, Δt) = (f′(i, j, k, grid, f_next, f_next_dfm) - f′(i, j, k, grid, f, f_dfm)) / Δt

@inline f′(i, j, k, grid, f, f_dfm) = @inbounds f[i, j, k] - f_dfm[i, j, k]
@inline f′g′(i, j, k, grid, f, f_dfm, g, g_dfm) = f′(i, j, k, grid, f, f_dfm) * f′(i, j, k, grid, g, g_dfm)
@inline f′Gg′(i, j, k, grid, f, f_dfm, G, g, g_dfm) = f′(i, j, k, grid, f, f_dfm) * G(i, j, k, grid, f′, g, g_dfm)
