using Oceananigans
using Oceananigans.Operators
using Oceananigans.BoundaryConditions: fill_halo_regions!
# Defines some functions for use in dfm.jl

# Down-front mean
@inline dfm(a) = Field(Average(a; dims=2))

function write_output(u_balance, v_balance, w_balance, b_balance, frame, path)
    jldopen(path, "a") do file
        for (k, v) in pairs(u_balance)
            file["$k/$frame"] = Array(v.data[-2:end, 1:1, -2:end])
        end
        for (k, v) in pairs(v_balance)
            file["$k/$frame"] = Array(v.data[-2:end, 1:1, -2:end])
        end
        for (k, v) in pairs(w_balance)
            file["$k/$frame"] = Array(v.data[-2:end, 1:1, -2:end])
        end
        for (k, v) in pairs(b_balance)
            file["$k/$frame"] = Array(v.data[-2:end, 1:1, -2:end])
        end
    end
    return nothing
end
