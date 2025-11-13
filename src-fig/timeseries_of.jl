# Take a vector of nd arrays and turn it into a single n+1d array
function expand(v::V) where {T, n, A<:Array{T, n}, V<:Vector{A}}
    l = length(v)
    s = size(v[1])
    new_size = (l, s...)
    a = zeros(T, new_size)
    
    inds = [Colon() for x in s] # There must be some syntax I'm missing
    for i in 1:l
        a[i, inds...] .= v[i]
    end
    a
end

expand(v::V) where {T, V<:Vector{T}} = v

# Reduce Nd field using func (Nd->Md), then find it for all times
@inline function timeseries_of(args...)
    timeseries_of(identity, args...)
end

@inline function timeseries_of(func::Function, file, field, iterations)
    data = map(iteration->get_field(func, file, field, iteration), iterations)
    expand(data)
end

# From a string
@inline function timeseries_of(func::Function, filename::String, args...)
    jldopen(filename) do file
        timeseries_of(func, file, args...)
    end
end
