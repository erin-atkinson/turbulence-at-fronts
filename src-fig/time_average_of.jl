# Time average of a function
@inline function time_average_of(func, file, field, iterations)
    mapreduce(.+, iterations) do iteration
        get_field(func, file, field, iteration)
    end ./ length(iterations)
end

# From a string
@inline function time_average_of(func, filename::String, args...)
    jldopen(filename) do file
        time_average_of(func, file, args...)
    end
end