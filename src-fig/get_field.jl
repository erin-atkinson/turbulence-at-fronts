nov = no_offset_view

# Get field directly from the JLD2 file
@inline function get_field(file, field, iteration; halo=false)
    
    jldpath = "timeseries/$field/$iteration"
    data = nov(file[jldpath])
    
    # Remove singleton dimensions and halo points
    Hs = halo ? zeros(Int, 3) : halos(file)
    indices = map(Hs, size(data)) do H, s
        s==1 ? 1 : (1+H):(s-H)
    end
    
    data[indices...]
end

@inline function get_field(files::AbstractVector, fields::AbstractVector, iteration; halo=false)
    map(files, fields) do file, field
        get_field(file, field, iteration; halo)
    end
end

# Get from a string
@inline function get_field(filename::String, args...; kwargs...)
    jldopen(file->get_field(file, args...; kwargs...), filename)
end

# Apply a function to the field after getting
@inline function get_field(f::Function, args...; kwargs...)
    f(get_field(args...; kwargs...))
end