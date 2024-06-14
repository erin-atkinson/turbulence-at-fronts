#=
analysis_functions.jl
    Defines some functions that help with reading data from simulation output

=#

using JLD2

function get_simulation_parameters(foldername)
    # Get the saved simulation parameters
    jldopen("$foldername/parameters.jld2") do file
        NamedTuple(map(k->Pair(Symbol(k), file[k]), keys(file)))
    end
end

function make_GPU_grid(grid)
    # Create a GPU grid from a CPU grid
    # (is there really no way to do this directly?)
    return grid
end

# Function to update input fields
function update_fields!(input_fields, frame, path)
    jldopen(path) do file
        input_fields.u .= file["timeseries/u/$frame"]
        input_fields.v .= file["timeseries/v/$frame"]
        input_fields.w .= file["timeseries/w/$frame"]
        input_fields.b .= file["timeseries/b/$frame"]
        input_fields.φ .= file["timeseries/φ/$frame"]
        map(fill_halo_regions!, input_fields)
    end
    return nothing
end