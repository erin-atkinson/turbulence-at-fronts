#= 
dfm.jl
    Computes the down-front (y) mean and first-order correlations of the simulation output

    Arguments
        [1]: Folder containing simulation output, will also output to this folder
    Pass options in last argument
        i: use initialisation run
=#

include("analysis_functions.jl")
include("dfm_functions.jl")

foldername = ARGS[1]

opts = length(ARGS) > 1 ? ARGS[2] : ""
filename = 'i' ∈ opts ? "initialisation" : "output"

# Get grid and iterations
grid, frames, ts = jldopen("$input_folder/$filename.jld2") do file
    frames = parse.(Int, keys(file["timeseries/t"]))
    file["serialized/grid"], frames, [file["timeseries/t/$frame"] for frame in frames]
end
# Get simulation parameters
sp = get_simulation_parameters(foldername)

# input: u, v, w, b, ϕ
# Create a field for each by loading the first iteration in file

u = FieldTimeSeries("$input_folder/$filename.jld2", "u"; iterations=frames[51])[1]
v = FieldTimeSeries("$input_folder/$filename.jld2", "v"; iterations=frames[51])[1]
w = FieldTimeSeries("$input_folder/$filename.jld2", "w"; iterations=frames[51])[1]
b = FieldTimeSeries("$input_folder/$filename.jld2", "b"; iterations=frames[51])[1]
φ = FieldTimeSeries("$input_folder/$filename.jld2", "φ"; iterations=frames[51])[1]

input_fields = (; u, v, w, b, φ)

# ᶜⁿᶠ
# Mean fields
uᶠⁿᶜ = Field(Average(u; dims=2))
vᶜⁿᶜ = Field(Average(v; dims=2))
wᶜⁿᶠ = Field(Average(w; dims=2))
bᶜⁿᶜ = Field(Average(b; dims=2))
φᶜⁿᶜ = Field(Average(φ; dims=2))

mean_fields = (; uᶠⁿᶜ, vᶜⁿᶜ, wᶜⁿᶠ, bᶜⁿᶜ, φᶜⁿᶜ)

# Correlation fields
#  Self
u′u′ᶠⁿᶜ = dfm(KernelFunctionOperation{Face, Center, Center}(variance_kernel_func, grid, u, u_dfm))
v′v′ᶜⁿᶜ = dfm(KernelFunctionOperation{Center, Face, Center}(variance_kernel_func, grid, u, u_dfm))
w′w′ᶜⁿᶠ = dfm(KernelFunctionOperation{Center, Center, Face}(variance_kernel_func, grid, u, u_dfm))

#  Cross terms
u′v′ᶜⁿᶜ = dfm(KernelFunctionOperation{Center, Center, Center}(u′v′ᶜᶜᶜ_kernel_func, grid, u, u_dfm, v, v_dfm))
v′w′ᶜⁿᶜ = dfm(KernelFunctionOperation{Center, Center, Center}(v′w′ᶜᶜᶜ_kernel_func, grid, v, v_dfm, w, w_dfm))
w′u′ᶜⁿᶜ = dfm(KernelFunctionOperation{Center, Center, Center}(w′u′ᶜᶜᶜ_kernel_func, grid, w, w_dfm, u, u_dfm))

#  Buoyancy
u′b′ᶠⁿᶜ = dfm(KernelFunctionOperation{Face, Center, Center}(u′b′ᶠᶜᶜ_kernel_func, grid, u, u_dfm, b, b_dfm))
v′b′ᶜⁿᶜ = dfm(KernelFunctionOperation{Center, Face, Center}(v′b′ᶜᶠᶜ_kernel_func, grid, v, v_dfm, b, b_dfm))
w′b′ᶜⁿᶠ = dfm(KernelFunctionOperation{Center, Center, Face}(w′b′ᶜᶜᶠ_kernel_func, grid, w, w_dfm, b, b_dfm))

correlation_fields = (; u′u′ᶠⁿᶜ, v′v′ᶜⁿᶜ, w′w′ᶜⁿᶠ)
# Now we have defined all our operations, we loop over the frames

output_path = "$input_folder/dfm.jld2"

for frame in frames
    # Get the data
    update_fields!(input_fields, frame, "$input_folder/$filename.jld2")
    # Compute 
    map(compute!, mean_fields)
    map(compute!, correlation_fields)
    # Save to file
    write_output(mean_fields, correlation_fields, frame, output_path)
end
write_grid_times(grid, frames, ts, output_path)


