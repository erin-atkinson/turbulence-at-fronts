#= 
dfm.jl
    Computes the down-front (y) mean and first-order correlations of the simulation output

    Arguments
        [1]: Folder containing simulation output, will also output to this folder
    Pass options in last argument
        i: use initialisation run
=#
ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0
include("analysis_functions.jl")
include("dfm_functions.jl")

foldername = ARGS[1]
opts = length(ARGS) > 1 ? ARGS[2] : ""
filename = 'i' ∈ opts ? "initialisation" : "output"

# Get grid and iterations
frames, ts = jldopen("$foldername/$filename.jld2") do file
    frames = parse.(Int, keys(file["timeseries/t"]))
    frames, [file["timeseries/t/$frame"] for frame in frames]
end
grid = jldopen("$foldername/$filename.jld2") do file
    file["serialized/grid"]
end

@info grid
# Get simulation parameters
sp = get_simulation_parameters(foldername)
@info sp
# input: u, v, w, b, ϕ
# Create a field for each by loading the first iteration in file

u = FieldTimeSeries("$foldername/$filename.jld2", "u"; iterations=frames[1])[1]
v = FieldTimeSeries("$foldername/$filename.jld2", "v"; iterations=frames[1])[1]
w = FieldTimeSeries("$foldername/$filename.jld2", "w"; iterations=frames[1])[1]
b = FieldTimeSeries("$foldername/$filename.jld2", "b"; iterations=frames[1])[1]
φ = FieldTimeSeries("$foldername/$filename.jld2", "φ"; iterations=frames[1])[1]

input_fields = (; u, v, w, b, φ)
@info "Created input fields"
# ᶜⁿᶠ
# Mean fields
uᶠⁿᶜ = dfm(u)
vᶜⁿᶜ = dfm(v)
wᶜⁿᶠ = dfm(w)
bᶜⁿᶜ = dfm(b)
φᶜⁿᶜ = dfm(φ)

mean_fields = (; uᶠⁿᶜ, vᶜⁿᶜ, wᶜⁿᶠ, bᶜⁿᶜ, φᶜⁿᶜ)
@info "Created mean fields"
# Correlation fields
#  Self
u′u′ᶠⁿᶜ = dfm(KernelFunctionOperation{Face, Center, Center}(variance_kernel_func, grid, u, uᶠⁿᶜ))
v′v′ᶜⁿᶜ = dfm(KernelFunctionOperation{Center, Face, Center}(variance_kernel_func, grid, v, vᶜⁿᶜ))
w′w′ᶜⁿᶠ = dfm(KernelFunctionOperation{Center, Center, Face}(variance_kernel_func, grid, w, wᶜⁿᶠ))

#  Cross terms
u′v′ᶜⁿᶜ = dfm(KernelFunctionOperation{Center, Center, Center}(u′v′ᶜᶜᶜ_kernel_func, grid, u, uᶠⁿᶜ, v, vᶜⁿᶜ))
v′w′ᶜⁿᶜ = dfm(KernelFunctionOperation{Center, Center, Center}(v′w′ᶜᶜᶜ_kernel_func, grid, v, vᶜⁿᶜ, w, wᶜⁿᶠ))
w′u′ᶜⁿᶜ = dfm(KernelFunctionOperation{Center, Center, Center}(w′u′ᶜᶜᶜ_kernel_func, grid, w, wᶜⁿᶠ, u, uᶠⁿᶜ))

#  Buoyancy
u′b′ᶜⁿᶜ = dfm(KernelFunctionOperation{Center, Center, Center}(u′b′ᶜᶜᶜ_kernel_func, grid, u, uᶠⁿᶜ, b, bᶜⁿᶜ))
v′b′ᶜⁿᶜ = dfm(KernelFunctionOperation{Center, Center, Center}(v′b′ᶜᶜᶜ_kernel_func, grid, v, vᶜⁿᶜ, b, bᶜⁿᶜ))
w′b′ᶜⁿᶜ = dfm(KernelFunctionOperation{Center, Center, Center}(w′b′ᶜᶜᶜ_kernel_func, grid, w, wᶜⁿᶠ, b, bᶜⁿᶜ))

correlation_fields = (; u′u′ᶠⁿᶜ, v′v′ᶜⁿᶜ, w′w′ᶜⁿᶠ, u′v′ᶜⁿᶜ, v′w′ᶜⁿᶜ, w′u′ᶜⁿᶜ, u′b′ᶜⁿᶜ, v′b′ᶜⁿᶜ, w′b′ᶜⁿᶜ)
@info "Created correlation fields"
# Now we have defined all our operations, we loop over the frames

# see #2024
output_path = "$foldername/dfm.jld2"
println()
for frame in frames
    # Get the data
    print("Calculating... $frame\r")
    update_fields!(input_fields, frame, "$foldername/$filename.jld2")
    # Compute
    map(compute!, mean_fields)
    map(compute!, correlation_fields)
    
    # Save to file
    write_output(mean_fields, correlation_fields, frame, output_path)
end
@info "Writing grid"
write_grid_times(grid, frames, ts, output_path)
@info "Finished!"
