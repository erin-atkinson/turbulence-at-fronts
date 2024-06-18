#= 
balance.jl
    Computes the down-front (y) mean and first-order correlations of the simulation output

    Arguments
        [1]: Folder containing simulation output, will also output to this folder
    Pass options in last argument
        i: use initialisation run
=#
ENV["JULIA_SCRATCH_TRACK_ACCESS"] = 0
include("analysis_functions.jl")
include("balance_functions.jl")

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
ν = FieldTimeSeries("$foldername/$filename.jld2", "νₑ"; iterations=frames[1])[1]
κ = FieldTimeSeries("$foldername/$filename.jld2", "κₑ"; iterations=frames[1])[1]

input_fields = (; u, v, w, b, φ, ν, κ)
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
# Coriolis
fvᶠⁿᶜ = Field(@at (Face, Nothing, Center) sp.f * vᶜⁿᶜ)
fuᶜⁿᶜ = Field(@at (Center, Nothing, Center) sp.f * uᶠⁿᶜ)

# Mean advection
uu_xᶠⁿᶜ = Field(@at (Face, Nothing, Center) uᶠⁿᶜ * ∂x(uᶠⁿᶜ))
wu_zᶠⁿᶜ = Field(@at (Face, Nothing, Center) wᶜⁿᶠ * ∂z(uᶠⁿᶜ))

uv_xᶜⁿᶜ = Field(@at (Center, Nothing, Center) uᶠⁿᶜ * ∂x(vᶜⁿᶜ))
wv_zᶜⁿᶜ = Field(@at (Center, Nothing, Center) wᶜⁿᶠ * ∂z(vᶜⁿᶜ))

uw_xᶜⁿᶠ = Field(@at (Center, Nothing, Face) uᶠⁿᶜ * ∂x(wᶜⁿᶠ))
ww_zᶜⁿᶠ = Field(@at (Center, Nothing, Face) wᶜⁿᶠ * ∂z(wᶜⁿᶠ))

ub_xᶜⁿᶜ = Field(@at (Center, Nothing, Center) uᶠⁿᶜ * ∂x(bᶜⁿᶜ))
wb_zᶜⁿᶜ = Field(@at (Center, Nothing, Center) wᶜⁿᶠ * ∂z(bᶜⁿᶜ))
@info "Created mean advection fields"
u′ = u - uᶠⁿᶜ
# We want v to be on Center to prevent interpolation from borking
v′ = @at((Center, Center, Center), v) - vᶜⁿᶜ
w′ = w - wᶜⁿᶠ
b′ = b - bᶜⁿᶜ

# Turbulence
u′u′_xᶠⁿᶜ = Field(@at (Face, Nothing, Center) Average(∂x(u′ * u′); dims=2))
w′u′_zᶠⁿᶜ = Field(@at (Face, Nothing, Center) Average(∂z(w′ * u′); dims=2))

u′v′_xᶜⁿᶜ = Field(@at (Center, Nothing, Center) Average(∂x(u′ * v′); dims=2))
w′v′_zᶜⁿᶜ = Field(@at (Center, Nothing, Center) Average(∂z(w′ * v′); dims=2))

u′w′_xᶜⁿᶠ = Field(@at (Center, Nothing, Face) Average(∂x(u′ * w′); dims=2))
w′w′_zᶜⁿᶠ = Field(@at (Center, Nothing, Face) Average(∂z(w′ * w′); dims=2))

u′b′_xᶜⁿᶜ = Field(@at (Center, Nothing, Center) Average(∂x(u′ * b′); dims=2))
w′b′_zᶜⁿᶜ = Field(@at (Center, Nothing, Center) Average(∂z(w′ * b′); dims=2))
@info "Created turbulence fields"
# Pressure
φ_xᶠⁿᶜ = Field(∂x(φᶜⁿᶜ))
φ_zᶜⁿᶠ = Field(∂z(φᶜⁿᶜ))
@info "Created pressure fields"
# Buoyancy
bᶜⁿᶠ = Field(@at (Center, Nothing, Face) bᶜⁿᶜ)
@info "Created buoyancy field"

@inline ∇_dot_b∇a(a, b) = ∂x(b * ∂x(a)) + ∂z(b * ∂z(a))

∇_dot_ν∇u = Field(@at (Face, Nothing, Center) Average(∇_dot_b∇a(u, ν); dims=2))
∇_dot_ν∇v = Field(@at (Center, Nothing, Center) Average(∇_dot_b∇a(v, ν); dims=2))
∇_dot_ν∇w = Field(@at (Center, Nothing, Face) Average(∇_dot_b∇a(w, ν); dims=2))

∇_dot_κ∇b = Field(@at (Center, Nothing, Center) Average(∇_dot_b∇a(b, κ); dims=2))

@info "Created sgs fields"
u_balance = (; fvᶠⁿᶜ, uu_xᶠⁿᶜ, wu_zᶠⁿᶜ, u′u′_xᶠⁿᶜ, w′u′_zᶠⁿᶜ, φ_xᶠⁿᶜ, ∇_dot_ν∇u)
v_balance = (; fuᶜⁿᶜ, uv_xᶜⁿᶜ, wv_zᶜⁿᶜ, u′v′_xᶜⁿᶜ, w′v′_zᶜⁿᶜ, ∇_dot_ν∇v)
w_balance = (; bᶜⁿᶠ, uw_xᶜⁿᶠ, ww_zᶜⁿᶠ, u′w′_xᶜⁿᶠ, w′w′_zᶜⁿᶠ, φ_zᶜⁿᶠ, ∇_dot_ν∇w)
b_balance = (; ub_xᶜⁿᶜ, wb_zᶜⁿᶜ, u′b′_xᶜⁿᶜ, w′b′_zᶜⁿᶜ, ∇_dot_κ∇b)

@info "Created correlation fields"
# Now we have defined all our operations, we loop over the frames

# see #2024
output_path = "$foldername/balance.jld2"
println()
for frame in frames
    # Get the data
    print("Updating... $frame\r")
    update_fields!(input_fields, frame, "$foldername/$filename.jld2")
    # Compute
    print("Calculating mean... $frame\r")
    map(compute!, mean_fields)
    
    print("Calculating u balance... $frame\r")
    map(compute!, u_balance)
    
    print("Calculating v balance... $frame\r")
    map(compute!, v_balance)
    
    print("Calculating w balance... $frame\r")
    map(compute!, w_balance)
    
    print("Calculating b balance... $frame\r")
    map(compute!, b_balance)
    
    # Save to file
    write_output(u_balance, v_balance, w_balance, b_balance, frame, output_path)
end
@info "Writing grid"
write_grid_times(grid, frames, ts, output_path)
@info "Finished!"
