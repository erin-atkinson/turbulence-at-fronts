include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/diffusion.jl")

include("terms/bbalance/advection.jl")
include("terms/bbalance/turbulence.jl")
include("terms/bbalance/laplacian.jl")


u_dfm = dfm(fields.u)
v_dfm = dfm(fields.v)
w_dfm = dfm(fields.w)
b_dfm = dfm(fields.b)

mean_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)

# Turbulence
turbulence_h = Field(KernelFunctionOperation{Center, Center, Center}(turbulence_h_func, grid, clock, fields, mean_fields, sp))
turbulence_z = Field(KernelFunctionOperation{Center, Center, Center}(turbulence_z_func, grid, clock, fields, mean_fields, sp))

adv_x = Field(KernelFunctionOperation{Center, Center, Center}(adv_x_func, grid, clock, fields, mean_fields, sp))
adv_α = Field(KernelFunctionOperation{Center, Center, Center}(adv_α_func, grid, clock, fields, mean_fields, sp))
adv_y = Field(KernelFunctionOperation{Center, Center, Center}(adv_y_func, grid, clock, fields, mean_fields, sp))
adv_z = Field(KernelFunctionOperation{Center, Center, Center}(adv_z_func, grid, clock, fields, mean_fields, sp))

laplacian = Field(KernelFunctionOperation{Center, Nothing, Center}(laplacian_func, grid, clock, fields, mean_fields, sp))

turbulence_h_dfm = dfm(turbulence_h)
turbulence_z_dfm = dfm(turbulence_z)

adv_x_dfm = dfm(adv_x)
adv_α_dfm = dfm(adv_α)
adv_y_dfm = dfm(adv_y)
adv_z_dfm = dfm(adv_z)

b_fields = (;
    laplacian,
    turbulence_h,
    turbulence_z,
    adv_x,
    adv_α,
    adv_y,
    adv_z,
    turbulence_h_dfm,
    turbulence_z_dfm,
    adv_x_dfm,
    adv_α_dfm,
    adv_y_dfm,
    adv_z_dfm
)

dependency_fields = merge(mean_fields, b_fields)

output_fields = (;
    laplacian,
    turbulence_h_dfm,
    turbulence_z_dfm,
    adv_x_dfm,
    adv_α_dfm,
    adv_y_dfm,
    adv_z_dfm
)

frames = 801:1001
