# The down-front (y) average of each field
include("terms/terms.jl")

u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)

dependency_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)
output_fields = dependency_fields
