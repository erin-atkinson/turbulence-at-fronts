# The down-front (y) average of each field
include("terms/terms.jl")

u_dfm = dfm(rawfields.u)
v_dfm = dfm(rawfields.v)
w_dfm = dfm(rawfields.w)
b_dfm = dfm(rawfields.b)

dependency_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)
output_fields = dependency_fields
