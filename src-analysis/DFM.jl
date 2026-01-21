# The down-front (y) average of each field
include("terms/terms.jl")

u_dfm = dfm(input_fields.u)
v_dfm = dfm(input_fields.v)
w_dfm = dfm(input_fields.w)
b_dfm = dfm(input_fields.b)

ψ = Field(CumulativeIntegral(-u_dfm; dims=3))

skip_update = (:pNHS, :u_next, :v_next, :w_next, :b_next, :pNHS_next)

dependency_fields = (; u_dfm, v_dfm, w_dfm, b_dfm, ψ)
output_fields = dependency_fields
