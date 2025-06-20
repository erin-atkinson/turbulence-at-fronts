# The down-front (y) average of each field

u_dfm = dfm(fields.u)
v_dfm = dfm(fields.v)
w_dfm = dfm(fields.w)
v_dfm = dfm(fields.b)

dependency_fields = (; u_dfm, v_dfm, w_dfm, b_dfm)
output_fields = dependency_fields
