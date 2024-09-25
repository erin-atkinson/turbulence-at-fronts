@inline dfm(a) = Field(Average(a; dims=2))

δ = dfm(∂x(u))
ζ = dfm(∂x(v))

outputs = (; δ, ζ)

function update_outputs!(outputs)
    map(compute!, outputs)
end