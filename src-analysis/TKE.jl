using Oceananigans.Operators

@inline dfm(a) = Field(Average(a; dims=2))

u_dfm = dfm(u)
v_dfm = dfm(v)
w_dfm = dfm(w)

mean_fields = (; u_dfm, v_dfm, w_dfm)

function KE_func(i, j, k, grid, u, v, w)
    
    uᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u)
    vᶜᶜᶜ = ℑyᵃᶜᵃ(i, j, k, grid, v)
    wᶜᶜᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w)
    
    KE = (
        uᶜᶜᶜ^2
        + vᶜᶜᶜ^2
        + wᶜᶜᶜ^2
    ) / 2
    
    return KE
end

function TKE_func(i, j, k, grid, u, v, w)
    U, V, W = 0, 0, 0
    KE = 0
    for _j in 1:grid.Ny
        
        uᶜᶜᶜ = ℑxᶜᵃᵃ(i, _j, k, grid, u)
        vᶜᶜᶜ = ℑyᵃᶜᵃ(i, _j, k, grid, v)
        wᶜᶜᶜ = ℑzᵃᵃᶜ(i, _j, k, grid, w)
        
        U += uᶜᶜᶜ
        V += vᶜᶜᶜ
        W += wᶜᶜᶜ
        
        KE += uᶜᶜᶜ^2
        KE += vᶜᶜᶜ^2
        KE += wᶜᶜᶜ^2
    end
    U = U / grid.Ny
    V = V / grid.Ny
    W = W / grid.Ny
    TKE = KE / (2grid.Ny) - (U^2 + V^2 + W^2) / 2
    
    return TKE
end


#=
function MKE_func(i, j, k, grid, u_dfm, v_dfm, w_dfm)
    
    uᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, u_dfm)
    vᶜᶜᶜ = v_dfm
    wᶜᶜᶜ = ℑzᵃᵃᶜ(i, j, k, grid, w_dfm)
    
    MKE = (
        uᶜᶜᶜ^2
        + vᶜᶜᶜ^2
        + wᶜᶜᶜ^2
    ) / 2
    
    return MKE
end

KE_op = KernelFunctionOperation{Center, Center, Center}(KE_func, grid, u, v, w)

KE = dfm(KE_op)

MKE_op = KernelFunctionOperation{Center, Nothing, Center}(MKE_func, grid, u_dfm, v_dfm, w_dfm)

TKE = Field(KE - MKE_op)
=#

TKE_op = KernelFunctionOperation{Center, Nothing, Center}(TKE_func, grid, u, v, w)
TKE = Field(TKE_op)
outputs = (; TKE)

#update_outputs!(outputs) = map(compute!, (; mean_fields, outputs))

function update_outputs!(outputs)
    #map(compute!, mean_fields)
    #compute!(KE)
    map(compute!, outputs)
end