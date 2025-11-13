# Shortcut for filtering...
@inline filt(a, σ::Tuple) = imfilter(a, gaussian(σ))
@inline filt(a, σ...) = filt(a, σ)

# If only one σ, repeat for length of a
@inline filt(a::A, σ::Integer) where {T, n, A<:Array{T, n}} = filt(a, (repeat([σ], n)..., ))
