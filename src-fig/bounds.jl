for i in [:x, :y, :z], j in [:ᶜ, :ᶠ]
    fname = Symbol(i, j, :bounds)
    nodes = Symbol(i, :s, j)
    @eval function $fname(file, args::Tuple)
        xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ = grid_nodes(file)
        map(arg->argmin(abs.($nodes .- arg)), args)
    end
    @eval $fname(file, args...) = $fname(file, args)
    @eval $fname(file, i) = $fname(file, (i, ))[1]
    @eval $fname(filename::String, args...) = jldopen(file->$fname(file, args...), filename)
end

function tbounds(file, args::Tuple)
    ts = times(file)
    map(arg->argmin(abs.(times .- arg)), args)
end
tbounds(file, args...) = tbounds(foldername, args)
tbounds(file, i) = tbounds(foldername, (i, ))[1]
tbounds(filename::String, args...) = jldopen(file->tbounds(file, args...), filename)