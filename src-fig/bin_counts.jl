function bin_counts(arr, bins)
    ns = zeros(Integer, length(bins) - 1)
    for i in 1:length(ns)
        ns[i] = count(arr) do a
            bins[i] <= a < bins[i+1]
        end
    end
    return ns
end