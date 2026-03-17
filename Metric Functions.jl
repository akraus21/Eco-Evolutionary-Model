# --- Metric Functions --- 

function geometric_mean_ic50(bins, counts)
    s = 0.0
    n = 0.0
    @inbounds @simd for i in eachindex(bins, counts)
        c = copy(counts[i])
        if c > 0
            s += log10(bins[i]) * c
            n += c
        end
    end
    return n == 0 ? 0 : 10.0^(s / n)
end

function shannon_diversity(counts)
    p = vec(sum(counts, dims=2))
    p = p[p .> 0]
    p = p / sum(p)
    return -sum(p .* log2.(p))
end

function count_establishments(prev_counts, counts; threshold=50)
    est = 0
    for i in eachindex(prev_counts)
        if prev_counts[i] < threshold && counts[i] >= threshold
            est += 1
        end
    end
    return est
end

function count_extinctions(prev_counts, counts)
    ext = 0
    for i in eachindex(prev_counts)
        if prev_counts[i] > 0 && counts[i] == 0
            ext += 1
        end
    end
    return ext
end



