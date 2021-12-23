# Paweł Polerowicz 254626

# helper function
function wektorRazyDwumian(x::Vector{Float64}, x1)
    n = length(x)
    a = Vector{Float64}(undef, n + 1)
    a[1] = 0
    for i in 2:n+1
        a[i] = x[i - 1]
    end
    for i in 1:n
        a[i] = a[i] + x[i] * x1
    end
    return a
end

function naturalna(x::Vector{Float64}, fx::Vector{Float64})
    n = length(x)
    # array of natural form factors
    a = Vector{Float64}(undef, n)
    # helper polynomial
    w = Vector{Float64}(undef, 1)
    w[1] = 1
    for i in 1:n
        a[i] = fx[i]
    end
    for i in 1:n-1
        w = wektorRazyDwumian(w, -x[i])
        for j in 1:i
            a[j] += fx[i+1] * w[j]
        end
    end
    return a
end
