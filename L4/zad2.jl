# Paweł Polerowicz 254626

function warNewton(x::Vector{Float64}, fx::Vector{Float64}, t::Float64)
    n = length(x)
    w = fx[n]
    for i in n-1:-1:1
        w = fx[i] + (t - x[i])*w
    end
    nt = w
    return nt
end
