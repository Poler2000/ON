# Paweł Polerowicz 254626

function ilorazyRoznicowe(x::Vector{Float64}, f::Vector{Float64})
    n = length(x)
    fx = Vector{Float64}(undef, n)
    for i in 1:n
        fx[i] = f[i]
    end
    for j in 1:n
        for i in n:-1:j+1
            fx[i] = (fx[i] - fx[i - 1]) / (x[i] - x[i - j])
        end
    end
    return fx
end
