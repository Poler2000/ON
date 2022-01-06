# Paweł Polerowicz 254626

using SparseArrays
using LinearAlgebra

function luPartialPivoting(A::Dict{Int64, Float64}, b::Vector{Float64}, n::Int64, l::Int64)
    P = zeros(Int64, n + 1)

    for i in 1:(n+1)
        P[i] = i
    end

    for i in 1:n
        maxA = 0.0
        i_max = i
        upper = min(n, i + l + 1)

        for k in i:upper
            absA = abs(get(A, k + (i - 1) * n, 0.0))
            if (absA > maxA)
                maxA = absA
                i_max = k
            end
        end
        if (maxA < eps(Float64))
            println("matrix is degenerate")
            return -1
        end
        upper = min(n, i + 2 * l + 1)
        if (i_max != i)
            j = P[i]
            P[i] = P[i_max]
            P[i_max] = j

            # O(n) * n

            for j in 1:upper
                tmp = get(A, i + (j - 1) * n, 0.0)
                if (tmp != 0 || get(A, (i_max + (j - 1) * n), 0.0) != 0)
                    A[i + (j - 1) * n] = get(A, (i_max + (j - 1) * n), 0.0)
                    A[i_max + (j - 1) * n] = tmp
                end
            end

            P[n+1] += 1
        end

        for j in i+1:upper
            A[j + (i - 1) * n] = get(A, (j + (i - 1) * n), 0.0) / get(A, (i + (i - 1) * n), 0.0)
            for k in i+1:upper
                A[j + (k - 1) * n] = get(A, (j + (k - 1) * n), 0.0) - get(A, (j + (i - 1) * n), 0.0) * get(A, (i + (k - 1) * n), 0.0)
            end
        end

        #=ii = P[i]
        for j in i+1:n
            jj = P[j]
            A[jj + (i - 1) * n] /= get(A, (ii + (i - 1) * n), 0.0)
            for k in i+1:n
                jj = P[j]
                A[jj + (k - 1) * n] -= get(A, (jj + (i - 1) * n), 0.0) * get(A, (ii + (k - 1) * n), 0.0)
            end
        end=#
    end
    println("Solving!")
    x = zeros(Float64, n)
    for i in 1:n
        x[i] = b[P[i]]

        for k in 1:i-1
            x[i] -= get(A, i + (k - 1) * n, 0.0) * x[k]
        end
    end
    for i in n:-1:1
        for k in i+1:n
            x[i] -= get(A, i + (k - 1) * n, 0.0) * x[k]
        end

        x[i] /= get(A, i + (i - 1) * n, 0.0)
    end
    println(x)
    #=for j in 1:n
        A[j + (j - 1) * n] = 1
    end=#

end

function luNoPivoting(A::Dict{Int64, Float64}, b::Vector{Float64}, n::Int64, l::Int64)
    for i in 1:n
        upper = min(n, i + l + 1)
        for j in (i+1):upper
            A[j + (i - 1) * n] = get(A, j + (i - 1) * n, 0.0) / get(A, i + (i - 1) * n, 0.0)

            for k in (i+1):upper
                A[j + (k - 1) * n] = get(A, j + (k - 1) * n, 0.0) - get(A, j + (i - 1) * n, 0.0) * get(A, i + (k - 1) * n, 0.0)
            end
        end
    end

    y = zeros(Float64, n)
    for i in 1:n
        sum = 0;
        lower = max(1, floor(Int64, (i - 1) / l) * l - 1)
        for k in lower:(i-1)
            sum = sum + get(A, i + (k - 1) * n, 0.0) * y[k]
        end
        y[i] = b[i] - sum
    end

    x = zeros(Float64, n)
    for i in n:-1:1
        sum = 0;

        upper = min(n, i + l + 1)
        for k in (i+1):upper
            sum = sum + get(A, i + (k - 1) * n, 0.0) * x[k]
        end
        x[i] = (1 / get(A, i + (i - 1) * n, 0.0)) * (y[i] - sum)
    end

    println(x)
end

function helper1(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, n::Int64, l::Int64, i, P)
    maxA = 0.0
    i_max = i
    upper = min(n, i + l + 1)

    for k in i:upper
        absA = abs(A[k, i])
        if (absA > maxA)
            maxA = absA
            i_max = i
        end
    end
    if (maxA < eps(Float64))
        println("matrix is degenerate")
        return -1
    end

    if (i_max != i)
        j = P[i]
        P[i] = P[i_max]
        P[i_max] = j

        # O(n) * n
        for j in 1:upper
            tmp = A[i, j]
            A[i, j] = A[i_max, j]
            A[i_max, j] = tmp
        end

        P[n+1] += 1
    end
end

function helper1M(A::Dict{Int64, Float64}, b::Vector{Float64}, n::Int64, l::Int64, i, P)
    maxA = 0.0
    i_max = i
    upper = min(n, i + l + 1)

    for k in i:n
        absA = abs(get(A, k + (i - 1) * n, 0.0))
        if (absA > maxA)
            maxA = absA
            i_max = k
        end
    end
    if (maxA < eps(Float64))
        println("matrix is degenerate")
        return -1
    end

    if (i_max != i)
        j = P[i]
        P[i] = P[i_max]
        P[i_max] = j

        # O(n) * n
        for j in 1:n
            tmp = get(A, i + (j - 1) * n, 0.0)
            A[i + (j - 1) * n] = get(A, (i_max + (j - 1) * n), 0.0)
            A[i_max + (j - 1) * n] = tmp
        end

        P[n+1] += 1
    end
end

function helper2(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, n::Int64, l::Int64, i)
    # O(c) * n
    upper = min(n, i + l + 1)

    for j in (i+1):n
        A[j, i] = A[j, i] / A[i, i]

        for k in (i+1):n
            A[j, k] = A[j, k] - A[j, i] * A[i, k]
        end
    end

end



function helper2M(A::Dict{Int64, Float64}, b::Vector{Float64}, n::Int64, l::Int64, i)
    upper = min(n, i + l + 1)

    for j in (i+1):n
        A[j + (i - 1) * n] = get(A, j + (i - 1) * n, 0.0) / get(A, i + (i - 1) * n, 0.0)

        for k in (i+1):n
            A[j + (k - 1) * n] = get(A, j + (k - 1) * n, 0.0) - get(A, j + (i - 1) * n, 0.0) * get(A, i + (k - 1) * n, 0.0)
        end
    end
end

function luDecompositionM(A::Dict{Int64, Float64}, b::Vector{Float64}, n::Int64, l::Int64)
    ops = 0
    P = zeros(Int64, n + 1)
    t1 = 0
    t2 = 0
    for i in 1:(n+1)
        P[i] = i
    end
    for i in 1:n
        t1 += @elapsed helper1M(A, b, n, l, i, P)
        t2 += @elapsed helper2M(A, b, n, l, i)
    end

    luSolveM(A, b, n, l, P)
    #=for i in 1:n
        for j in 1:n
            print(A[i, j])
            print(" ")
        end
        println(" ")
    end=#
    println("main part done")
    #println(ops)
    println(t1)
    println(t2)



    #luSolve(A, b, n, l)
end

function luDecomposition(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, n::Int64, l::Int64)
    ops = 0
    P = zeros(Int64, n + 1)
    t1 = 0
    t2 = 0
    for i in 1:(n+1)
        P[i] = i
    end
    for i in 1:n
        t1 += @elapsed helper1(A, b, n, l, i, P)
        t2 += @elapsed helper2(A, b, n, l, i)
    end
    #=for i in 1:n
        for j in 1:n
            print(A[i, j])
            print(" ")
        end
        println(" ")
    end=#
    println("main part done")
    #println(ops)
    println(t1)
    println(t2)



    #luSolve(A, b, n, l)
end

function luSolveM(A::Dict{Int64, Float64}, b::Vector{Float64}, n::Int64, l::Int64, P)
    x = zeros(Float64, n)
    for i in 1:n
        x[i] = b[P[i]]

        for k in 1:i
            x[i] -= get(A, i + (k - 1) * n, 0.0) * x[k]
        end
    end
    for i in n:-1:1
        for k in i+1:n
            x[i] -= get(A, i + (k - 1) * n, 0.0) * x[k]
        end

        x[i] /= get(A, i + (i - 1) * n, 0.0)
    end
    #=v = n / l
    y = zeros(Float64, n)
    ops = 0
    for i in 1:n
        sum = 0;
        lower = max(1, floor(Int64, (i - 1) / l) * l - 1)
        #println(lower)
        for k in lower:(i-1)
            sum = sum + get(A, i + (k - 1) * n, 0.0) * y[k]
            ops += 1
        end
        y[i] = b[i] - sum
        ops += 1
    end

    x = zeros(Float64, n)
    for i in n:-1:1
        sum = 0;
        ops += 1

        upper = min(n, i + l + 1)
        for k in (i+1):upper
            sum = sum + get(A, i + (k - 1) * n, 0.0) * x[k]
            ops += 1
        end
        x[i] = (1 / get(A, i + (i - 1) * n, 0.0)) * (y[i] - sum)
        ops += 1
    end=#
    println("second part done")
    #println(ops)
    println(x)
end

function luSolve(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, n::Int64, l::Int64)
    v = n / l
    y = zeros(Float64, n)
    ops = 0
    for i in 1:n
        sum = 0;
        lower = max(1, floor(Int64, (i - 1) / l) * l - 1)
        #println(lower)
        for k in lower:(i-1)
            sum = sum + A[i, k] * y[k]
            ops += 1
        end
        y[i] = b[i] - sum
        ops += 1
    end

    x = zeros(Float64, n)
    for i in n:-1:1
        sum = 0;
        ops += 1

        upper = min(n, i + l + 1)
        for k in (i+1):upper
            sum = sum + A[i, k] * x[k]
            ops += 1
        end
        x[i] = (1 / A[i, i]) * (y[i] - sum)
        ops += 1
    end
    println("second part done")
    #println(ops)
    #println(x)
end

function luTestM(A::Dict{Int64, Float64}, b::Vector{Float64}, n::Int64, l::Int64)
    @time luDecompositionM(A, b, n, l)
    #@time luSolveM(A, b, n, l)
end

function luTest(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, n::Int64, l::Int64)
    @time luDecomposition(A, b, n, l)
    @time luSolve(A, b, n, l)
end


f = open("Dane10000_1_1/A.txt", "r")

args = readline(f)
nl = parse.(Int64, split(args))
n = nl[1]
l = nl[2]
A = spzeros(Float64, n, n)
args = readlines(f)
A_Map = Dict{Int64, Float64}()

for arg in args
    entry = parse.(Float64, split(arg))
    x = floor(Int32, entry[1])
    y = floor(Int32, entry[2])
    v = entry[3]
    A[x,y] = v
    A_Map[x + (y-1) * n] = v
end

g = open("Dane10000_1_1/b.txt", "r")
m = readline(g)

b = ones(Float64, n)

for i in 1:n
    line = readline(g)
    v = parse(Float64, line)
    b[i] = v
    #A[i,n+1] = b[i]
end

C = spzeros(Float64, n, n)
for i in 1:n
    for j in 1:n
        if (A[i, j] != 0)
            C[i,j] = A[i,j]
        end
    end
end
#luTestM(A_Map, b, n, l)
#luTest(A, b, n, l)

println("Hello")
#luNoPivoting(A_Map, b, n, l)
luPartialPivoting(A_Map, b, n, l)

#@time lu(A)
