# Paweł Polerowicz 254626

using SparseArrays
using LinearAlgebra

function swapRows(A::Matrix{Float64}, h::Int64, i_max::Int64, n::Int64)
    for i in 1:(n+1)
        tmp = A[h, i]
        A[h, i] = A[i_max, i]
        A[i_max, i] = tmp
    end
end

function gaussPartialPivoting(A::Dict{Int64, Float64},
    b::Vector{Float64}, n::Int64, l::Int64)
    k = 1

    while k <= n
        i_max = k
        v_max = abs(get(A, i_max + (k - 1) * n, 0.0))
        upper = min(n, k + l + 1)

        for i in (k+1):upper
            if (abs(get(A, i + (k - 1) * n, 0.0)) > v_max)
                i_max = i
                v_max = abs(get(A, i + (k - 1) * n, 0.0))
            end
        end
        if (v_max < eps(Float64))
            println("matrix is degenerate")
            return -1
        end

        upper = min(n, k + 2 * l)
        lower = k - l
        if (i_max != k)
            tmp = b[k]
            b[k] = b[i_max]
            b[i_max] = tmp

            for j in lower:upper
                tmp = get(A, k + (j - 1) * n, 0.0)
                if (tmp != 0 || get(A, (i_max + (j - 1) * n), 0.0) != 0)
                    A[k + (j - 1) * n] = get(A, (i_max + (j - 1) * n), 0.0)
                    A[i_max + (j - 1) * n] = tmp
                end
            end
        end

        for i in (k+1):upper
            f = get(A, i + (k - 1) * n, 0.0) / get(A, k + (k - 1) * n, 0.0)
            A[i + (k - 1) * n] = 0
            for j in (k+1):upper
                A[i + (j - 1) * n] = get(A, i + (j - 1) * n, 0.0) - get(A, k + (j - 1) * n, 0.0) * f
            end
            b[i] -= b[k] * f
        end

        k = k + 1
    end
    x = zeros(Float64, n)
    for i in n:-1:1
        upper = min(n, i + 2 * l)
        x[i] = b[i]
        for j in (i+1):upper
            x[i] = x[i] - get(A, i + (j - 1) * n, 0.0)*x[j]
        end
        x[i] = x[i]/get(A, i + (i - 1) * n, 0.0)
    end
    return x
end

function gaussNoPivoting(A::Dict{Int64, Float64},
    b::Vector{Float64}, n::Int64, l::Int64)
    k = 1

    while k <= n
        i_max = k
        v_max = abs(get(A, i_max + (k - 1) * n, 0.0))

        if (v_max < eps(Float64))
            println("matrix is degenerate")
            return -1
        end

        upper = min(n, k + 2 * l)
        lower = k - l

        for i in (k+1):upper
            f = get(A, i + (k - 1) * n, 0.0) / get(A, k + (k - 1) * n, 0.0)
            A[i + (k - 1) * n] = 0
            for j in (k+1):upper
                A[i + (j - 1) * n] = get(A, i + (j - 1) * n, 0.0) - get(A, k + (j - 1) * n, 0.0) * f
            end
            b[i] -= b[k] * f
        end

        k = k + 1
    end
    x = zeros(Float64, n)
    for i in n:-1:1
        upper = min(n, i + 2 * l)
        x[i] = b[i]
        for j in (i+1):upper
            x[i] = x[i] - get(A, i + (j - 1) * n, 0.0)*x[j]
        end
        x[i] = x[i]/get(A, i + (i - 1) * n, 0.0)
    end
    return x
end

function gaussElimination(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, n::Int64)
    k = 1
    A_star = zeros(Float64, n, n + 1)
    A_star[1:n, 1:n] = A[1:n, 1:n]
    A_star[1:n, n + 1] = b[1:n]

    # n times
    while k <= n
        i_max = k
        v_max = abs(A[i_max, k])
        # O(n)
        for i in (k+1):n
            if (abs(A_star[i, k]) > v_max)
                i_max = i
                v_max = abs(A_star[i,k])
            end
        end
        if (A_star[i_max, k] == 0)
            k = k + 1
            continue
        else
            if (i_max != k)
                swapRows(A_star, k, i_max, n)
            end

            for i in (k+1):n
                f = A_star[i, k] / A_star[k, k]
                # Fill with zeros the lower part of pivot column: */
                A_star[i, k] = 0
                # Do for all remaining elements in current row: */
                for j in (k+1):(n+1)
                    A_star[i, j] = A_star[i, j] - A_star[k, j] * f
                end
            end
        end
        k = k + 1
    end

    for i in 1:n
        for j in 1:n
            print(A_star[i, j])
            print(" ")
        end
        println(" ")
    end

    x = zeros(Float64, n)
    for i in n:-1:1
        x[i] = A_star[i,n+1]
        for j in (i+1):n
            x[i] = x[i] - A_star[i,j]*x[j]
        end
        x[i] = x[i]/A_star[i,i]
    end
    return x
end

function gaussThomasMap(A::Dict{Int64, Float64},
    b::Vector{Float64}, n::Int64, l::Int64)
    v = floor(Int64, n / l)
    b = reshape(b,l,v)
    x = zeros(Float64, l, v)
    c = zeros(Float64, l, v)
    D = zeros(Float64, l, l, v)
    Q = zeros(Float64, l, l, v)
    G = zeros(Float64, l, l, v)
    C = zeros(Float64, l, l, v - 1)
    B = zeros(Float64, l, l, v - 1)
    for k in 1:(v - 1)
        for i in 1:l
            for j in 1:l
                D[i, j, k] = get(A, ((k-1)*l)+i + (((k-1)*l)+j) * n, 0.0)
                B[i, j, k] = get(A, ((k-1)*l)+i + l + n * (((k-1)*l) + j), 0.0)
                C[i, j, k] = get(A, ((k-1)*l)+i + n * (((k-1)*l)+j + l), 0.0)
            end
        end
    end

    for i in 1:l
        for j in 1:l
            D[i, j, v] = get(A, (v-1)*l + i + ((v-1)*l + j) * n, 0.0)
        end
    end

    Q[:, :, 1] = D[:, :, 1]
    G[:, :, 1] = Q[:, :, 1] \ C[:, :, 1]

    for i in 1:l
        for j in 1:l
            print(Q[i, j, 1])
            print(" ")
        end
        println(" ")
    end
    for i in 1:l
        for j in 1:l
            print(G[i, j, 1])
            print(" ")
        end
        println(" ")
    end

    for k in 2:(v-1)
        Q[:,:,k]=D[:,:,k]-B[:,:,k-1]*G[:,:,k-1]
        G[:, :, k] = Q[:, :, k] \ C[:, :, k]
    end

    Q[:, :, v] = D[:, :, v] - B[:, :, v - 1] * G[:, :, v-1]

    #=for k in 1:v
        for i in 1:l
            for j in 1:l
                print(G[i, j, k])
                print(" ")
            end
            println(" ")
        end
        println(" ")
    end=#

    c[:, 1] = gaussElimination(Q[:, :, 1], b[:, 1], l)


    for k in 2:v
        c[:, k] = Q[:, :, k] \ (b[:,k]- B[:,:,k-1] * c[:,k-1])
    end


    x[:, v] = c[:, v]

    for k in (v-1):-1:1
        x[:,k]=c[:,k]-G[:,:,k]*x[:,k+1];
    end
    println(x)
end

function loadA(inputFileA)
    f = open(inputFileA, "r")
    args = readline(f)
    nl = parse.(Int64, split(args))
    n = nl[1]
    l = nl[2]
    args = readlines(f)
    A_Map = Dict{Int64, Float64}()

    for arg in args
        entry = parse.(Float64, split(arg))
        x = floor(Int32, entry[1])
        y = floor(Int32, entry[2])
        v = entry[3]
        A_Map[x + (y-1) * n] = v
    end

    return A_Map, n, l
end

function loadb(inputFileb)
    g = open(inputFileb, "r")
    mstr = readline(g)
    m = parse(Int64, mstr)
    b = ones(Float64, m)

    for i in 1:m
        line = readline(g)
        v = parse(Float64, line)
        b[i] = v
    end
    return b
end
#=
testRange = [16, 10000, 50000, 100000]

for i in eachindex(testRange)
    n = testRange[i]
    path = "Dane$n" * "_1_1/A.txt"
    A, n, l = loadA(path)
    path = "Dane$n" * "_1_1/b.txt"
    c = loadb(path)
    b = loadb(path)
    A_copy = Dict{Int64, Float64}()
    for (key, value) in A
        A_copy[key] = value
    end

    x = gaussNoPivoting(A, b, n, l)

    errors = Float64(0.0)

    for i in eachindex(x)
        errors += abs(1.0 - x[i])
    end
    print("$errors & ")
    errors /= n
    print("$errors & ")

    x = gaussPartialPivoting(A_copy, c, n, l)

    errors = Float64(0.0)

    for i in eachindex(x)
        errors += abs(1.0 - x[i])
    end
    print("$errors & ")
    errors /= n
    print("$errors")
    println(" ")
end=#

for i in 5000:5000:100000
    path = "A$i.txt"
    A, n, l = loadA(path)
    path = "b$i.txt"
    c = loadb(path)
    b = loadb(path)
    A_copy = Dict{Int64, Float64}()
    for (key, value) in A
        A_copy[key] = value
    end

    t1 = @elapsed gaussNoPivoting(A, b, n, l)
    t2 = @elapsed gaussPartialPivoting(A_copy, c, n, l)
    println("$i; $t1; $t2")
end

#=
A, n, l = loadA("Dane100000_1_1/A.txt")
b = loadb("Dane100000_1_1/b.txt")
println("to arms!")
@time gaussEliminationOptimized(A, b, n, l)=#
#=x = gaussEliminationOptimized(A, b, n, l)

errors = Float64(0.0)

for i in eachindex(x)
    global errors += abs(1.0 - x[i])
end
println(errors)
errors /= 100000
println(errors)
println(x)=#


#=
f = open("Dane16_1_1/A.txt", "r")
#f = open("A1000.txt", "r")

args = readline(f)
nl = parse.(Int64, split(args))
n = nl[1]
l = nl[2]
A = spzeros(Float64, n, n + 1)
A_Map = Dict{Int64, Float64}()
args = readlines(f)

for arg in args
    entry = parse.(Float64, split(arg))
    x = floor(Int32, entry[1])
    y = floor(Int32, entry[2])
    v = entry[3]
    #println(v)
    A[x,y] = v
    A_Map[x + (y - 1) * n] = v
    #println(A[x,y])
end

g = open("Dane16_1_1/b.txt", "r")
#g = open("Dane16_1_1/b.txt", "r")
m = readline(g)

b = ones(Float64, n)

for i in 1:n
    line = readline(g)
    v = parse(Float64, line)
    b[i] = v
    A[i,n+1] = b[i]
end
#println(A)


for i in 1:n
    for j in 1:n
        print(get(A_Map, i + (j - 1) * n, 0.0))
        print(" ")
    end
    println(" ")
end

#println(b)
#luTest(A, b, n, l)
#println(gaussEliminationOptimized(A_Map, b, n, l))
#println(gaussElimination(A, b, n))
#luDecompositionOptimized(A, b, n, l)
#@time gaussThomasMap(A_Map, b, n, l)
=#
