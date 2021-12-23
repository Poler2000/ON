# Paweł Polerowicz 254626

using SparseArrays
using LinearAlgebra

function swapRows(A::SparseMatrixCSC{Float64, Int64}, h::Int64, i_max::Int64, n::Int64)
    for i in 1:(n+1)
        tmp = A[h, i]
        A[h, i] = A[i_max, i]
        A[i_max, i] = tmp
    end
end

function swapRowsLazy(A::SparseMatrixCSC{Float64, Int64}, h::Int64, i_max::Int64, n::Int64, l::Int64)
    for i in 1:(n+1)
        tmp = A[h, i]
        A[h, i] = A[i_max, i]
        A[i_max, i] = tmp
    end
end

function gaussElimination(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, n::Int64)
    k = 1

    # n times
    while k <= n
        i_max = k
        v_max = abs(A[i_max, k])
        # O(n)
        for i in (k+1):n
            if (abs(A[i, k]) > v_max)
                i_max = i
                v_max = abs(A[i,k])
            end
        end
        if (A[i_max, k] == 0)
            k = k + 1
            continue
        else
            if (i_max != k)
                swapRows(A, k, i_max, n)
            end

            for i in (k+1):n
                f = A[i, k] / A[k, k]
                # Fill with zeros the lower part of pivot column: */
                A[i, k] = 0
                # Do for all remaining elements in current row: */
                for j in (k+1):(n+1)
                    A[i, j] = A[i, j] - A[k, j] * f
                end
            end
        end
        k = k + 1
    end
    x = zeros(Float64, n)
    for i in n:-1:1
        x[i] = A[i,n+1]
        for j in (i+1):n
            x[i] = x[i] - A[i,j]*x[j]
        end
        x[i] = x[i]/A[i,i]
    end
    return x
end

function gaussEliminationOptimized(A::SparseMatrixCSC{Float64, Int64},
    b::Vector{Float64}, n::Int64, l::Int64)
    k = 1
    v = n / l

    while k <= n
        i_max = k
        v_max = abs(A[i_max, k])
        upper = min(n, k + l)
        for i in (k+1):upper
            if (abs(A[i, k]) > v_max)
                i_max = i
                v_max = abs(A[i,k])
            end
        end
        if (A[i_max, k] == 0)
            k = k + 1
            continue
        else
            if (i_max != k)
                swapRowsLazy(A, k, i_max, n, l)
            end

            for i in (k+1):n
                f = A[i, k] / A[k, k]
                # Fill with zeros the lower part of pivot column: */
                A[i, k] = 0
                # Do for all remaining elements in current row: */
                upper = min(n, i + l + 2)
                for j in (k+1):upper
                    A[i, j] = A[i, j] - A[k, j] * f
                end
                A[i, n+1] = A[i, n+1] - A[k, n+1] * f
            end
        end
        k = k + 1
    end
    x = zeros(Float64, n)
    for i in n:-1:1
        x[i] = A[i,n+1]
        for j in (i+1):n
            x[i] = x[i] - A[i,j]*x[j]
        end
        x[i] = x[i]/A[i,i]
    end
    return x
end


function luDecomposition(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, n::Int64, l::Int64)
    lower = spzeros(Float64, n, n)
    upper = spzeros(Float64, n, n)
    P = zeros(Int64, n + 1)

    for i in 1:(n+1)
        P[i] = i
    end
    for i in 1:n
        maxA = 0.0
        i_max = i

        for k in i:n
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

            for j in 1:n
                tmp = A[i, j]
                A[i, j] = A[i_max, j]
                A[i_max, j] = tmp
            end

            P[n+1] += 1
        end


        for j in (i+1):n
            A[j, i] = A[j, i] / A[i, i]

            for k in (i+1):n
                A[j, k] = A[j, k] - A[j, i] * A[i, k]
            end
        end
    end

    println("lower")
    for i in 1:n
        for j in 1:n
            print(A[i, j])
            print(" ")
        end
        println(" ")
    end
    #=
    x = zeros(n)
    for i in 1:n
        x[i] = b[P[i]]

        for k in 1:(i-1)
            x[i] = x[i] - A[i, k] * x[k]
        end
    end

    for i in n:-1:1
        for k in (i+1):n
            x[i] = A[i, k] * x[k]
        end

        x[i] = x[i] / A[i, i]
    end

    println("Results:")
    println(x)
    =#

    y = zeros(Float64, n)
    for i in 1:n
        sum = 0;
        for k in 1:(i-1)
            sum = sum + A[i, k] * y[k]
        end
        y[i] = b[i] - sum
    end

    x = zeros(Float64, n)
    for i in n:-1:1
        sum = 0;
        for k in (i+1):n
            sum = sum + A[i, k] * x[k]
        end
        x[i] = (1 / A[i, i]) * (y[i] - sum)
    end

    println(x)
end


#=
function luDecomposition(A::SparseMatrixCSC{Float64, Int64}, b::Vector{Float64}, n::Int64, l::Int64)
    lower = spzeros(Float64, n, n)
    upper = spzeros(Float64, n, n)

    for i in 1:n
        for k in i:n
            sum = 0
            for j in 1:(i-1)
                sum = sum + (lower[i, j] * upper[j, k])
            end
            upper[i, k] = A[i, k] - sum
        end

        for k in i:n
            if (i == k)
                lower[i, i] = 1
            else
                sum = 0
                for j in 1:(i-1)
                    sum = sum + (lower[k, j] * upper[j, i])
                end
                lower[k, i] = (A[k, i] - sum) / upper[i, i]
            end
        end
    end
    println("lower")
    for i in 1:n
        for j in 1:n
            print(lower[i, j])
            print(" ")
        end
        println(" ")
    end
    println("upper")
    for i in 1:n
        for j in 1:n
            print(upper[i, j])
            print(" ")
        end
        println(" ")
    end
end
=#
f = open("Dane16_1_1/A.txt", "r")

args = readline(f)
nl = parse.(Int64, split(args))
n = nl[1]
l = nl[2]
A = spzeros(Float64, n, n + 1)
args = readlines(f)

for arg in args
    entry = parse.(Float64, split(arg))
    x = floor(Int32, entry[1])
    y = floor(Int32, entry[2])
    v = entry[3]
    #println(v)
    A[x,y] = v
    #println(A[x,y])
end

g = open("Dane16_1_1/b.txt", "r")
m = readline(g)

b = zeros(Float64, n)

for i in 1:n
    line = readline(g)
    v = parse(Float64, line)
    b[i] = v
    A[i,n+1] = b[i]
end

#println(A)
#=
for i in 1:n
    for j in 1:n
        print(A[i, j])
        print(" ")
    end
    println(" ")
end
=#
#println(b)
luDecomposition(A, b, n, l)
C = spzeros(Float64, n, n)
for i in 1:n
    for j in 1:n
        if (A[i, j] != 0)
            C[i,j] = A[i,j]
        end
    end
end
println("gauss:")
result = gaussEliminationOptimized(A, b, n, l)
for i in 1:n
    for j in 1:n
        print(A[i, j])
        print(" ")
    end
    println(" ")
end

#println(C)
println(result)
x = C\b
println(x)
