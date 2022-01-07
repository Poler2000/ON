# Paweł Polerowicz 254626

using SparseArrays
using LinearAlgebra

function luPartialPivoting(A::Dict{Int64, Float64}, b::Vector{Float64}, n::Int64, l::Int64)
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
        upper = min(n, i + 2 * l)
        lower = i - l
        if (i_max != i)
            tmp = b[i]
            b[i] = b[i_max]
            b[i_max] = tmp

            for j in lower:upper
                tmp = get(A, i + (j - 1) * n, 0.0)
                if (tmp != 0 || get(A, (i_max + (j - 1) * n), 0.0) != 0)
                    A[i + (j - 1) * n] = get(A, (i_max + (j - 1) * n), 0.0)
                    A[i_max + (j - 1) * n] = tmp
                end
            end

        end

        for j in i+1:upper
            A[j + (i - 1) * n] = get(A, (j + (i - 1) * n), 0.0) / get(A, (i + (i - 1) * n), 0.0)
            for k in i+1:upper
                A[j + (k - 1) * n] = get(A, (j + (k - 1) * n), 0.0) - get(A, (j + (i - 1) * n), 0.0) * get(A, (i + (k - 1) * n), 0.0)
            end
        end
    end
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
end

function luSolve(A::Dict{Int64, Float64}, b::Vector{Float64}, n::Int64, l::Int64)
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

        upper = min(n, i + 2*l + 1)
        for k in (i+1):upper
            sum = sum + get(A, i + (k - 1) * n, 0.0) * x[k]
        end
        x[i] = (1 / get(A, i + (i - 1) * n, 0.0)) * (y[i] - sum)
    end
    return x
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
    b = loadb(path)
    c = loadb(path)
    A_copy = Dict{Int64, Float64}()
    for (key, value) in A
        A_copy[key] = value
    end

    luNoPivoting(A, b, n, l)
    x = luSolve(A, b, n, l)

    errors = Float64(0.0)

    for i in eachindex(x)
        errors += abs(1.0 - x[i])
    end
    print("$errors & ")
    errors /= n
    print("$errors & ")

    luPartialPivoting(A_copy, c, n, l)
    x = luSolve(A_copy, c, n, l)

    errors = Float64(0.0)

    for i in eachindex(x)
        errors += abs(1.0 - x[i])
    end
    print("$errors & ")
    errors /= n
    print("$errors")

    println(" ")
end=#
#=

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

    #=t1 = @elapsed luNoPivoting(A, b, n, l)
    t3 = @elapsed luSolve(A, b, n, l)
    t1 += t3
    t2 = @elapsed luPartialPivoting(A_copy, c, n, l)
    t4 = @elapsed luSolve(A_copy, c, n, l)
    t2 += t4
    println("$t1; $t2")=#
end=#

for i in 5000:5000:100000
    path = "A$i.txt"
    A, n, l = loadA(path)
    path = "b$i.txt"
    b = loadb(path)
    A_copy = spzeros(Float64, n, n)
    for (key, value) in A
        if key != i * i
            A_copy[((key - 1) % n) + 1, floor(Int64, key / n) + 1] = value
        else
            A_copy[i, i] = value
        end
    end
    L,U,p = lu(A_copy)
    #U\(L\b[p])
    #=t1 = @elapsed for j in 1:1
        L,U,p = lu(A_copy)
        U\(L\b[p])
    end=#

    #t1 = @elapsed A_copy \ b
    #println("$t1")

    #=t1 = @elapsed luNoPivoting(A, b, n, l)
    t3 = @elapsed luSolve(A, b, n, l)
    t1 += t3
    t2 = @elapsed luPartialPivoting(A_copy, c, n, l)
    t4 = @elapsed luSolve(A_copy, c, n, l)
    t2 += t4
    println("$t1; $t2")=#
end
