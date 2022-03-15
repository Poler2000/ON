# Paweł Polerowicz 254626

module blocksys

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

# load A matrix from file
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

# load b vector from file
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

function countb(A::Dict{Int64, Float64}, n::Int64)
    b = zeros(Float64, n)
    for (key, value) in A
        x = ((key - 1) % n) + 1
        b[x] += value
    end

    return b
end

# write x vector to file
function writeX(outputFile, x)
    f = open(outputFile, "w")
    n = length(x)
    for i in 1:n
        println(f, x[i])
    end
end

# write x vector to file with error value
function writeXWithError(outputFile, x)
    f = open(outputFile, "w")
    n = length(x)

    errors = Float64(0.0)

    for i in eachindex(x)
        errors += abs(1.0 - x[i])
    end
    println(f, errors)
    for i in 1:n
        println(f, x[i])
    end
end

export gaussNoPivoting, gaussPartialPivoting, luNoPivoting, luPartialPivoting, luSolve
export loadA, loadb, writeX, writeXWithError, countb
end
