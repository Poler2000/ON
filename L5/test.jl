# Paweł Polerowicz 254626

include("./blocksys.jl")

using .blocksys


function testGauss(pathA, pathb, outputFile, pivoting::Bool)
    A, n, l = loadA(pathA)
    b = loadb(pathb)
    x = []
    if pivoting
        println("Gauss partial pivoting")
        @time for j in 1:1
            x = gaussPartialPivoting(A, b, n, l)
        end
        writeX(outputFile, x)
    else
        println("Gauss no pivoting")
        @time for j in 1:1
            x = gaussNoPivoting(A, b, n, l)
        end
        writeX(outputFile, x)
    end
    errors = Float64(0.0)

    for i in eachindex(x)
        errors += abs(1.0 - x[i])
    end
    println(errors)
    errors /= n
    println(errors)
    println("")

end

function testGauss(pathA, outputFile, pivoting::Bool)
    A, n, l = loadA(pathA)
    b = countb(A, n)
    x = []
    if pivoting
        println("Gauss partial pivoting")
        @time for j in 1:1
            x = gaussPartialPivoting(A, b, n, l)
        end
        writeX(outputFile, x)
    else
        println("Gauss no pivoting")
        @time for j in 1:1
            luPartialPivoting(A, b, n, l)
            x = luSolve(A, b, n, l)
        end
        writeXWithError(outputFile, x)
    end

end

function testLU(pathA, pathb, outputFile, pivoting::Bool)
    A, n, l = loadA(pathA)
    b = loadb(pathb)
    x = []
    if pivoting
        println("LU partial pivoting")
        @time luPartialPivoting(A, b, n, l)
    else
        println("LU no pivoting")
        @time luNoPivoting(A, b, n, l)
    end
    println("lu solve")
    @time for j in 1:1
        x = luSolve(A, b, n, l)
    end
    errors = Float64(0.0)

    for i in eachindex(x)
        errors += abs(1.0 - x[i])
    end
    println(errors)
    errors /= n
    println(errors)
    println("")

    writeX(outputFile, x)
end

function testLU(pathA, outputFile, pivoting::Bool)
    A, n, l = loadA(pathA)
    b = countb(A, n)
    x = []
    if pivoting
        println("LU partial pivoting")
        @time luPartialPivoting(A, b, n, l)
    else
        println("LU no pivoting")
        @time luNoPivoting(A, b, n, l)
    end
    println("lu solve")
    @time for j in 1:1
        x = luSolve(A, b, n, l)
    end

    writeXWithError(outputFile, x)
end

testGauss("Dane10000_1_1/A.txt", "Dane10000_1_1/b.txt", "result1.txt", true)
testGauss("Dane10000_1_1/A.txt", "Dane10000_1_1/b.txt", "result2.txt", false)
testLU("Dane10000_1_1/A.txt", "Dane10000_1_1/b.txt", "result3.txt", true)
testLU("Dane10000_1_1/A.txt", "Dane10000_1_1/b.txt", "result4.txt", false)


testLU("Dane100000_1_1/A.txt", "Dane100000_1_1/b.txt", "result5.txt", true)

#=
for i in 5000:5000:100000
    path = "A$i.txt"
    A, n, l = loadA(path)
    path = "b$i.txt"
    c = countb(A, n)
    b = countb(A, n)
    A_copy = Dict{Int64, Float64}()
    for (key, value) in A
        A_copy[key] = value
    end
    x = []
    t1 = @elapsed for j in 1:1
        luPartialPivoting(A_copy, c, n, l)
        x = luSolve(A_copy, c, n, l)
    end
    writeXWithError("x$i.txt", x)
end
=#

#=
for i in 5000:5000:100000
    blockmat(i, 4 ,10.0, "A$i.txt")
    f = open("b$i.txt", "w")
    println(f, i)
    for j in 1:i
        println(f, "1.0")
    end
end
=#
