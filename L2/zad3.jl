# Paweł Polerowicz 254626

using LinearAlgebra
include("hilb.jl")
include("matcond.jl")

#Gauss
println("Hilbert's matrix: Gauss elimination & inv algorithm")
for i in 1:20
    A = hilb(i)
    x = ones(i,1)
    b = A * x
    onesnorm = norm(x)

    x = A\b
    gaussnorm = norm(x)
    gausserror = abs(onesnorm - gaussnorm) / onesnorm

    x = inv(A)*b
    invnorm = norm(x)
    inverror = abs(onesnorm - invnorm) / onesnorm

    println("$i & $gausserror & $inverror \\\\ \\hline")
end

nrange = [5,10,20]
crange = [1,10,10^3, 10^7, 10^12, 10^16]
println("Random matrix Gauss elimination & inv algorithm")
for i in nrange
    for j in crange
        A = matcond(i, Float64(j))
        x = ones(i,1)
        b = A * x
        onesnorm = norm(x)

        x = A\b
        gaussnorm = norm(x)
        gausserror = abs(onesnorm - gaussnorm) / onesnorm

        x = inv(A)*b
        invnorm = norm(x)
        inverror = abs(onesnorm - invnorm) / onesnorm

        println("$i & $j & $gausserror & $inverror \\\\ \\hline")
    end
end
