include("hilb.jl")

for i in 1:10
    A = hilb(i)
    x = ones(i,1)
    b = A * x

    x = A\b

    println("Gauss: ")
    println(x)
    println("Inv: ")

    x = inv(A)*b

    println(x)

end
