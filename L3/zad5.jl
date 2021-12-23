include("./Rootfinding.jl")

using .Rootfinding

myr, myv, myit, myerr = mbisekcji(x->exp(x) - 3x, 1.0, 2.0, 0.0001, 0.0001)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = mbisekcji(x->exp(x) - 3x, 0.0, 1.0, 0.0001, 0.0001)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = mbisekcji(x->exp(x) - 3x, 0.0, 2.0, 0.0001, 0.0001)

println("$myr, $myv, $myit, $myerr")
