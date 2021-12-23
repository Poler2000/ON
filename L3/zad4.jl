include("./Rootfinding.jl")

using .Rootfinding

myr, myv, myit, myerr = mbisekcji(x->sin(x) - (0.5x)^2, 1.5, 2.0, 0.000005, 0.000005)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = mstycznych(x->sin(x) - (0.5x)^2, x-> cos(x) - 0.5x, 1.5, 0.000005, 0.000005, 100)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = msiecznych(x->sin(x) - (0.5x)^2, 1.5, 2.0, 0.000005, 0.000005, 100)

println("$myr, $myv, $myit, $myerr")
