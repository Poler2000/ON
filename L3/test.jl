include("./Rootfinding.jl")

using .Rootfinding

println("Bisekcja - f(x) = 2x, x0 = -1, x1 = 2")
myr, myv, myit, myerr = mbisekcji(x->2x, -1.0, 2.0, 0.0001, 0.0001)

println("$myr, $myv, $myit, $myerr")

println("Bisekcja - f(x) = 2x, x0 = 1, x1 = 2")
myr, myv, myit, myerr = mbisekcji(x->2x, 1.0, 2.0, 0.0001, 0.0001)

println("$myr, $myv, $myit, $myerr")

println("Newtona - f(x) = 2x^2 - 1, x0 = 2")
myr, myv, myit, myerr = mstycznych(x->2(x^2) - 1, x-> 4x, 2.0, 0.000005, 0.000005, 100)

println("$myr, $myv, $myit, $myerr")

println("Siecznych - f(x) = 2x^2 - 1, x0 = 0, x1 = 2")
myr, myv, myit, myerr = msiecznych(x->2(x^2) - 1, 0.0, 2.0, 0.000005, 0.000005, 100)

println("$myr, $myv, $myit, $myerr")
