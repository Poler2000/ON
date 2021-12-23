include("./Rootfinding.jl")
using .Rootfinding

myr, myv, myit, myerr = mbisekcji(x->exp(1-x) - 1, 0.5, 2.0, 0.00001, 0.00001)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = mstycznych(x->exp(1-x) - 1, x->(-exp(1-x)), 0.5, 0.00001, 0.00001, 100)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = msiecznych(x->exp(1-x) - 1,  0.5, 2.0, 0.00001, 0.00001, 100)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = mbisekcji(x->x*exp(-x), -0.5, 2.0, 0.00001, 0.00001)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = mstycznych(x-> x * exp(-x), x->(x - 1) * exp(-x), 2.0, 0.00001, 0.00001, 100)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = msiecznych(x->x*exp(-x),  -0.5, 2.0, 0.00001, 0.00001, 100)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = mstycznych(x->exp(1-x) - 1, x->(-exp(1-x)), 2.0, 0.00001, 0.00001, 100)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = mstycznych(x->exp(1-x) - 1, x->(-exp(1-x)), 10.0, 0.00001, 0.00001, 100)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = mstycznych(x-> x * exp(-x), x->(x - 1) * exp(-x), 0.5, 0.00001, 0.00001, 100)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = mstycznych(x-> x * exp(-x), x->(x - 1) * exp(-x), 0.5, 0.00001, 0.00001, 100)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = mstycznych(x-> x * exp(-x), x->(x - 1) * exp(-x), 1.5, 0.00001, 0.00001, 100)

println("$myr, $myv, $myit, $myerr")

myr, myv, myit, myerr = mstycznych(x-> x * exp(-x), x->(x - 1) * exp(-x), 1.0, 0.00001, 0.00001, 100)

println("$myr, $myv, $myit, $myerr")
