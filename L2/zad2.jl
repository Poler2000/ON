# Paweł Polerowicz 254626

using Plots

function test(x)
    return Float64(exp(x) * log(1.0 + exp(-x)))
end

x = -20:0.5:100
plot(-20:0.5:100, test.(x), label = "e^xln(1+e^-x)")

savefig("newplot.png")
