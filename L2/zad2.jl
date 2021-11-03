using Plots

function test(x)
    return exp(x)
end

x = 1:100
plot(test.(x))

savefig("myplot1.png")
