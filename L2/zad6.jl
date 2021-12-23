# Paweł Polerowicz 254626

using Plots

function nextX(x, c)
    return x^2 + c
end

currentX = Float64(0.75)
currentC = Float64(-1)

function performIterations(c, x0, filename)
    global currentX = x0
    global currentC = c
    funcrange = -3:0.1:3

    values = Array{Float64}(undef, 80)
    arguments = Array{Float64}(undef, 80)

    for i in 1:40
        arguments[2*i - 1] = currentX
        values[2*i - 1] = test(currentX)
        values[2*i] = currentX
        arguments[2*i] = currentX
    end

    x = nextX.(funcrange, c); y = -3:0.1:3
    plot(funcrange, x, label="y=x^2-1")
    plot!(funcrange, y, label="y=x")
    plot!(arguments, values, xlims=(-3,3), label="graphical iteration")


    savefig(filename)


    println("c = $c, x0 = $x0")
    x = x0;
    for i in 1:40
        x = Float64(nextX(x, c))
        #println("x$i = $x")
        println(x)
    end
end

function test(x)
    global currentX = nextX(currentX, currentC)
    return currentX
end

performIterations(Float64(-2), Float64(1), "plot1.png");
performIterations(Float64(-2), Float64(2), "plot2.png");
performIterations(Float64(-2), Float64(1.99999999999999), "plot3.png");
performIterations(Float64(-1), Float64(1), "plot4.png");
performIterations(Float64(-1), Float64(-1), "plot5.png");
performIterations(Float64(-1), Float64(0.75), "plot6.png");
performIterations(Float64(-1), Float64(0.25), "plot7.png");
