using Plots

function nextX(x, c)
    return x^2 + c
end

function performIterations(c, x0)
    global currentX = x0
    global currentC = c
    myrange = 1:40
    plot(test.(myrange))

    savefig("myplot6.png")


    println("c = $c, x0 = $x0")
    x = x0;
    for i in 1:40
        x = Float64(nextX(x, c))
        println("x$i = $x")
    end
end

currentX = Float64(1)
currentC = Float64(-2)

function test(x)
    global currentX = nextX(currentX, currentC)
    return currentX
end

#performIterations(Float64(-2), Float64(1));
#performIterations(Float64(-2), Float64(2));
#performIterations(Float64(-2), Float64(1.99999999999999));
#performIterations(Float64(-1), Float64(1));
#performIterations(Float64(-1), Float64(-1));
#performIterations(Float64(-1), Float64(0.75));
performIterations(Float64(-1), Float64(0.25));
