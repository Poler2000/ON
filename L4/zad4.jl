# Paweł Polerowicz 254626

include("./zad1.jl")
include("./zad2.jl")

using Plots

function rysujNnfx(f, a::Float64, b::Float64, n::Int)
    h = (b - a) / n
    # arguments
    argrange = Vector{Float64}(undef, n + 1)
    # values
    valrange = Vector{Float64}(undef, n + 1)
    # for performance we are using fixed number of points on the plot
    points = Vector{Float64}(undef, 51)

    for i in 0:n
        argrange[i+1] = a + h * i
        valrange[i+1] = f(a + h * i)
    end

    # we calculate differential quotients
    fx = ilorazyRoznicowe(argrange, valrange)
    step = (b - a) / 50
    nrange = a:step:b
    for i in 1:51
        points[i] = warNewton(argrange, fx, a + i * step)
    end

    plot(a:step:b, points, label = "interpolacja")
    plot!(a:step:b, f.(nrange), label = "f(x)")
    savefig("myplot.png")
end

# Helper for testing
function rysujNnfxDoPliku(f, a::Float64, b::Float64, n::Int, filename)
    h = (b - a) / n
    argrange = Vector{Float64}(undef, n + 1)
    valrange = Vector{Float64}(undef, n + 1)
    points = Vector{Float64}(undef, 51)

    for i in 0:n
        argrange[i+1] = a + h * i
        valrange[i+1] = f(a + h * i)
    end
    println(argrange)
    println(valrange)

    fx = ilorazyRoznicowe(argrange, valrange)
    step = (b - a) / 50
    nrange = a:step:b
    for i in 1:51
        points[i] = warNewton(argrange, fx, a + i * step)
    end

    plot(a:step:b, points, label = "interpolacja")
    plot!(a:step:b, f.(nrange), label = "f(x)")
    savefig(filename)
end
