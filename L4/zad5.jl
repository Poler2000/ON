# Paweł Polerowicz 254626

include("./Interpolation.jl")

using .Interpolation

rysujNnfxDoPliku(x -> exp(x), 0.0, 1.0, 5, "5a_5.png")
rysujNnfxDoPliku(x -> exp(x), 0.0, 1.0, 10, "5a_10.png")
rysujNnfxDoPliku(x -> exp(x), 0.0, 1.0, 15, "5a_15.png")

rysujNnfxDoPliku(x -> x^2 * sin(x), -1.0, 1.0, 5, "5b_5.png")
rysujNnfxDoPliku(x -> x^2 * sin(x), -1.0, 1.0, 10, "5b_10.png")
rysujNnfxDoPliku(x -> x^2 * sin(x), -1.0, 1.0, 15, "5b_15.png")
