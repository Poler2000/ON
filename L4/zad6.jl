# Paweł Polerowicz 254626

include("./Interpolation.jl")

using .Interpolation

rysujNnfxDoPliku(x -> abs(x), -1.0, 1.0, 5, "6a_5.png")
rysujNnfxDoPliku(x -> abs(x), -1.0, 1.0, 10, "6a_10.png")
rysujNnfxDoPliku(x -> abs(x), -1.0, 1.0, 15, "6a_15.png")

rysujNnfxDoPliku(x -> 1 / (1 + x^2), -5.0, 5.0, 5, "6b_5.png")
rysujNnfxDoPliku(x -> 1 / (1 + x^2), -5.0, 5.0, 10, "6b_10.png")
rysujNnfxDoPliku(x -> 1 / (1 + x^2), -5.0, 5.0, 15, "6b_15.png")
