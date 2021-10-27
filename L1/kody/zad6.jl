# Paweł Polerowicz 254625

# argument
x = Float64(1/8)

#iteration
i = 1

while x > 0
    f = Float64(sqrt(x^2 + 1.0) - 1.0)
    g = Float64(x^2 / (sqrt(x^2 + 1.0) + 1.0))

    println("x = $x, i = $i, f(x) = $f, g(x) = $g")
    global x /= 8
    global i += 1
end
