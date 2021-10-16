# Paweł Polerowicz 254626

mydelta = Float64(2^-52)
x = Float64(1.0)

while Float64(x*(1.0/x)) == Float64(1.0) && x < 2
    global x += mydelta
end

println("found number: $x")

# check
result = Float64(x*(1.0/x))
println(result)
println(result == Float64(1.0))
