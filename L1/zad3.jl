# Paweł Polerowicz 254626

step = Float64(2^-52)
i = Float64(1.0)

while i <= 2.0
    println(bitstring(i))
    global i += step
end
