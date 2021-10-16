# Paweł Polerowicz 254626

# Float16
i = Float16(1.0)
macheps = Float16(1.0)

while Float16(1.0 + i) > 1.0 && Float16(1.0 + i) == (1.0 + i)
    global macheps = i
    global i /= 2
end

builtinEps = eps(Float16)
println("macheps16: $macheps eps: $builtinEps")

# Float32
i = Float32(1.0)
macheps = Float32(1.0)

while Float32(1.0 + i) > 1.0 && Float32(1.0 + i) == (1.0 + i)
    global macheps = i
    global i /= 2
end
builtinEps = eps(Float32)
println("macheps32: $macheps eps: $builtinEps")

# Float64
i = Float64(1.0)
macheps = Float64(1.0)

while Float64(1.0 + i) > 1.0 && Float64(1.0 + i) == (1.0 + i)
    global macheps = i
    global i /= 2
end
builtinEps = eps(Float64)
println("macheps64: $macheps eps: $builtinEps")
