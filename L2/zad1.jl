# Paweł Polerowicz 254626

x =  [2.718281828, -3.141592654, 1.414213562, 0.577215664, 0.301029995]
y =  [1486.2497, 878366.9879, -22.37492, 4773714.647, 0.000185049]
S32 = Float32(0.0)
S64 = Float64(0.0)

# forward
for i in 1:5
    global S32 += Float32(x[i]) * Float32(y[i])
    global S64 += x[i] * y[i]
end
println("forward")

println(S32)
println(S64)
S32 = Float32(0.0)
S64 = Float64(0.0)

# backwards
for i in 5:-1:1
    global S32 += Float32(x[i]) * Float32(y[i])
    global S64 += x[i] * y[i]
end
println("backwards")
println(S32)
println(S64)


# the biggest to the smallestDokładność
S32 = Float32(0.0)
S64 = Float64(0.0)
partialSums32 = Array{Float32}(undef, 5)
partialSums64 = Array{Float64}(undef, 5)
S32N = Float32(0.0)
S64N = Float64(0.0)
S32P = Float32(0.0)
S64P = Float64(0.0)

for i in 1:5
    global partialSums32[i] = Float32(x[i]) * Float32(y[i])
    global partialSums64[i] = x[i] * y[i]
end
sort!(partialSums32)
sort!(partialSums64)

i = 5
while i > 0 && partialSums32[i] >= 0.0
    global S32P += partialSums32[i]
    global i -= 1
end

i = 5
while i > 0 && partialSums64[i] >= 0.0
    global S64P += partialSums64[i]
    global i -= 1
end

i = 1
while i < 6 && partialSums32[i] < 0.0
    global S32N += partialSums32[i]
    global i += 1
end

i = 1
while i < 6 && partialSums64[i] < 0.0
    global S64N += partialSums64[i]
    global i += 1
end

S32 = S32N + S32P
S64 = S64N + S64P

println("big to small")
println(S32)
println(S64)


# the smallest to the biggest

S32 = Float32(0.0)
S64 = Float64(0.0)
S32N = Float32(0.0)
S64N = Float64(0.0)
S32P = Float32(0.0)
S64P = Float64(0.0)

i = 5
while i > 0
    if partialSums32[i] < 0
        global S32N += partialSums32[i]
    end
    global i -= 1
end

i = 5
while i > 0
    if partialSums64[i] < 0
        global S64N += partialSums64[i]
    end
    global i -= 1
end

i = 1
while i < 6
    if partialSums32[i] >= 0
        global S32P += partialSums32[i]
    end
    global i += 1
end

i = 1
while i < 6
    if partialSums64[i] >= 0
        global S64P += partialSums64[i]
    end
    global i += 1
end

S32 = S32N + S32P
S64 = S64N + S64P

println("small to big")
println(S32)
println(S64)
