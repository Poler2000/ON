function nextP(p, r)
    return p + r * p * (1.0 - p)
end

p0 = Float32(0.01)
r = Float32(3.0)
p = p0;


println("Float32, no round down")
for i in 1:40
    global p = Float32(nextP(p, r))
    println("p$i = $p")
end

p = p0;

println("Float32, round down after 10 iterations")
for i in 1:40
    global p = Float32(nextP(p, r))
    println("p$i = $p")
    if (i == 10)
        p = round(p, digits=3, RoundDown)
    end
end

println("Float64, round down after 10 iterations")
p0 = Float64(0.01)
r = Float64(3.0)
p = p0;
for i in 1:40
    global p = Float64(nextP(p, r))
    println("p$i = $p")
end
