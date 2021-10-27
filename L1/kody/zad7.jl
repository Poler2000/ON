# Paweł Polerowicz 254626

x = Float64(1.0)
result = Float64(0.0)

# aprrox of derivative
function f(x)
    return Float64(sin(x) + cos(3.0x))
end

# actual derivative
function fDerivative(x)
    return Float64(cos(x) - 3.0(sin(3.0x)))
end

derivative = fDerivative(x)
for n in 0:54
    h = Float64(2.0^-n)
    result = (f(x + h) - f(x)) / h
    myerror = Float64(abs(result - derivative))
    println("n: $n, result: $result, error: $myerror")
end
