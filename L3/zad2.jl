# Paweł Polerowicz 254626

function mstycznych(f,pf,x0::Float64, delta::Float64, epsilon::Float64, maxit::Int)
    v = f(x0)
    x1 = 0
    if abs(v) < epsilon
        return (x0,v,0,0)
    end
    for k in 1:maxit
        x1 = x0 - v / pf(x0)
        v = f(x1)
        if abs(pf(x0)) <= eps(Float64)
            return (x1,v,k,2)
        end
        if abs(x1 - x0) < delta || abs(v) < epsilon
            return (x1,v,k,0)
        end
        x0 = x1
    end
    return (x1, v, maxit, 1)
end
