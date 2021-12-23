# Paweł Polerowicz 254626

function msiecznych(f, x0::Float64, x1::Float64, delta::Float64,
    epsilon::Float64, maxit::Int)
    a = x0
    b = x1
    fa = f(a)
    fb = f(b)

    for k in 1:maxit
        if abs(fa) > abs(fb)
            tmp = a
            a = b
            b = tmp
            tmp = fa
            fa = fb
            fb = tmp
        end
        s = (b - a) / (fb - fa)
        b = a
        fb = fa
        a = a - fa * s
        fa = f(a)
        if abs(b - a) < delta || abs(fa) < epsilon
            return (a,fa,k,0)
        end
    end
    return (a,fa,maxit,1)
end
