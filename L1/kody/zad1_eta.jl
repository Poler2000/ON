# Paweł Polerowicz 254626
using Printf

# Float16
i = Float16(1.0)
eta = Float16(0.0)
while i > 0
    global eta = i
    global i /= 2
end
min16 = nextfloat(Float16(0))
@printf "eta16 %.2e nextfloat: %.2e\n" eta min16

# Float32
i = Float32(1.0)
eta = Float32(0.0)
while i > 0
    global eta = i
    global i /= 2
end
min32 = nextfloat(Float32(0.0))
@printf "eta32 %.2e nextfloat: %.2e\n" eta min32

# Float64
i = Float64(1.0)
eta = Float64(0.0)
while i > 0
    global eta = i
    global i /= 2
end

min64 = nextfloat(Float64(0))
@printf "eta64 %.2e nextfloat: %.2e\n" eta min64
