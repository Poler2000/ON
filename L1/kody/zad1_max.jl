# Paweł Polerowicz

max16 = floatmax(Float16);
max32 = floatmax(Float32);
max64 = floatmax(Float64);

# helper variable
i = Float16(1)
# will hold max value at the end
mymax = Float16(1)

# we search for the biggest number that is a power of 2
while !isinf(i)
    global mymax = i
    global i *= 2
end

# number to add
addition = Float16(mymax);
# prev value of mymax
prev = Float16(0)

# we try to add smaller and smaller numbers to i
while addition > 0.0
    global prev = 0
    # add addition as long as i < inf and
    # addition change the number
    while !isinf(i) && i != prev
        global mymax = i
        global prev = i
        global i += addition
    end
    global i = mymax
    global addition /= 2
end

println("mymax: $mymax, floatmax: $max16")

# helper variable
i = Float32(1)
# will hold max value at the end
mymax = Float32(1)

# we search for the biggest number that is a power of 2
while !isinf(i)
    global mymax = i
    global i *= 2
end

# number to add
addition = Float32(mymax);
# prev value of mymax
prev = Float32(0)

# we try to add smaller and smaller numbers to i
while addition > 0.0
    global prev = 0
    # add addition as long as i < inf and
    # addition change the number
    while !isinf(i) && i != prev
        global mymax = i
        global prev = i
        global i += addition
    end
    global i = mymax
    global addition /= 2
end

println("mymax: $mymax, floatmax: $max32")

# helper variable
i = Float64(1)
# will hold max value at the end
mymax = Float64(1)

# we search for the biggest number that is a power of 2
while !isinf(i)
    global mymax = i
    global i *= 2
end

# number to add
addition = Float64(mymax);
# prev value of mymax
prev = Float64(0)

# we try to add smaller and smaller numbers to i
while addition > 0.0
    global prev = 0
    # add addition as long as i < inf and
    # addition change the number
    while !isinf(i) && i != prev
        global mymax = i
        global prev = i
        global i += addition
    end
    global i = mymax
    global addition /= 2
end

println("mymax: $mymax, floatmax: $max64")
