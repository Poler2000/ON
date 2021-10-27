# Paweł Polerowicz 254626

kahanEps16 = Float16(3.0*(Float16(4.0/3.0) - 1.0) - 1.0)
kahanEps32 = Float32(3.0*(Float32(4.0/3.0) - 1.0) - 1.0)
kahanEps64 = Float64(3.0*(Float64(4.0/3.0) - 1.0) - 1.0)

println(kahanEps16)
println(kahanEps32)
println(kahanEps64)
