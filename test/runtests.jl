using Test
using GeophysicalModelGenerator
using OpenSWPC

#using .OpenSWPC
@testset "sourceCF read/write" begin
    # 1) Read sample llmwdc file
    sample = joinpath(@__DIR__, "input_tests", "sourceCF.dat")
    @show sample pwd()
    @test isfile(sample)
    sources = read_sourceCF(sample; format=:auto)
    @test length(sources) >= 1
    s = sources[1]
    @test s isa SourceLLMWDC
    @test isapprox(s.lon, 14.1338; atol=1e-4)
    @test isapprox(s.lat, 40.8262; atol=1e-6)
    @test isapprox(s.z, 1.879; atol=1e-6)
    @test isapprox(s.tbeg, 0.1; atol=1e-6)
    @test isapprox(s.trise, 4.0; atol=1e-6)
    @test isapprox(s.mag, 3.9; atol=1e-6)

    # 2) Round-trip llmwdc
    tmp = tempname()*".dat"
    write_sourceCF!(tmp, s)
    @test isfile(tmp)
    sources2 = read_sourceCF(tmp; format=:auto)
    @test length(sources2) == 1
    s2 = sources2[0+1]
    @test isapprox(s2.lon, s.lon; atol=1e-4)
    @test isapprox(s2.lat, s.lat; atol=1e-4)
    @test isapprox(s2.mag, s.mag; atol=1e-4)

    # 3) XY variant write + read
    sxy = SourceXYMWDC(0.0, 0.0, 1.879, 0.0, 2.5, 3.9, 243.3627, 31.597, -73.886)
    tmpxy = tempname()*".dat"
    write_sourceCF!(tmpxy, [sxy, sxy])
    sxys = read_sourceCF(tmpxy; format=:auto)
    @test length(sxys) == 2
    @test all(isa.(sxys, SourceXYMWDC))
    @test isapprox(sxys[1].z, 1.879; atol=1e-6)

    # 4) Header includes the chosen format
    txt = read(tmp, String)
    @test occursin("format 'llmwdc'", txt)
    txtxy = read(tmpxy, String)
    @test occursin("format 'xymwdc'", txtxy)
end

@testset "3D layered model" begin
    include(joinpath(@__DIR__, "test_3D_1.jl"))
end
@testset "3D mode with CartData input" begin
    include(joinpath(@__DIR__, "test_3D_2.jl"))
end
