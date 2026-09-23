using Test
using GeophysicalModelGenerator
using OpenSWPC
using NCDatasets

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

@testset "StationXY read/write" begin
    stations = [
        StationXY(10.5, 20.3, 0.5, "ST01", "dep"),
        StationXY(15.0, 25.0, 1.0, "ST02", "obb"),
    ]

    # write and check file exists
    tmp = tempname() * ".xy"
    write_stations_xy!(tmp, stations)
    @test isfile(tmp)

    # check format keyword in header
    txt = read(tmp, String)
    @test occursin("stloc.xy", txt)
    @test occursin("Cartesian", txt)

    # check all station values appear in file
    @test occursin("ST01", txt)
    @test occursin("ST02", txt)
    @test occursin("10.5000", txt)
    @test occursin("'dep'", txt)
    @test occursin("'obb'", txt)

    # single-station convenience form
    tmp2 = tempname() * ".xy"
    write_stations_xy!(tmp2, StationXY(0.0, 0.0, 0.0, "S1", "fsb"))
    @test isfile(tmp2)
    @test occursin("S1", read(tmp2, String))

    rm(tmp, force=true)
    rm(tmp2, force=true)
end

@testset "3D layered model" begin
    include(joinpath(@__DIR__, "test_3D_1.jl"))
end
@testset "3D model with CartData input" begin
    include(joinpath(@__DIR__, "test_3D_2.jl"))
end
@testset "3D model with topography input" begin
    include(joinpath(@__DIR__, "test_3D_3.jl"))
end
@testset "2D netCDF, random medium and show" begin
    include(joinpath(@__DIR__, "test_2D_io.jl"))
end
@testset "2D P-SV model with CartData input" begin
    include(joinpath(@__DIR__, "test_2D_psv.jl"))
end
@testset "2D SH model with CartData input" begin
    include(joinpath(@__DIR__, "test_2D_sh.jl"))
end
