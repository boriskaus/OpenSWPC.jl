#using Test, OpenSWPC
#using GeophysicalModelGenerator

# test netCDF round-trip, random media and printing for 2D models

for N in ((20,1,10), (20,8,10))
    local x = range(-5,5,length=N[1])
    local z = range(-5,0,length=N[3])
    local data = CartData(xyz_grid(x, N[2] == 1 ? 0 : range(-2,2,length=N[2]), z))
    local rho = rand(N...)
    local vp  = rand(N...)
    data = addfield(data, (;rho, vp))

    local file = tempname()*".nc"
    write_netcdf(data, file)
    local data2 = read_netcdf(file)
    @test size(data2) == N
    @test data2.x.val == Float32.(data.x.val)
    @test data2.fields.rho == Float32.(rho)
    @test data2.fields.vp == Float32.(vp)
    rm(file)
end

N = 32,1,16
vmod = CartData(xyz_grid(range(0,3.1,length=N[1]), 0, range(-1.5,0,length=N[3])))
fields = (; rho=fill(2.5,N...), Qp=fill(200.0,N...), Qs=fill(200.0,N...), mu=fill(14.4,N...), lambda=fill(0.1,N...))
cfg = OpenSWPCConfig(addfield(vmod, fields))

rfile = tempname()*".nc"
generate_random_medium(cfg, outfile=rfile)
rmed = read_netcdf(rfile)
@test size(rmed) == N
@test all(isfinite, rmed.fields.var"random media")
rm(rfile)

txt = sprint(show, MIME"text/plain"(), cfg)
@test occursin("nx=$(N[1])", txt) && occursin("xz_ps_sw", txt)
@test !occursin(r"\bny=|dy=|ybeg=|nproc_y=|xy_v_sw|yz_v_sw|fs_v_sw|vol_v_sw", txt)
@test !occursin("xz_ps_sw", sprint(show, MIME"text/plain"(), OpenSWPCConfig(solver="sh")))
@test occursin("ny=", sprint(show, MIME"text/plain"(), OpenSWPCConfig()))
