#using Test, OpenSWPC, NCDatasets
#using GeophysicalModelGenerator

# test 2D P-SV: define a 2D velocity model (CartData of size (nx,1,nz)) and run it with swpc_psv

N = 200,1,100
x = range(-5,5,length=N[1])
z = range(-5,0,length=N[3])
vmod = CartData(xyz_grid(x,0,z))

vp  = fill(3.4,N...);
vs  = fill(2.4,N...);
rho = fill(2.5,N...);
Qp  = fill(200.0,N...);
Qs  = fill(200.0,N...);

Temp = zero(rho)
add_sphere!(rho,Temp, vmod, cen=(0,0,-3), radius=1.0, phase=ConstantPhase(3.0))

mu = rho.* vs .^2
lambda = rho.* (vp .^2 .- 2 .* vs .^2)
vmod = addfield(vmod,(;rho,Qp,Qs,mu,lambda, vp, vs))

s1 = SourceXYMWDC( 0.0, 0.0, 2.0, 0.1, 0.5, 3.0, 0.0, 45.0, 90.0 )

stations = [
    StationXY(-2.0, 0.0, 0.1, "ST01", "dep"),
    StationXY( 2.0, 0.0, 0.1, "ST02", "dep"),
]

cfg = OpenSWPCConfig(vmod,
                        input_file="input_psv.dat",
                        odir="cf_swp_psv",
                        nproc_x=2,
                        dt = 0.005, vcut=0.1, nt=200,
                        xz_ps_sw=false, xz_v_sw=true, xz_u_sw=false,
                        stftype="triangle", ntdec_s=50, idec=1, kdec=1,
                        sw_wav_v=true, ntdec_w=10,
                        source=[s1], stations=stations)

@test cfg.solver == "psv"
@test (cfg.nx, cfg.nz) == (N[1], N[3])

write_input!(cfg)
inp = read(cfg.input_file, String)
@test occursin("nx               = $(N[1])", inp)
@test !occursin(r"^\s*(ny|dy|ybeg|nproc_y)\s*="m, inp)
@test occursin("xz_v%sw", inp) && !occursin("xy_v%sw", inp)

vfile = joinpath(cfg.dir_velm, cfg.fn_velm)
NCDataset(vfile) do ds
    @test collect(keys(ds.dim)) == ["x","z"]
    @test ds["rho"][:,:] == Float32.(reverse(rho[:,1,:], dims=2))    # z flipped: k=1 is the top
end

run_swpc(cfg)

dat,t = read_xz_slice(joinpath(cfg.odir, "swpc.psv.xz.v.nc"), timestep=cfg.nt ÷ cfg.ntdec_s)
V = dat.fields.var"V_m/s"
@test size(V[1]) == (N[1], 1, N[3])
@test all(isfinite, V[1]) && all(isfinite, V[3])
@test maximum(abs, V[3]) > 0

clean(cfg)
rm(vfile)
