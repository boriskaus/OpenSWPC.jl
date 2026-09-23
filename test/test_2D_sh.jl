#using Test, OpenSWPC
#using GeophysicalModelGenerator

# test 2D SH: same 2D velocity model setup as the P-SV test, run with swpc_sh

N = 200,1,100
x = range(-5,5,length=N[1])
z = range(-5,0,length=N[3])
vmod = CartData(xyz_grid(x,0,z))

vp  = fill(3.4,N...);
vs  = fill(2.4,N...);
rho = fill(2.5,N...);
Qp  = fill(200.0,N...);
Qs  = fill(200.0,N...);

mu = rho.* vs .^2
lambda = rho.* (vp .^2 .- 2 .* vs .^2)
vmod = addfield(vmod,(;rho,Qp,Qs,mu,lambda, vp, vs))

s1 = SourceXYMWDC( 0.0, 0.0, 2.0, 0.1, 0.5, 3.0, 0.0, 90.0, 0.0 )

stations = [
    StationXY(-2.0, 0.0, 0.1, "ST01", "dep"),
    StationXY( 2.0, 0.0, 0.1, "ST02", "dep"),
]

cfg = OpenSWPCConfig(vmod, solver="sh",
                        input_file="input_sh.dat",
                        odir="cf_swp_sh",
                        nproc_x=2,
                        dt = 0.005, vcut=0.1, nt=200,
                        xz_v_sw=true, xz_u_sw=false,
                        stftype="triangle", ntdec_s=50, idec=1, kdec=1,
                        sw_wav_v=true, ntdec_w=10,
                        source=[s1], stations=stations)

write_input!(cfg)
@test !occursin("xz_ps%sw", read(cfg.input_file, String))

run_swpc(cfg)

dat,t = read_xz_slice(joinpath(cfg.odir, "swpc.sh.xz.v.nc"), timestep=cfg.nt ÷ cfg.ntdec_s)
Vy = dat.fields.var"Vy_m/s"
@test size(Vy) == (N[1], 1, N[3])
@test all(isfinite, Vy)
@test maximum(abs, Vy) > 0

clean(cfg)
rm(joinpath(cfg.dir_velm, cfg.fn_velm))
