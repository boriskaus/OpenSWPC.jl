# 2D model from a cross-section through a 3D GMG model, run with swpc_psv

N = 60,60,40
x = range(-6,6,length=N[1])
y = range(-6,6,length=N[2])
z = range(-5,0,length=N[3])
vmod3D = CartData(xyz_grid(x,y,z))

vp  = fill(3.4,N...);
vs  = fill(2.4,N...);
rho = fill(2.5,N...);
Qp  = fill(200.0,N...);
Qs  = fill(200.0,N...);
Temp = zero(rho)
add_sphere!(rho,Temp, vmod3D, cen=(1,1,-3), radius=1.0, phase=ConstantPhase(3.0))
mu = rho.* vs .^2
lambda = rho.* (vp .^2 .- 2 .* vs .^2)
vmod3D = addfield(vmod3D,(;rho,Qp,Qs,mu,lambda, vp, vs))

# oblique profile: (n, nz, 1) -> (n, 1, nz) with x = distance along the profile
cs = cross_section(vmod3D, Start=(-4,-4), End=(4,4), dims=(160,80))
@test size(cs) == (160,80,1)
m2d = to_2D(cs)
@test size(m2d) == (160,1,80)
dxs = diff(m2d.x.val[:,1,1])
@test all(isapprox.(dxs, dxs[1]; rtol=1e-6)) && dxs[1] > 0
@test m2d.x.val[1,1,1] ≈ 0 && m2d.x.val[end,1,1] ≈ 8*sqrt(2)
@test maximum(m2d.fields.rho) ≈ 3.0        # the sphere is on the diagonal

# fixed-x section: (1, ny, nz) -> (ny, 1, nz)
cs_x = cross_section(vmod3D, Lon_level=1.0)
@test size(to_2D(cs_x)) == (N[2],1,N[3])

# a plain 3D model is rejected
@test_throws ErrorException to_2D(vmod3D)

# OpenSWPCConfig converts the section itself; source and stations are in profile coordinates
L = m2d.x.val[end,1,1]
s1 = SourceXYMWDC( L/2, 0.0, 2.0, 0.1, 0.5, 3.0, 0.0, 45.0, 90.0 )
stations = [StationXY(L/4, 0.0, 0.1, "ST01", "dep"), StationXY(3L/4, 0.0, 0.1, "ST02", "dep")]
cfg = OpenSWPCConfig(cs,
                        input_file="input_section.dat",
                        odir="cf_swp_section",
                        nproc_x=2,
                        dt = 0.005, vcut=0.1, nt=100,
                        xz_ps_sw=false, xz_v_sw=true, xz_u_sw=false,
                        stftype="triangle", ntdec_s=50, idec=1, kdec=1,
                        sw_wav_v=true, ntdec_w=10,
                        source=[s1], stations=stations)
@test cfg.solver == "psv"
@test (cfg.nx, cfg.nz) == (160, 80)
@test cfg.xbeg ≈ 0

write_input!(cfg)
run_swpc(cfg)

dat,t = read_xz_slice(joinpath(cfg.odir, "swpc.psv.xz.v.nc"), timestep=cfg.nt ÷ cfg.ntdec_s)
V = dat.fields.var"V_m/s"
@test size(V[1]) == (160, 1, 80)
@test all(isfinite, V[1]) && all(isfinite, V[3])
@test maximum(abs, V[3]) > 0

clean(cfg)
rm(joinpath(cfg.dir_velm, cfg.fn_velm))
