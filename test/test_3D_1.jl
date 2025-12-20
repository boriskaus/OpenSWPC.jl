using Test, OpenSWPC
using GeophysicalModelGenerator

# --------------------------------------------------------------------
# test 1: 3D layered velocity model with source and stations
#=
# Define a layered velocity model
m = LHMModel([Layer1D(0;   rho=2.0, vp=1.2, vs=0.8, qp=50, qs=30),
                Layer1D(2;   rho=2.2, vp=2.0, vs=1.2, qp=80, qs=50),
                Layer1D(4;   rho=2.2, vp=3.0, vs=2.0, qp=80, qs=50),
                Layer1D(5;   rho=2.4, vp=4.0, vs=2.5, qp=80, qs=50),
               ])

s1 = SourceLLMWDC( 14.1338, 40.826221,  1.879,  0.1,    4,   3.9, 243.3627, 31.597,  -73.886 )

stations = [
    StationLL(14.1420, 40.8201, 0.170, "CAAM", "dep"),
    StationLL(14.1493, 40.8294, 0.110, "CSTH", "oba"),
]

# Specify model
cfg = OpenSWPCConfig(   odir="cf_swp_layers",
                        nx = 128,   ny=128,  nz=128,
                        dx = 0.095, dy=0.095,dz=0.095,   
                        xbeg=-5,ybeg=-5,zbeg=-3.0,tbeg=0.0,
                        dt = 0.01, vcut=0.1, nt=300,
                        clon = 14.14, clat = 40.83,
                        xy_ps_sw=true,xz_ps_sw=true,yz_ps_sw=true,
                        z0_xy=2.5,stftype="triangle", ntdec_s=1,
                        vmodel=m, source=[s1], stations=stations)

# Run model
run_swpc(cfg)
=#
# --------------------------------------------------------------------


# test 2, define 3D velocity model from file and run it

# Define a 3D velocity model
nx,ny,nz = 128, 128, 128
dx,dy,dz = 0.095, 0.095, 0.095   
x        = -5:dx:(nx-1)*dx-5
y        = -5:dy:(ny-1)*dy-5
z        = range(-((nz-1)*dz-3),3.0,length=nz)  # ensure proper spacing

uni = UniformVelocityModel(3.4,2.4,2.5)

vp = fill(uni.vp0,nx,ny,nz);
vs = fill(uni.vs0,nx,ny,nz);
rho = fill(uni.rho0,nx,ny,nz);
Qp = fill(uni.qp0,nx,ny,nz);
Qs = fill(uni.qs0,nx,ny,nz);

vmod = CartData(xyz_grid(x,y,z))

Temp = zero(rho)
add_sphere!(rho,Temp, vmod, cen=(0,0,-4), radius=2.5, phase=ConstantPhase(3.0))

mu = rho.* vs .^2
lambda = rho.* (vp .^2 .- 2 .* vs .^2)
vmod = addfield(vmod,(;rho,Qp,Qs,mu,lambda, vp, vs))

write_paraview(vmod, "vmodel")

s1 = SourceLLMWDC( 14.1338, 40.826221,  1.879,  0.1,    4,   3.9, 243.3627, 31.597,  -73.886 )

stations = [
    StationLL(14.1420, 40.8201, 0.170, "CAAM", "dep"),
    StationLL(14.1493, 40.8294, 0.110, "CSTH", "oba"),
]

# Specify model
cfg = OpenSWPCConfig(   input_file="input_velm.dat",
                        odir="cf_swp_layers",
                        nx = 128,   ny=128,  nz=128,
                        dx = 0.095, dy=0.095,dz=0.095,   
                        xbeg=-5,ybeg=-5,zbeg=-3.0,tbeg=0.0,
                        dt = 0.0075, vcut=0.1, nt=300,
                        clon = 14.14, clat = 40.83,
                        xy_ps_sw=true,xz_ps_sw=true,yz_ps_sw=true,
                        z0_xy=2.5,stftype="triangle", ntdec_s=1,
                        vmodel=vmod,
                        source=[s1], stations=stations)

cfg.vmodel_type = "velm"  # to be added to OpenSWPCConfig


# Run model
write_netcdf(vmod, "velocity_model.nc")
write_input!(cfg)
#run_swpc(cfg)
