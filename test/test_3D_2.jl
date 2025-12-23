#using Test, OpenSWPC
#using GeophysicalModelGenerator


# test 2, define 3D velocity model from file and run it

# Define a 3D velocity model
N = 128,128,128
x = range(-5,7,length=N[1])
y = range(-5,7,length=N[2])
z = range(-9,3,length=N[3])
vmod = CartData(xyz_grid(x,y,z))

# Define seismic velocities
vp = fill(3.4,N...);
vs = fill(2.4,N...);
rho = fill(2.5,N...);
Qp = fill(200.0,N...);
Qs = fill(200.0,N...);

# Set a high density spherical anomaly
Temp = zero(rho)
add_sphere!(rho,Temp, vmod, cen=(0,0,-4), radius=2.5, phase=ConstantPhase(3.0))

# Define elastic properties from seismic velocities:
mu = rho.* vs .^2
lambda = rho.* (vp .^2 .- 2 .* vs .^2)

# Add fields to the CartData structure 
vmod = addfield(vmod,(;rho,Qp,Qs,mu,lambda, vp, vs))


s1 = SourceLLMWDC( 14.1338, 40.826221,  1.879,  0.1,    4,   3.9, 243.3627, 31.597,  -73.886 )

stations = [
    StationLL(14.1420, 40.8201, 0.170, "CAAM", "dep"),
    StationLL(14.1493, 40.8294, 0.110, "CSTH", "oba"),
]

# Specify model
cfg = OpenSWPCConfig(vmod,
                        input_file="input_velm.dat",
                        odir="cf_swp_layers",
                        tbeg=0.0,
                        dt = 0.0075, vcut=0.1, nt=300,
                        clon = 14.14, clat = 40.83,
                        xy_ps_sw=true,xz_ps_sw=true,yz_ps_sw=true,
                        vol_v_sw=true, vol_u_sw=true, vol_ps_sw=true,
                        z0_xy=2.5,stftype="triangle", ntdec_s=1,
                        source=[s1], stations=stations)

# Run model
write_input!(cfg)
run_swpc(cfg)

# Read some output back
dat,t=read_volume("cf_swp_layers/swpc.3d.vol.v.nc",timestep=100)
V = dat.fields.var"V_m/s"
@test all(extrema(V[2]) .≈ (-0.009220209f0, 0.008867903f0))

clean(cfg)

