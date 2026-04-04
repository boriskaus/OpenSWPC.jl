#using Test, OpenSWPC
#using GeophysicalModelGenerator

# test 3; this uses the LaPalma volcano example of thhe GeophysicalModelGenerator.jl package as input 
# to define a 3D velocity model with a free surface
#

# Read topography
# This can in general be dome by using GMG and GMT:
#
# using GeophysicalModelGenerator, GMT
# Topo = import_topo(lon = [-18.2, -17.5], lat=[28.4, 29.0], file="@earth_relief_15s")
#
# yet, since GMT frequently fails, we have downloaded the topography and saved it to disk
# with save_GMG("topo_LaPalma", Topo)
# load this with:
Topo = load_GMG("topo_LaPalma")

# project to cartesian coordinates
proj = ProjectionPoint(Lon=-17.84, Lat=28.56)
Topo_model = CartData(xyz_grid(-35:.1:30,-15:.2:45,0));
Topo_model = project_CartData(Topo_model, Topo, proj)
write_paraview(Topo_model, "topo_LaPalma.vts")

# Define a 3D velocity model
#N = 128,128,128
N = 256,256,256

x = range(-20,20,length=N[1])
y = range(-15,35,length=N[2])
z = range(-20,5,length=N[3])
vmod = CartData(xyz_grid(x,y,z))

# Define seismic velocities
vp = fill(3.4,N...);
vs = fill(2.4,N...);
rho = fill(2.5,N...);
Qp = fill(200.0,N...);
Qs = fill(200.0,N...);

# Set a high density spherical anomaly
Temp = zero(rho)
add_sphere!(rho,Temp, vmod, cen=(0,20,-7), radius=5.0, phase=ConstantPhase(3.0))

# Set air properties for cells above the free surface (topography).
# OpenSWPC automatically detects the free surface from where vp=vs=0;
# without this step the model would use a flat free surface at z=zbeg.
air = above_surface(vmod, Topo_model)
vp[air]  .= 0.0
vs[air]  .= 0.0
rho[air] .= 0.001   # non-zero to avoid division by zero in OpenSWPC
Qp[air]  .= 10.0    # strong attenuation: waves do not propagate through air
Qs[air]  .= 10.0

# Define elastic properties from seismic velocities:
mu = rho.* vs .^2
lambda = rho.* (vp .^2 .- 2 .* vs .^2)

# Add fields to the CartData structure
vmod = addfield(vmod,(;rho,Qp,Qs,mu,lambda, vp, vs))
#write_paraview(vmod, "vmod_LaPalma")        # save to disk



#s1 = SourceLLMWDC( -17.84, 28.56,  10.0,  0.1,    2,   2.9, 243.3627, 31.597,  -73.886 )
s1 = SourceXYMWDC( 0, 0.1,  10.0,  0.1,    2,   2.9, 243.3627, 31.597,  -73.886 )

#"dep" sets the station at the 3D location
stations = [
    StationXY( 0.8,  -0.3,  0.170, "STA1", "dep"),
    StationXY( 1.5,   14.0, 1.210, "STA2", "dep"),
]

# Specify model
cfg = OpenSWPCConfig(vmod,
                        input_file="input_velm.dat",
                        odir="cf_swp_lapalma",
                        tbeg=0.0,
                        nproc_x=2, nproc_y=2,
                        dt = 0.01, vcut=0.1, nt=1000,
                        clon = proj.Lon, clat = proj.Lat,
                        xy_ps_sw=true,xz_ps_sw=true,yz_ps_sw=true,
                        vol_v_sw=false, vol_u_sw=false, vol_ps_sw=false,
                        z0_xy=2.5,stftype="triangle", 
                        
                        # set the horizontal reference depth above the free surface (topography); the actual air parameters are already set in 3d
                        topo0 = -maximum(Topo_model.z.val), 

                        # save every grid point, and every 20 timesteps:
                        ntdec_s=20,idec=1, jdec=1, kdec=1,  

                        # waveform output
                        sw_wav_u=true,sw_wav_v=true, ntdec_w=10,

                        source=[s1], stations=stations)

# Run model
write_input!(cfg)
run_swpc(cfg)

# Read some output back
ntsteps_saved = cfg.nt ÷ cfg.ntdec_s

dat,t=read_xz_slice(joinpath(cfg.odir, "swpc.3d.xz.v.nc"), timestep=ntsteps_saved)
V = dat.fields.var"V_m/s"
@test  sum(V[1]) ≈ -0.00026084436f0 rtol=1e-4

# Convert all NetCDF output files to ParaView format (*.pvd + *.vts)
movie_slice(cfg)


clean(cfg)
rm("topo_LaPalma.vts")  
rm("vmod_LaPalma.vts")  

