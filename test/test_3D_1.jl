#using Test, OpenSWPC

# --------------------------------------------------------------------
# test 1: 3D layered velocity model with source and stations
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
cfg = OpenSWPCConfig(   odir="cf_swp_layers_1",
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

# Read some output back
name = joinpath(@__DIR__, "cf_swp_layers_1", "swpc.3d.xz.v.nc")
dat,t= read_xz_slice(name,timestep=100)
V    = dat.fields.var"V_m/s"


@test all(extrema(V[2]) .≈ (-0.014120221f0, 0.017427806f0))

clean(cfg)

# --------------------------------------------------------------------

