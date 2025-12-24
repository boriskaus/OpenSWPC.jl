[![CI](https://github.com/boriskaus/OpenSWPC.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/boriskaus/OpenSWPC.jl/actions/workflows/ci.yml)

# OpenSWPC.jl 

`OpenSWPC.jl` is a julia interface to the [OpenSWPC code](https://openswpc.github.io), a MPI-parallel finite difference code to simulate seismic wave propagation developed by Takuto Maeda.

The Julia interface comes with precompiled binaries for linux and Mac, so you can directly run the code in parallel. It also provides a seamless integration with the [GeophysicalModelGenerator.jl](https://github.com/JuliaGeodynamics/GeophysicalModelGenerator.jl) julia package, such that you can directly run wave propagation simulations through any GMG setup. 


### 1. Installing and running

Installing and testing the code is easy: Install julia (version 1.10 or newer) and write:
```julia
julia>]
pkg> add https://github.com/boriskaus/OpenSWPC.jl
pkg> instantiate
pkg> test OpenSWPC
```
Note that both `OpenSWPC.jl` and `OpenSWPC_jll` (the precompiled binaries) are not yet registered, which is why the `instantiate` is necessary.

Also note that it only works on linux and mac. Windows users can install WSL and install julia in there to use the linux version of this code.


### 2. Example of using it

First, we create a 3D model setup using the GeophysicalModelGenerator 
```julia
# load the package
julia> using GeophysicalModelGenerator

# Create a simple 3D CartData model 
julia> nx,ny,nz=128,128,128
julia> x = range(-5,7,length=nx)
julia> y = range(-5,7,length=ny)
julia> z = range(-9,3,length=nz)
julia> data=  CartData(xyz_grid(x,y,z))

# Add a heterogeneity
julia> Temp = zero(rho)
julia> add_sphere!(rho,Temp, vmod, cen=(0,0,-4), radius=2.5, phase=ConstantPhase(3.0))

# Compute elastic properties:
julia> mu = rho.* vs .^2
julia> lambda = rho.* (vp .^2 .- 2 .* vs .^2)

# Add all fields to the model (using these names!)
julia> vmod = addfield(vmod,(;rho,Qp,Qs,mu,lambda, vp, vs))
```
Next, define an earthquake with focal mechanism:
```julia
julia> s1 = SourceLLMWDC( 14.1338, 40.826221,  1.879,  0.1,    4,   3.9, 243.3627, 31.597,  -73.886 )
```
and some seismic stations:
```julia
julia> stations = [
    StationLL(14.1420, 40.8201, 0.170, "CAAM", "dep"),
    StationLL(14.1493, 40.8294, 0.110, "CSTH", "oba"),
]
```

We can create an OpenSWPC model setup from that with:
```julia
# Specify model
julia> cfg = OpenSWPCConfig(vmod,
                        input_file="input_velm.dat",
                        odir="cf_swp_layers",
                        tbeg=0.0,
                        dt = 0.0075, vcut=0.1, nt=300,
                        clon = 14.14, clat = 40.83,
                        xy_ps_sw=true,xz_ps_sw=true,yz_ps_sw=true,
                        vol_v_sw=true, vol_u_sw=true, vol_ps_sw=true,
                        z0_xy=2.5,stftype="triangle", ntdec_s=1,
                        source=[s1], stations=stations)
OpenSWPCConfig:
  Control             : input_file="input_velm.dat", title="swpc", odir="cf_swp_layers", ntdec_r=50, strict_mode=false
  Grid Size           : nx=128, ny=128, nz=128
  Model Domain        : dx=0.09448818897637778, dy=0.09448818897637778, dz=0.09448818897637778, vcut=0.1, xbeg=-5.0, ybeg=-5.0, zbeg=-9.0
  Model Coord         : clon=14.14, clat=40.83, phi=0.0
  Frequency           : fq_min=0.02, fq_max=2.0, fq_ref=1.0
  Timestepping        : nt=300, dt=0.0075, tbeg=0.0
  Parallelisation     : nproc_x=2, nproc_y=2
  Output common       : ntdec_s=1, idec=2, jdec=2, kdec=2
  Output Free Surface : fs_v_sw=false, fs_u_sw=false, fs_ps_sw=false
  Output Ocean Bottom : ob_v_sw=true, ob_u_sw=true, ob_ps_sw=true
  Output XY Slice     : z0_xy=2.5, xy_v_sw=false, xy_u_sw=false, xy_ps_sw=true
  Output XZ Slice     : y0_xz=0.0, xz_ps_sw=true, xz_v_sw=true, yz_v_sw=false
  Output YZ Slice     : x0_yz=0.0, yz_ps_sw=true, yz_v_sw=false, yz_u_sw=false
  Output 3D           : vol_v_sw=true, vol_u_sw=true, vol_ps_sw=true
  Body Force Mode     : bf_mode=false
  Plane Wave Mode     : pw_mode=false, pw_ztop=100.0, pw_zlen=30.0, pw_ps="p", pw_strike=0.0, pw_dip=0.0, pw_rake=0.0
  Absorbing BCs       : abc_type="pml", na=20, stabilize_pml=false
  GMT Grid            : dir_grd="/vmodel/ejivsm/", fn_grdlst="grd.lst", node_grd=0
  Random Medium       : dir_rmed="./in/", fn_grdlst_rmed="grd.lst", rhomin=1.0, fn_rmed0="dummy.ns"
  Checkpoint/Restart  : is_ckp=false, ckpdir="./out/ckp", ckp_interval=1000000, ckp_time=1.0e6, ckp_seq=true
  Reciproc. Green's F.: green_mode=false, green_stnm="st01", green_cmp="z", green_trise=1.0, green_bforce=false, green_maxdist=550.0, green_fmt="llz", fn_glst="example/green.lst"
  MISC                : stopwatch_mode=false, benchmark_mode=false, ipad=0, jpad=0, kpad=0
  Waveform Output     : fn_stloc="stations.ll", wav_format="sac",ntdec_w_prg=0, sw_wav_v/u=(true, false), sw_wav_stress/strain=(false, false)
                         StationLL(stnm=CAAM, lon=14.142, lat=40.8201, dep=0.17 km, zsw=dep)
                         StationLL(stnm=CSTH, lon=14.1493, lat=40.8294, dep=0.11 km, zsw=oba)
  Earthquake Sources  : stf_format="llmwdc", stftype="triangle", fn_stf="sourceCF.dat", sdep_fit="asis", number of sources=1
                         SourceLLMWDC(lon=14.1338, lat=40.8262, z=1.879 km, tbeg=0.1 s, trise=4 s, Mw=3.9, DC[strike=243.4°, dip=31.6°, rake=-73.9°])
  Velocity Model      : vmodel_type="velm", is_ocean=false, topo_flatten=false, munk_profile=false, earth_flattening=false
                        Velocity model from file, dir_velm="./", fn_velm="velocity_model.nc" 
                        FullVelocityModel(3D CartDat velocity model, size=(128, 128, 128), fields=(:Z, :rho, :Qp, :Qs, :mu, :lambda, :vp, :vs))

```

And run this with:
```julia
julia> run_swpc(cfg)
```

Once the simulation is ready, you can create VTK (paraview) files from the snapshots with:
```julia
julia> movie_slice(cfg)
...
Saved PVD file
Movie saved to swpc_3d_vol_ps_nc.pvd
Movie-file (*.pvd) created for: cf_swp_layers/swpc_3d_xy_ps_nc.pvd
Movie-file (*.pvd) created for: cf_swp_layers/swpc_3d_xz_ps_nc.pvd
Movie-file (*.pvd) created for: cf_swp_layers/swpc_3d_yz_ps_nc.pvd
Movie-file (*.pvd) created for: cf_swp_layers/swpc_3d_ob_ps_nc.pvd
Movie-file (*.pvd) created for: cf_swp_layers/swpc_3d_xz_v_nc.pvd
Movie-file (*.pvd) created for: cf_swp_layers/swpc_3d_xz_u_nc.pvd
Movie-file (*.pvd) created for: cf_swp_layers/swpc_3d_vol_v_nc.pvd
Movie-file (*.pvd) created for: cf_swp_layers/swpc_3d_vol_u_nc.pvd
Movie-file (*.pvd) created for: cf_swp_layers/swpc_3d_vol_ps_nc.pvd
```
Open any of the `*.pvd` files with Paraview and play the movie.

### 3. Differences to OpenSWPC
OpenSWPC.jl ships with a precompiled binary version of the code. We have modified the fortran source code such that we can use 3D velocity models created with GMG in the models. Likewise, you can also

### 4. Getting help
All options of the `OpenSWPC` input file are configurable from the julia REPL. See the online documentation of [OpenSWPC](https://openswpc.github.io) tp understand the meaning of it all.

### 5. Citing
This package ofcourse relies on the great [OpenSWPC](https://github.com/OpenSWPC/OpenSWPC) code by Takuto Maeda. Please cite that if you use this in your work. 


### 6. Funding
Funding for this Julia interface was provided by the SAKURA project and by the DEGREE project funded by the German ministry of Science, Education and Space.