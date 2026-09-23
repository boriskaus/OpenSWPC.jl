# generates a random medium velocity model 
using OpenSWPC

export generate_random_medium

"""
    generate_random_medium(;nx=100,ny=100,nz=100, 
                            dx=0.1, dy=0.1, dz=0.1,
                            ax=1,ay=1,az=1,
                            epsil=0.1, kappa=1.1,
                            ptype=1,
                            outfile="test_random.nc")

Generates a random medium 3D velocity model using the OpenSWPC `gen_rmed3d` utility.                                    
- nx, ny, nz: number of grid points in x, y, z directions
- dx, dy, dz: grid spacing in x, y, z directions
- ax, ay, az: correlation lengths in x, y, z directions
- epsil: standard deviation of velocity fluctuations
- kappa: spectral index
- ptype: type of random medium (1: Gaussian, 2: Exponential, 3: von Karman
- outfile: name of the output NetCDF file

"""
function generate_random_medium(;   nx=64,ny=64,nz=64, 
                                    dx=0.1, dy=0.1, dz=0.1,
                                    ax=1,ay=1,az=1,
                                    epsil=0.1, kappa=1.1,
                                    ptype::Int=1,
                                    outfile="random_field.nc")


    
    exe = OpenSWPC.OpenSWPC_jll.gen_rmed3d()
    cmd = `$exe -o $outfile -nx $nx -ny $ny -nz $nz -dx $dx -dy $dy -dz $dz  -ax $ax -ay $ay -az $az -epsil $epsil -kappa $kappa -ptype $ptype`
    run(cmd);

    println("Random medium velocity model written to: $(outfile)")

    return nothing
end

"""
    generate_random_medium(cfg::OpenSWPCConfig; kwargs...)
Generates a random medium consistent with the model parameters in `cfg`.
For a 2D `cfg` the medium has `ny = 1`, so `read_netcdf` returns a `CartData` of size `(nx, 1, nz)`.
"""
function generate_random_medium(cfg::OpenSWPCConfig; kwargs...)
    ny, dy = cfg.solver == "3d" ? (cfg.ny, cfg.dy) : (1, cfg.dx)    # gen_rmed3d returns NaN for dy=0
    generate_random_medium(  nx=cfg.nx, ny=ny, nz=cfg.nz,
                                dx=cfg.dx, dy=dy, dz=cfg.dz; kwargs...)
end



