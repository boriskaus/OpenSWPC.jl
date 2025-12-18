module OpenSWPC

using OpenSWPC_jll 
import MPI


if isdefined(OpenSWPC_jll,:MPICH_jll)
    const mpiexec = OpenSWPC_jll.MPICH_jll.mpiexec()
    const MPI_LIBPATH = OpenSWPC_jll.MPICH_jll.LIBPATH
elseif isdefined(OpenSWPC_jll,:MicrosoftMPI_jll) 
    const mpiexec = OpenSWPC_jll.MicrosoftMPI_jll.mpiexec()
    const MPI_LIBPATH = OpenSWPC_jll.MicrosoftMPI_jll.LIBPATH
elseif isdefined(OpenSWPC_jll,:OpenMPI_jll) 
    const mpiexec = OpenSWPC_jll.OpenMPI_jll.mpiexec()
    const MPI_LIBPATH = OpenSWPC_jll.OpenMPI_jll.LIBPATH
elseif isdefined(OpenSWPC_jll,:MPItrampoline_jll) 
    const mpiexec = OpenSWPC_jll.MPItrampoline_jll.mpiexec()
    const MPI_LIBPATH = OpenSWPC_jll.MPItrampoline_jll.LIBPATH
else
    println("Be careful! No MPI library detected; parallel runs won't work")
    const mpiexec = nothing
    const MPI_LIBPATH = Ref{String}("")
end

if Sys.iswindows()
    pathsep = ';'
elseif Sys.isapple()
    pathsep = ':'
else
    pathsep = ':'
end

include("utils.jl")

# Source file writer (sourceCF.dat)
include("source_cf.jl")
export AbstractSource, SourceLLMWDC, SourceXYMWDC, write_sourceCF!, read_sourceCF

# Station locations (station.ll)
include("stations.jl")
export StationLL, write_stations_ll!

# Velocity models --------------------------------
abstract type AbstractVelocityModel end

# Layered homogeneous medium (lhm)
include("lhm.jl")
export Layer1D, LHMModel, write_lhm!, read_lhm

# Uniform velocity model (vmodel = "uni")
include("vmodel_uni.jl")
export UniformVelocityModel

# ------------------------------------------------

# Random medium field
include("random_medium.jl")
export generate_random_medium

# Setup configuration of a model run
include("openswpc_input.jl")
export OpenSWPCConfig, write_input!

# run the code 
include("run.jl")   
export swpc_path, swpc_cmd, run_swpc

# read output files
include("read_output.jl")
export read_xy_slice, read_yz_slice, read_xz_slice, movie_slice

# 3D velocity model (vmodel = "user")
include("vmodel_3D.jl")
export VelocityModel3D


end # module OpenSWPC
