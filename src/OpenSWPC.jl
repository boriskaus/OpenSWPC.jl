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

# Source file writer (sourceCF.dat)
include("source_cf.jl")
using .OpenSWPCSource: AbstractSource, SourceLLMWDC, SourceXYMWDC, write_sourceCF!, read_sourceCF
export AbstractSource, SourceLLMWDC, SourceXYMWDC, write_sourceCF!, read_sourceCF

# Station locations (station.ll)
include("stations.jl")
using .OpenSWPCStations: StationLL, write_stations_ll!
export StationLL, write_stations_ll!

# Layered homogeneous medium (lhm)
include("lhm.jl")
using .OpenSWPCLHM: Layer1D, LHMModel, write_lhm!, read_lhm
export Layer1D, LHMModel, write_lhm!, read_lhm

# Setup configuration of a model run
include("openswpc_input.jl")
using .OpenSWPCInput: OpenSWPCConfig, write_input!
export OpenSWPCConfig, write_input!

# run the code 
include("run.jl")   
export swpc_path, swpc_cmd, run_swpc


end # module OpenSWPC
