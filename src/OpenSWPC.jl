module OpenSWPC

using OpenSWPC_jll 
using NetCDF_jll
import MPI

include("openswpc_input.jl")
using .OpenSWPCInput: OpenSWPCConfig, write_input!
export OpenSWPCConfig, write_input!


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



"""
	swpc_path() -> String

Return the absolute path to the OpenSWPC `swpc_3d` executable provided by `OpenSWPC_jll`.
"""
swpc_path() = OpenSWPC_jll.swpc_3d()

"""
	swpc_cmd(np::Integer; input::AbstractString="input.dat", extra_args::Vector{String}=String[], mpiexec_cmd::Union{Nothing,AbstractString}=nothing, env::Dict{String,String}=Dict()) -> Cmd

Build a `Cmd` that runs `swpc_3d` under MPI with `np` ranks using the text input file `input`.
If `mpiexec_cmd` is not provided, it tries `MPI.mpiexec()` when `MPI.jl` is installed, then falls back to `ENV["MPIEXEC"]`, otherwise `mpiexec` on PATH.
You can pass additional CLI arguments via `extra_args` and environment variables via `env`.
"""
function swpc_cmd(np::Integer; input::AbstractString="input.dat", extra_args::Vector{String}=String[], mpiexec_cmd::Union{Nothing,AbstractString}=nothing, env::Dict{String,String}=Dict{String,String}())
	isfile(input) || error("Input file not found: $(input)")

	# Resolve mpiexec command
	if mpiexec_cmd === nothing
		mpiexec_cmd = try
			MPI.mpiexec()
		catch
			get(ENV, "MPIEXEC", something(Sys.which("mpiexec"), "mpiexec"))
		end
	end
    key = OpenSWPC.OpenSWPC_jll.JLLWrappers.JLLWrappers.LIBPATH_env
    # Merge OpenSWPC, MPI, and NetCDF library paths so dyld finds dependencies
    netcdf_lib = "/Users/kausb/.julia/artifacts/8b2c8bc12a13efbc55a3e6c0334604e4c367931d/lib"
    merged = vcat(OpenSWPC.OpenSWPC_jll.LIBPATH[], MPI_LIBPATH[], NetCDF_jll.LIBPATH[], [netcdf_lib])
    
    mpirun = addenv(mpiexec, key => join(merged, pathsep))
  
	exe = swpc_path().exec
    @show exe mpiexec_cmd
	#cmd = `$(mpiexec_cmd) -n $np $(exe) $input $(extra_args...)`
    
    cmd = `$(mpirun) -n $np $(OpenSWPC_jll.swpc_3d().exec)  -i $(input) $(extra_args...)`

end

"""
	run_swpc(np::Integer; kwargs...) -> Base.Process

Run `swpc_3d` under MPI with `np` ranks. Keyword args are passed through to `swpc_cmd`.
Returns the process object after successful completion.
"""
function run_swpc(np::Integer; kwargs...)
	cmd = swpc_cmd(np; kwargs...)
    @show cmd
	run(cmd)
end

export swpc_path, swpc_cmd, run_swpc

end # module OpenSWPC
