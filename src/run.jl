# Run OpenSWPC simulations via MPI

"""
	swpc_path() -> String

Return the absolute path to the OpenSWPC `swpc_3d` executable provided by `OpenSWPC_jll`.
"""
swpc_path() = OpenSWPC_jll.swpc_3d()

"""
	cmd = swpc_cmd( np::Integer; 
                    input::AbstractString="input.dat", 
                    swpc_psv=false,
                    swpc_sh= false,
                    swpc_3d=true,
                    extra_args::Vector{String}=String[], 
                    mpiexec_cmd::Union{Nothing,AbstractString}=nothing)

Build a `Cmd` that runs `swpc_3d` under MPI with `np` ranks using the text input file `input`.
If `mpiexec_cmd` is not provided, it tries `MPI.mpiexec()` when `MPI.jl` is installed, then falls back to `ENV["MPIEXEC"]`, otherwise `mpiexec` on PATH.
You can pass additional CLI arguments via `extra_args` and environment variables via `env`.
"""
function swpc_cmd(np::Integer; 
                input::AbstractString="input.dat", 
                swpc_psv=false,
                swpc_sh= false,
                swpc_3d=true,
                extra_args::Vector{String}=String[], 
                mpiexec_cmd::Union{Nothing,AbstractString}=nothing)

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
    merged = vcat(OpenSWPC.OpenSWPC_jll.LIBPATH[], MPI_LIBPATH[])
    mpirun = addenv(mpiexec, key => join(merged, pathsep))

    # error catching
    if !swpc_3d && !swpc_psv && !swpc_sh
        error("At least one of swpc_3d, swpc_psv, or swpc_sh must be true")
    elseif sum(Bool[swpc_3d, swpc_psv, swpc_sh]) > 1
        error("Only one of swpc_3d, swpc_psv, or swpc_sh can be true at a time")
    end

    if swpc_3d
        # 3D code
        cmd = `$(mpirun) -n $np $(OpenSWPC_jll.swpc_3d().exec)  -i $(input) $(extra_args...)`
    elseif swpc_psv
        cmd = `$(mpirun) -n $np $(OpenSWPC_jll.swpc_psv().exec)  -i $(input) $(extra_args...)`
    elseif swpc_sh
        cmd = `$(mpirun) -n $np $(OpenSWPC_jll.swpc_sh().exec)  -i $(input) $(extra_args...)`
    end

    return cmd
end

"""
	cmd = run_swpc(np::Integer; kwargs...) 

Run `swpc_3d` under MPI with `np` ranks. Keyword args are passed through to `swpc_cmd`.
Returns the process object after successful completion.

Optional keyword arguments
===
- `input::AbstractString`: path to the SWPC input file (default=`"input.dat"`)
- `swpc_3d::Bool`        : run the 3D code (default=`true`)
- `swpc_ps::Bool`        : run the 2D PSV code (default=`false`)
- `swpc_sh::Bool`        : run the 2D SH code (default=`false`)


"""
function run_swpc(np::Integer; kwargs...)
    @show kwargs
	cmd = swpc_cmd(np; kwargs...)
	run(cmd)
end

