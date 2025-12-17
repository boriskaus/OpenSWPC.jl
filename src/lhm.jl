using Printf
export Layer1D, LHMModel, write_lhm!, read_lhm


"""
    Layer1D(z; vp, vs, rho, qp=200, qs=200)

One-dimensional homogeneous layer properties at depth `z` (km).
- `vp`, `vs`: km/s
- `rho`: g/cm^3
- `qp`, `qs`: quality factors (dimensionless)§
"""
struct Layer1D <: AbstractVelocityModel
    z::Float64
    vp::Float64
    vs::Float64
    rho::Float64
    qp::Float64
    qs::Float64
    function Layer1D(z::Real; vp::Real, vs::Real, rho::Real, qp::Real=200, qs::Real=200)
        new(Float64(z), Float64(vp), Float64(vs), Float64(rho), Float64(qp), Float64(qs))
    end
end

"""
    LHMModel(layers::Vector{Layer1D})

Container for a stack of `Layer1D` from shallow to deep.

Example:
====
```julia
julia> m = LHMModel([
                      Layer1D(0.0; vp=5.5, vs=3.2, rho=2.7, qp=200, qs=150),
                      Layer1D(5.0; vp=6.2, vs=3.6, rho=2.8, qp=300, qs=200),
                    ])
```
"""
struct LHMModel <: AbstractVelocityModel
    layers::Vector{Layer1D}
    function LHMModel(layers::Vector{Layer1D})
        isempty(layers) && error("Provide at least one layer")
        # ensure ordered by depth
        new(layers)
    end
end

#_num(x::Real) = @sprintf("%g", x)

function Base.show(io::IO, ::MIME"text/plain", l::Layer1D)
    print(io, "Layer1D(")
    print(io, "z=", _num(l.z), " km, ")
    print(io, "vp=", _num(l.vp), " km/s, vs=", _num(l.vs), " km/s, rho=", _num(l.rho), " g/cm^3, ")
    print(io, "Qp=", _num(l.qp), ", Qs=", _num(l.qs))
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", m::LHMModel)
    println(io, "LHMModel:")
    for (i, l) in enumerate(m.layers)
        print(io, "                         [", i, "] ")
        show(io, MIME"text/plain"(), l)
        println(io)
    end
end

"""
    write_lhm!(path::AbstractString, model::LHMModel)

Write a layered homogeneous medium file. The format is:
  # z_top z_bot vp vs rho qp qs
for each layer, matching OpenSWPC 'lhm' expectations.
"""
function write_lhm!(path::AbstractString, model::LHMModel)
    open(path, "w") do io
        println(io, "# Layered homogeneous medium (lhm)")
        println(io, "# depth  rho(g/cm^3)  vp(km/s)   vs(km/s)     Qp      Qs")
        println(io, "# --------------------------------------------------------")
        for l in model.layers
            println(io, @sprintf("%7.3f    %6.3f      %6.3f     %6.3f    %5.0f    %5.0f",
                l.z, l.rho, l.vp, l.vs, l.qp, l.qs))
        end
    end
    return nothing
end

"""
    read_lhm(path::AbstractString) -> LHMModel

Parse an lhm-style file into a `LHMModel`.
"""
function read_lhm(path::AbstractString)
    layers = Layer1D[]
    for ln in eachline(path)
        s = strip(ln)
        isempty(s) && continue
        startswith(s, '#') && continue
        toks = split(s)
        length(toks) >= 7 || continue
        try
            zt, zb, vp, vs, rho, qp, qs = parse.(Float64, toks[1:7])
            push!(layers, Layer1D(zt, zb; vp=vp, vs=vs, rho=rho, qp=qp, qs=qs))
        catch
            # skip bad line
        end
    end
    return LHMModel(layers)
end
