using Printf
export UniformVelocityModel


"""
    UniformVelocityModel(vp0, vs0, rho0; qp0=200, qs0=200, topo0=0)

Uniform (homogeneous) velocity model parameters used when `vmodel = "uni"`.

- vp0: P-wave velocity [km/s]
- vs0: S-wave velocity [km/s]
- rho0: Mass density [g/cm^3]
- qp0: Qp (dimensionless)
- qs0: Qs (dimensionless)
- topo0: Topography depth [km]; if > 0, seawater is filled from z=0 down to this depth.
"""
struct UniformVelocityModel <: AbstractVelocityModel
    vp0::Float64
    vs0::Float64
    rho0::Float64
    qp0::Float64
    qs0::Float64
    topo0::Float64
    function UniformVelocityModel(vp0::Real, vs0::Real, rho0::Real; qp0::Real=200, qs0::Real=200, topo0::Real=0)
        new(Float64(vp0), Float64(vs0), Float64(rho0), Float64(qp0), Float64(qs0), Float64(topo0))
    end
end

function Base.show(io::IO, ::MIME"text/plain", m::UniformVelocityModel)
    print(io, "UniformVelocityModel(")
    print(io, "vp0=", _num(m.vp0), " km/s, vs0=", _num(m.vs0), " km/s, rho0=", _num(m.rho0), " g/cm^3, ")
    print(io, "Qp=", _num(m.qp0), ", Qs=", _num(m.qs0), "")
    print(io, ", topo0=", _num(m.topo0), " km")
    print(io, ")")
end
