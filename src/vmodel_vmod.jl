using Printf, GeophysicalModelGenerator
export FullVelocityModel


"""
    FullVelocityModel(vmod::CartData)

Sets a CartData velocity model for OpenSWPC.
"""
struct FullVelocityModel <: AbstractVelocityModel
    data::CartData
    function FullVelocityModel(vmod::CartData)
        new(vmod)
    end
end

# Allow seamless conversion from `CartData` to `AbstractVelocityModel`
Base.convert(::Type{AbstractVelocityModel}, vmod::CartData) = FullVelocityModel(vmod)

function Base.show(io::IO, ::MIME"text/plain", m::FullVelocityModel)
    print(io, "FullVelocityModel(")
    print(io, "3D CartData velocity model, size=", size(m.data), ", fields=", keys(m.data.fields),")")
end
