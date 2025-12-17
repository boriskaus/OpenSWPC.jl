export StationLL, write_stations_ll!

using Printf


"""
    StationLL(lon, lat, dep, stnm, zsw)

Station location in geographic coordinates (lon/lat) with depth and placement mode.
- `lon`, `lat`: degrees
- `dep`: depth (km; positive down)
- `stnm`: station code/name
- `zsw`: placement mode, one of `dep`, `fsb`, `obb`, `oba`, `bd0`..`bd9`
"""
struct StationLL
    lon::Float64
    lat::Float64
    dep::Float64
    stnm::String
    zsw::String
end

function Base.show(io::IO, ::MIME"text/plain", s::StationLL)
    print(io, "                         StationLL(")
    print(io, "stnm=", s.stnm, ", lon=", _num(s.lon), ", lat=", _num(s.lat), ", dep=", _num(s.dep), " km, zsw=", s.zsw)
    print(io, ")")
end

"""
    write_stations_ll!(path::AbstractString, stations)

Write a `station.ll`-style file with the `ll` format header. Accepts a single
`StationLL` or a vector of them.
"""
function write_stations_ll!(path::AbstractString, stations)
    if stations isa StationLL
        return write_stations_ll!(path, [stations])
    elseif !(stations isa AbstractVector) || isempty(stations)
        error("Provide a StationLL or a non-empty vector of StationLL")
    end
    all(st -> st isa StationLL, stations) || error("All entries must be StationLL")

    open(path, "w") do io
        println(io, "#                                                    -*- mode:sh -*-")
        println(io, "# stloc.ll")
        println(io)
        println(io, "# station location data by geographical format. ")
        println(io, "# lines starting from '#' and blank lines are omitted. ")
        println(io, "#")
        println(io, "# zsw: controls station depth")
        println(io, "#      'dep': use the depth")
        println(io, "#      'fsb': locate one-grid below from the free surface/sea surface")
        println(io, "#      'obb': locate one-grid below from the groud surface/seafloor")
        println(io, "#      'oba': locate one-grid above from the groud surface/seafloor")
        println(io, "#      'bd{i}' (i=0,...,9) i-th boundary interface")
        println(io, "#")
        println(io, "#     lon       lat    dep     stnm   zsw")
        println(io, "# --------------------------------------")
        for s in stations
            # Match spacing style: lon lat dep stnm 'zsw'
            println(io, @sprintf(" %-10.4f %-10.4f %-6.3f   %-4s   %s", s.lon, s.lat, s.dep, s.stnm, _qs(s.zsw)))
        end
    end
    return nothing
end

