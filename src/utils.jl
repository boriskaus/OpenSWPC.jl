using Printf
using GeophysicalModelGenerator

# few utility functions
_qs(s::AbstractString) = "'" * String(s) * "'"
_num(x::Real) = @sprintf("%g", x)
_deg(x::Real) = @sprintf("%.1f°", x)

# flip matrix upside down
flip_ud(A::AbstractArray) = reverse(A, dims=ndims(A))

# determine if file is netcdf
isnetcdf(file::AbstractString) = endswith(lowercase(file), ".nc")

"""
    to_2D(cs::CartData)

Re-lays a vertical GMG cross-section into the `(nx, 1, nz)` form used for 2D models, with `x`
the distance along the profile (from its first point) and `y = 0`. Accepts the three shapes
`cross_section` produces: fixed y (`Lat_level`, already `(nx, 1, nz)`), fixed x (`Lon_level`,
`(1, ny, nz)`) and oblique profiles (`Start`/`End`, `(n, nz, 1)`).
Sources and stations for such a model are given in profile coordinates.
"""
function to_2D(cs::CartData)
    n = size(cs)
    n[2] == 1 && return cs
    if n[1] == 1            # Lon_level: (1, ny, nz)
        perm = (2, 1, 3)
    elseif n[3] == 1        # Start/End profile: (n, nz, 1)
        perm = (1, 3, 2)
    else
        error("CartData of size $n is not a vertical cross-section")
    end
    p(A::AbstractArray) = permutedims(A, perm)
    p(t::Tuple) = map(p, t)
    x = p(flatten_cross_section(cs))
    fields = NamedTuple{keys(cs.fields)}(map(p, values(cs.fields)))
    return CartData(x, zero(x), p(cs.z.val), fields)
end


