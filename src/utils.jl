using Printf

# few utility functions
_qs(s::AbstractString) = "'" * String(s) * "'"
_num(x::Real) = @sprintf("%g", x)
_deg(x::Real) = @sprintf("%.1f°", x)

