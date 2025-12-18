using Printf

# few utility functions
_qs(s::AbstractString) = "'" * String(s) * "'"
_num(x::Real) = @sprintf("%g", x)
_deg(x::Real) = @sprintf("%.1f°", x)

# flip matrix upside down
flip_ud(A::AbstractArray) = reverse(A, dims=ndims(A))

# determine if file is netcdf
isnetcdf(file::AbstractString) = endswith(lowercase(file), ".nc")


