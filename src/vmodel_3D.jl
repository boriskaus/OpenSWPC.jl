using Printf
export VelocityModel3D, write_velocity_model_3d!, read_velocity_model_3d, velocity_model_to_cartdata, cartdata_to_velocity_model, interpolate_velocity_model, phase_to_velocity_model
using GeophysicalModelGenerator

"""
    VelocityModel3D(nx, ny, nz, dx, dy, dz, x0, y0, z0, vp, vs, ρ, μ, λ, qp, qs)

Three-dimensional velocity model container.
- `nx`, `ny`, `nz`: number of grid points in x, y, z directions
- `x`, `y`, `z`: coordinates in km
- `x0`, `y0`, `z0`: origin coordinates
- `vp`, `vs`: 3D arrays of P-wave and S-wave velocities (km/s)
- `ρ`: 3D array of density (g/cm^3)
- `μ`, `λ`: 3D arrays of Lamé parameters (GPa)
- `qp`, `qs`: 3D arrays of quality factors (dimensionless)
"""
struct VelocityModel3D <: AbstractVelocityModel
    nx::Int
    ny::Int
    nz::Int
    x::Array{Float64,3}
    y::Array{Float64,3}
    z::Array{Float64,3}
    vp::Array{Float64,3}
    vs::Array{Float64,3}
    ρ::Array{Float64,3}
    μ::Array{Float64,3}
    λ::Array{Float64,3}
    qp::Array{Float64,3}
    qs::Array{Float64,3}
    function velocityModel3D(nx::Int, ny::Int, nz::Int, x::Array{Float64,3}, y::Array{Float64,3}, z::Array{Float64,3}, vp::Array{Float64,3}, vs::Array{Float64,3}, ρ::Array{Float64,3}, μ::Array{Float64,3}, λ::Array{Float64,3}, qp::Array{Float64,3}, qs::Array{Float64,3})
        new(nx, ny, nz, x, y, z, vp, vs, ρ, μ, λ, qp, qs)
    end
end

# Constructor without μ and λ - calculates them from vp, vs, ρ
function VelocityModel3D(nx::Int, ny::Int, nz::Int, 
                        x::Array{Float64,3}, y::Array{Float64,3}, z::Array{Float64,3}, 
                        vp::Array{Float64,3}, vs::Array{Float64,3}, ρ::Array{Float64,3}, 
                        qp::Array{Float64,3}, qs::Array{Float64,3})
    μ = ρ .* vs.^2
    λ = ρ .* (vp.^2 .- 2.0 .* vs.^2)
    return VelocityModel3D(nx, ny, nz, x, y, z, vp, vs, ρ, μ, λ, qp, qs)
end

function Base.show(io::IO, ::MIME"text/plain", m::VelocityModel3D)
    println(io, "VelocityModel3D:")
    println(io, "  Dimensions: $(m.nx) x $(m.ny) x $(m.nz)")
    println(io, "  Extent: x: $(extrema(m.x)...) km, y: $(extrema(m.y)...) km, z: $(extrema(m.z)...) km")
    println(io, "  vp range: $(extrema(m.vp)...) km/s")
    println(io, "  vs range: $(extrema(m.vs)...) km/s")
    println(io, "  ρ range: $(extrema(m.ρ)...) g/cm^3")
    println(io, "  μ range: $(extrema(m.μ)...) GPa")
    println(io, "  λ range: $(extrema(m.λ)...) GPa")
    println(io, "  qp range: $(extrema(m.qp)...)")
    println(io, "  qs range: $(extrema(m.qs)...)")
end

# NOTE TO SELF: rewrite the OpenSWPC code to ignore the header lines

"""
    write_velocity_model_3d!(filename::String, model::VelocityModel3D)

    Write the 3D velocity model to a file that can be read by OpenSWPC using the 'user' defined input file.

"""
function write_velocity_model_3D!(model::VelocityModel3D; filename="coordinates.txt")
    open(filename, "w") do io
        # Write header information
        write(io, "# Velocity Model 3D\n")
        write(io, "# x y z vp vs ρ μ λ qp qs\n")
        write(io, "# Model dimensions: $(model.nx) $(model.ny) $(model.nz)\n")
        write(io, "model extent: $(extrema(model.x) ...) km, $(extrema(model.y)...) km, $(extrema(model.z)...) km\n")
        println(io, "# --------------------------------------------------------")

        # Write material parameter data
        for k in 1:model.nz
            for j in 1:model.ny
                for i in 1:model.nx
                    write(io, "$(model.x[i,j,k]) $(model.y[i,j,k]) $(model.z[i,j,k]) $(model.vp[i,j,k]) $(model.vs[i,j,k]) $(model.ρ[i,j,k]) $(model.μ[i,j,k]) $(model.λ[i,j,k]) $(model.qp[i,j,k]) $(model.qs[i,j,k])\n")
                end
            end
        end
    end
end


function read_velocity_model_3d(filename::String, nx::Int, ny::Int, nz::Int)
    data = readdlm(filename, Float64; skipstart=5)  # Skip header
    x  = reshape(data[:, 1], nx, ny, nz)
    y  = reshape(data[:, 2], nx, ny, nz)
    z  = reshape(data[:, 3], nx, ny, nz)
    vp = reshape(data[:, 4], nx, ny, nz)
    vs = reshape(data[:, 5], nx, ny, nz)
    ρ  = reshape(data[:, 6], nx, ny, nz)
    μ  = reshape(data[:, 7], nx, ny, nz)
    λ  = reshape(data[:, 8], nx, ny, nz)
    qp = reshape(data[:, 9], nx, ny, nz)
    qs = reshape(data[:, 10], nx, ny, nz)

    return VelocityModel3D(nx, ny, nz, x, y, z, vp, vs, ρ, μ, λ, qp, qs)
end

function velocity_model_to_cartdata(model::VelocityModel3D)
    return CartData(model.x, model.y, model.z, fields=(
        vp = model.vp,
        vs = model.vs,
        ρ  = model.ρ,
        μ  = model.μ,
        λ  = model.λ,
        qp = model.qp,
        qs = model.qs,
    ))
end

"""
    cartdata_to_velocity_model(cart::CartData)

Convert a `CartData` object to a `VelocityModel3D`.

"""
function cartdata_to_velocity_model(cart::CartData)
    nx, ny, nz = size(cart.x)
    @assert haskey(cart.fields, :vp) "CartData missing 'vp' field"
    @assert haskey(cart.fields, :vs) "CartData missing 'vs' field"
    @assert haskey(cart.fields, :ρ)  "CartData missing 'ρ' field"
    @assert haskey(cart.fields, :μ)  "CartData missing 'μ' field"
    @assert haskey(cart.fields, :λ)  "CartData missing 'λ' field"
    @assert haskey(cart.fields, :qp) "CartData missing 'qp' field"
    @assert haskey(cart.fields, :qs) "CartData missing 'qs' field"
    return VelocityModel3D(
        nx, ny, nz,
        cart.x,
        cart.y,
        cart.z,
        cart.fields.vp,
        cart.fields.vs,
        cart.fields.ρ,
        cart.fields.μ,
        cart.fields.λ,
        cart.fields.qp,
        cart.fields.qs,
    )
end

# These functions are probably overkill and a more efficient interpolate function could be used that doesn't require multiple conversions
"""
    interpolate_velocity_model(model::VelocityModel3D; input_file=nothing, input_dimensions=nothing)

    Interpolate the given `VelocityModel3D` to match the grid defined in `input_file` or `input_dimensions`.
    The new dimensions are determined from the input file if provided, otherwise from `input_dimensions` tuple.
    Returns a new `VelocityModel3D` interpolated to the target grid.
"""
function interpolate_velocity_model(model::VelocityModel3D, input_file::String)
    # Determine target grid dimensions
    nx, ny, nz = parse_input_dimensions(input_file)

    x_new = range(minimum(model.x), stop=maximum(model.x), length=nx)
    y_new = range(minimum(model.y), stop=maximum(model.y), length=ny)
    z_new = range(minimum(model.z), stop=maximum(model.z), length=nz)
    X, Y, Z = xyz_grid(x_new, y_new, z_new)

    # transform VelocityModel3D to CartData for interpolation
    vel_cart = velocity_model_to_cartdata(model)

    # interpolate to new grid
    model_int_cart = interpolate_datafields(vel_cart, X, Y, Z)

    # transform back to VelocityModel3D
    model_int = cartdata_to_velocity_model(model_int_cart)

    return model_int  # Return the model as is for now
end

function parse_input_dimensions(input_file::String)
    nx = ny = nz = 0
    open(input_file, "r") do io
        for line in eachline(io)
            if startswith(strip(line), "nx")
                nx = parse(Int, split(line, "=")[2])
            elseif startswith(strip(line), "ny")
                ny = parse(Int, split(line, "=")[2])
            elseif startswith(strip(line), "nz")
                nz = parse(Int, split(line, "=")[2])
            end
        end
    end
    return nx, ny, nz
end

"""
    interpolate_velocity_model(model::VelocityModel3D, input_dimensions::Tuple{Int,Int,Int})

    Interpolate the given `VelocityModel3D` to match the specified grid dimensions.
    The `input_dimensions` tuple should contain (nx, ny, nz).
    Returns a new `VelocityModel3D` interpolated to the target grid.
"""
function interpolate_velocity_model(model::VelocityModel3D, input_dimensions::Tuple{Int,Int,Int})
    nx, ny, nz = input_dimensions

    x_new = LinRange(minimum(model.x), maximum(model.x), nx)
    y_new = LinRange(minimum(model.y), maximum(model.y), ny)
    z_new = LinRange(minimum(model.z), maximum(model.z), nz)
    X, Y, Z = xyz_grid(x_new, y_new, z_new)

    # transform VelocityModel3D to CartData for interpolation
    vel_cart = velocity_model_to_cartdata(model)

    # interpolate to new grid
    model_int_cart = interpolate_datafields(vel_cart, X, Y, Z)

    # transform back to VelocityModel3D
    model_int = cartdata_to_velocity_model(model_int_cart)

    return model_int  # Return the model as is for now
end

"""
    interpolate_velocity_model(model::VelocityModel3D, cfg::OpenSWPCConfig)

    Interpolate the given `VelocityModel3D` to match the grid defined in the `OpenSWPCConfig`.
    Returns a new `VelocityModel3D` interpolated to the target grid.
"""
function interpolate_velocity_model(model::VelocityModel3D, cfg::OpenSWPCConfig)
    nx, ny, nz = cfg.nx, cfg.ny, cfg.nz

    x_new = LinRange(minimum(model.x), maximum(model.x), nx)
    y_new = LinRange(minimum(model.y), maximum(model.y), ny)
    z_new = LinRange(minimum(model.z), maximum(model.z), nz)
    X, Y, Z = xyz_grid(x_new, y_new, z_new)

    # transform VelocityModel3D to CartData for interpolation
    vel_cart = velocity_model_to_cartdata(model)

    # interpolate to new grid
    model_int_cart = interpolate_datafields(vel_cart, X, Y, Z)

    # transform back to VelocityModel3D
    model_int = cartdata_to_velocity_model(model_int_cart)

    return model_int  # Return the model as is for now
end

"""
phase_to_velocity_model(phases::AbstractArray{Int,3},
              x::AbstractArray{Float64,3}, 
              y::AbstractArray{Float64,3}, 
              z::AbstractArray{Float64,3}, 
              rho::AbstractArray{Float64,1}, 
              vp::AbstractArray{Float64,1}, 
              vs::AbstractArray{Float64,1}, 
              qp::AbstractArray{Float64,1}, 
              qs::AbstractArray{Float64,1};
              mu=nothing, lambda=nothing)
    
    Function to convert phase model from LaMEM or GMG to SWPC format.
    phases: 3D array of phase numbers
    x, y, z: 3D arrays of coordinates
    rho, vp, vs, qp, qs: 1D arrays mapping phase number to properties
    mu, lambda: optional 1D arrays for shear and Lamé's first parameter; if not provided, they are computed from vp, vs, and rho.
    
    IMPORTANT: Phase 0 is treated as air above surface (z >= 0) and water below surface (z < 0).
"""
function phase_to_velocity_model(phases::AbstractArray{Int,3}, x::AbstractArray{Float64,3}, y::AbstractArray{Float64,3}, z::AbstractArray{Float64,3}, rho::AbstractArray{Float64,1}, vp::AbstractArray{Float64,1}, vs::AbstractArray{Float64,1}, qp::AbstractArray{Float64,1}, qs::AbstractArray{Float64,1}; mu=nothing, lambda=nothing)

    # initialize arrays
    vp3d    = similar(phases, Float64)
    vs3d    = similar(phases, Float64)
    rho3d   = similar(phases, Float64)
    mu3d    = similar(phases, Float64)
    lam3d   = similar(phases, Float64)
    qp3d    = similar(phases, Float64)
    qs3d    = similar(phases, Float64)

    # Initialize 3D arrays for material properties
    for i in eachindex(phases)
        if phases[i] == 0 && z[i] >= 0.0
            # default air properties
            rho3d[i] = 0.001
            mu3d[i]  = 0.0
            lam3d[i] = 0.0
            qp3d[i]  = 10.0
            qs3d[i]  = 10.0
            vp3d[i]  = 0.0
            vs3d[i]  = 0.0
        elseif phases[i] == 0 && z[i] < 0.0
            # default water properties
            vp3d[i]  = 1.5
            vs3d[i]  = 0.0 
            rho3d[i] = 1.0
            mu3d[i]  = 0.0
            lam3d[i] = rho3d[i] * vp3d[i]^2
            qp3d[i]  = 1e6
            qs3d[i]  = 1e6
        else
            rho3d[i] = rho[phases[i]]
            qp3d[i]  = qp[phases[i]]
            qs3d[i]  = qs[phases[i]]
            vp3d[i]  = vp[phases[i]]
            vs3d[i]  = vs[phases[i]]

            if isnothing(mu)
                mu3d[i]  = rho3d[i] * vs3d[i]^2
            else
                mu3d[i]  = mu[phases[i]]
            end

            if isnothing(lambda)
                lam3d[i]  = rho3d[i] * (vp3d[i]^2 - 2.0* vs3d[i]^2)
            else
                lam3d[i]  = lambda[phases[i]]
            end
        end
    end
    return VelocityModel3D(size(phases, 1), size(phases, 2), size(phases, 3), x, y, z, vp3d, vs3d, rho3d, mu3d, lam3d, qp3d, qs3d)
end