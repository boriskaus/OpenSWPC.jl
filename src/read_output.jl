# reads output files of a simulation back into julia
# This assumes data to be saved in netcdf format
using NCDatasets, GeophysicalModelGenerator

export read_xy_slice, read_yz_slice, read_xz_slice

slice_xy(array::Array{T,2}) where T = reshape(array[:,:], size(array[:,:])...,1)
slice_xz(array::Array{T,2}) where T = reshape(array[:,:], size(array,1),1,size(array,2))
slice_yz(array::Array{T,2}) where T = reshape(array[:,:], 1,size(array,1),size(array,2))

isnetcdf(file::AbstractString) = endswith(lowercase(file), ".nc")
create_tuple_field(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,1}) where T = NamedTuple{(Symbol(name),)}((field[:],))
create_tuple_field(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,3}) where T = NamedTuple{(Symbol(name),)}((field[:,:,:],))

create_tuple_field_slice_xy(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,2}) where T = NamedTuple{(Symbol(name),)}((slice_xy(field[:,:]),))
create_tuple_field_slice_xy(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,3}, itime) where T = NamedTuple{(Symbol(name),)}((slice_xy(field[:,:,itime]),))
create_tuple_field_slice_xz(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,2}) where T = NamedTuple{(Symbol(name),)}((slice_xz(field[:,:]),))
create_tuple_field_slice_xz(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,3}, itime) where T = NamedTuple{(Symbol(name),)}((slice_xz(field[:,:,itime]),))
create_tuple_field_slice_yz(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,2}) where T = NamedTuple{(Symbol(name),)}((slice_yz(field[:,:]),))
create_tuple_field_slice_yz(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,3}, itime) where T = NamedTuple{(Symbol(name),)}((slice_yz(field[:,:,itime]),))


"""
    data, t  = read_xy_slice(file::AbstractString, depth=0.0; timestep::Int=1, cart_data=false) 

Read a horizontal (XY) slice from a NetCDF output file of OpenSWPC simulation.
`depth` specifies the depth (Z) value for the slice.
`timestep` specifies the time index to read for 3D fields (default: 1).
If `cart_data` is true, it returns `CartData` data instead of `GeoData`.

# Example
===
```julia
julia> dir = "cf_swp_layers"
julia> files = readdir(dir)
9-element Vector{String}:
 "swpc.3d.ob.ps.nc"
 "swpc.3d.ob.u.nc"
 "swpc.3d.ob.v.nc"
 "swpc.3d.xy.ps.nc"
 "swpc.3d.xz.ps.nc"
 "swpc.3d.xz.u.nc"
 "swpc.3d.xz.v.nc"
 "swpc.3d.yz.ps.nc"
 "wav"
julia> file = joinpath(dir, files[1])
julia> data = read_xy_slice(file, -2.5, timestep=5, cart_data=true)
```
"""
function read_xy_slice(file::AbstractString, depth=0.0; timestep::Int=1, cart_data=true)
    ds = NCDataset(file)
    @assert isnetcdf(file) "File $file is not a NetCDF file"
    @assert haskey(ds, "y") "File $file does not contain variable 'y'"

    x  = ds["x"][:]
    y  = ds["y"][:]
    t  = ds["t"][:]

    lon = ds["lon"][:,:]
    lat = ds["lat"][:,:]

    fields = extract_fields_slice(ds, timestep, ("x","y"))

    Lon = zeros(size(lon,1), size(lon,2),1)
    Lat = zeros(size(lon,1), size(lon,2),1)
    Dep = fill(depth, size(lon,1), size(lon,2),1)
    Lon[:,:,1]  = lon[:,:]
    Lat[:,:,1]  = lat[:,:]
    out = GeoData(Lon,Lat,Dep, fields)  
    if cart_data
        X,Y,Z = xyz_grid(x,y,depth);
        out = CartData(X,Y,Z, fields)  
    end

    return out, t[timestep]
end

"""
    data, t  = read_xz_slice(file::AbstractString, y0=0.0; timestep::Int=1)
extracts a x,z slides at given y0 from netcdf file
"""
function read_xz_slice(file::AbstractString, y0=0.0; timestep::Int=1)
    ds = NCDataset(file)
    @assert isnetcdf(file) "File $file is not a NetCDF file"
    @assert haskey(ds, "z") "File $file does not contain variable 'z'"

    t  = ds["t"][:]
    fields = extract_fields_slice(ds, timestep, ("x","z"))
    
    X,Y,Z = xyz_grid(ds["x"][:],y0,ds["z"][:]);
    out = CartData(X,Y,Z, fields)  
   
    return out, t[timestep]
end

"""
    data, t  = read_yz_slice(file::AbstractString, x0=0.0; timestep::Int=1)
extracts a y,z slides at given x0 from netcdf file
"""
function read_yz_slice(file::AbstractString, x0=0.0; timestep::Int=1)
    ds = NCDataset(file)
    @assert isnetcdf(file) "File $file is not a NetCDF file"
    @assert haskey(ds, "z") "File $file does not contain variable 'z'"
    @assert haskey(ds, "y") "File $file does not contain variable 'y'"

    t  = ds["t"][:]
    fields = extract_fields_slice(ds, timestep, ("y","z"))
    
    X,Y,Z = xyz_grid(x0, ds["y"][:], ds["z"][:]);
    out = CartData(X,Y,Z, fields)  
   
    return out, t[timestep]
end

# extract slices and converts them in a format that GMG needs
function extract_fields_slice(ds::NCDataset, timestep::Int, types=("x","y"))
    fields = (;)
    for varname in keys(ds)
        # add all 2D fields
        if !(varname in types)
            if length(size(ds[varname])) == 2
                if types == ("x","y")
                    fields = merge(fields, create_tuple_field_slice_xy(varname, ds[varname]))
                elseif types == ("x","z")
                    fields = merge(fields, create_tuple_field_slice_xz(varname, ds[varname]))
                elseif types == ("y","z")
                    fields = merge(fields, create_tuple_field_slice_yz(varname, ds[varname]))
                else
                    error("Unsupported variable dimension for types $types")
                end
            elseif length(size(ds[varname])) == 3
                if types == ("x","y")
                    fields = merge(fields, create_tuple_field_slice_xy(varname, ds[varname], timestep))
                elseif types == ("x","z")
                    fields = merge(fields, create_tuple_field_slice_xz(varname, ds[varname], timestep))
                elseif types == ("y","z")
                    fields = merge(fields, create_tuple_field_slice_yz(varname, ds[varname], timestep))
                end    
            end
        end
    end

    return fields
end