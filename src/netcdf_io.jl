using NCDatasets
export read_netcdf, write_netcdf

"""
    data = read_netcdf(ncfile::String; xbeg=0.0, ybeg=0.0, zbeg=0.0)

Reads a 2D or 3D NetCDF file used in OpenSWPC and returns a `CartData` object.
Note that the z-coordinates are flipped to match the depth convention in GMG (positive up).
A 2D file (dimensions `x`, `z`) is returned as a `CartData` of size `(nx, 1, nz)` with `y = ybeg`.
"""
function read_netcdf(ncfile::String; xbeg=0.0, ybeg=0.0, zbeg=0.0)
    ds = NCDataset(ncfile)
    is2d = !haskey(ds.dim, "y")
    
    x = ds["x"][:] .+ xbeg
    y = is2d ? [ybeg] : ds["y"][:] .+ ybeg
    z = ds["z"][:] .+ zbeg

    field_names = keys(ds)
    field_names = setdiff(field_names, ["x","y","z"])

    fields = (;)
    for varname in field_names
        field = OpenSWPC.create_tuple_field(varname, ds[varname])
        is2d && (field = map(f -> reshape(f, length(x), 1, length(z)), field))
        fields = merge(fields, field)
    end
    X,Y,Z = xyz_grid(x,y,-OpenSWPC.flip_ud(z));
    close(ds)

    return CartData(X,Y,Z, fields)
end


"""
    data = read_netcdf(ncfile::String, cfg::OpenSWPC.OpenSWPCConfig)   
Reads a 2D or 3D NetCDF file used in OpenSWPC and returns a `CartData` object, consistent with the model parameters in `cfg`.
"""
read_netcdf(ncfile::String, cfg::OpenSWPC.OpenSWPCConfig) = 
    read_netcdf(ncfile; xbeg=cfg.xbeg, ybeg=cfg.ybeg, zbeg=cfg.zbeg)

"""
    write_netcdf(data::CartData, filename="output.nc")

Writes a `CartData` object to a NetCDF file as used in OpenSWPC.
Note that the z-coordinates are flipped to match the depth convention in OpenSWPC (positive down).
A `CartData` of size `(nx, 1, nz)` is written as a 2D file with dimensions `x`, `z`, as read by `swpc_psv` and `swpc_sh`.
"""
function write_netcdf(data::CartData, filename="output.nc")
    ds = NCDataset(filename,"c")
    is2d = size(data)[2] == 1
    dims = is2d ? ("x","z") : ("x","y","z")

    x = Float32.(data.x.val[:,1,1])
    y = Float32.(data.y.val[1,:,1])
    z = Float32.(flip_ud(data.z.val[1,1,:]))

    # Define the dimension "x","y","z" with the size 100 and 110 resp.
    defDim(ds,"x",length(x))
    is2d || defDim(ds,"y",length(y))
    defDim(ds,"z",length(z))

    # Define the variables temperature with the attribute units
    defVar(ds,"x",x,("x",), attrib=Dict("units"=>"km","long_name"=>"x"))
    is2d || defVar(ds,"y",y,("y",), attrib=Dict("units"=>"km","long_name"=>"y"))
    defVar(ds,"z",z,("z",), attrib=Dict("units"=>"km","long_name"=>"z"))

    # save fields
    for (varname, field) in pairs(data.fields)
        field = is2d ? field[:,1,:] : field
        defVar(ds, String(varname), flip_ud(Float32.(field)), dims)
    end
    close(ds)

    return nothing
end
