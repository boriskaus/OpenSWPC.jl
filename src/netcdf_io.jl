using NCDatasets
export read_netcdf, write_netcdf

"""
    data = read_netcdf(ncfile::String; xbeg=0.0, ybeg=0.0, zbeg=0.0)

Reads a 3D NetCDF file used in OpenSWPC and returns a `CartData` object.
Note that the z-coordinates are flipped to match the depth convention in GMG (positive up).
"""
function read_netcdf(ncfile::String; xbeg=0.0, ybeg=0.0, zbeg=0.0)
    ds = NCDataset(ncfile)
    
    x = ds["x"][:] .+ xbeg
    y = ds["y"][:] .+ ybeg
    z = ds["z"][:] .+ zbeg

    field_names = keys(ds)
    field_names = setdiff(field_names, ["x","y","z"])

    fields = (;)
    for varname in field_names
        fields = merge(fields, OpenSWPC.create_tuple_field(varname, ds[varname]))
    end
    X,Y,Z = xyz_grid(x,y,-OpenSWPC.flip_ud(z));

    return CartData(X,Y,Z, fields)
end


"""
    data = read_netcdf(ncfile::String, cfg::OpenSWPC.OpenSWPCConfig)   
Reads a 3D NetCDF file used in OpenSWPC and returns a `CartData` object, consistent with the model parameters in `cfg`.
"""
read_netcdf(ncfile::String, cfg::OpenSWPC.OpenSWPCConfig) = 
    read_netcdf(ncfile; xbeg=cfg.xbeg, ybeg=cfg.ybeg, zbeg=cfg.zbeg)

"""
    write_netcdf(data::CartData, filename="output.nc")

Writes a `CartData` object to a NetCDF file as used in OpenSWPC.
Note that the z-coordinates are flipped to match the depth convention in OpenSWPC (positive down).
"""
function write_netcdf(data::CartData, filename="output.nc")
    ds = NCDataset(filename,"c")

    x = Float32.(data.x.val[:,1,1])
    y = Float32.(data.y.val[1,:,1])
    z = Float32.(flip_ud(data.z.val[1,1,:]))

    # Define the dimension "x","y","z" with the size 100 and 110 resp.
    defDim(ds,"x",length(x))
    defDim(ds,"y",length(y))
    defDim(ds,"z",length(z))

    # Define the variables temperature with the attribute units
    defVar(ds,"x",x,("x",), attrib=Dict("units"=>"km","long_name"=>"x"))
    defVar(ds,"y",y,("y",), attrib=Dict("units"=>"km","long_name"=>"y"))
    defVar(ds,"z",z,("z",), attrib=Dict("units"=>"km","long_name"=>"z"))

    # save fields
    for (varname, field) in pairs(data.fields)
        defVar(ds, String(varname), flip_ud(Float32.(field)), ("x","y","z"))
    end
    close(ds)

    return nothing
end

