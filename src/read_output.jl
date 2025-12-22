# reads output files of a simulation back into julia
# This assumes data to be saved in netcdf format
using NCDatasets, GeophysicalModelGenerator

export read_xy_slice, read_yz_slice, read_xz_slice, read_volume, movie_slice

 
slice_xy(array::Array{T,2}) where T = reshape(array[:,:], size(array[:,:])...,1)
slice_xz(array::Array{T,2}) where T = reshape(array[:,:], size(array,1),1,size(array,2))
slice_yz(array::Array{T,2}) where T = reshape(array[:,:], 1,size(array,1),size(array,2))

create_tuple_field(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,1}) where T = NamedTuple{(Symbol(name),)}((flip_ud(field[:]),))
create_tuple_field(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,2}) where T = NamedTuple{(Symbol(name),)}((flip_ud(field[:,:]),))
create_tuple_field(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,3}) where T = NamedTuple{(Symbol(name),)}((flip_ud(field[:,:,:]),))

create_tuple_field_slice_xy(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,2}) where T = NamedTuple{(Symbol(name),)}((slice_xy(field[:,:]),))
create_tuple_field_slice_xy(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,3}, itime) where T = NamedTuple{(Symbol(name),)}((slice_xy(field[:,:,itime]),))
create_tuple_field_slice_xz(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,2}) where T = NamedTuple{(Symbol(name),)}((flip_ud(slice_xz(field[:,:])),))
create_tuple_field_slice_xz(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,3}, itime) where T = NamedTuple{(Symbol(name),)}((flip_ud(slice_xz(field[:,:,itime])),))
create_tuple_field_slice_yz(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,2}) where T = NamedTuple{(Symbol(name),)}((flip_ud(slice_yz(field[:,:])),))
create_tuple_field_slice_yz(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,3}, itime) where T = NamedTuple{(Symbol(name),)}((flip_ud(slice_yz(field[:,:,itime])),))

create_tuple_field_tstep(name::String, field::NCDatasets.CommonDataModel.CFVariable{T,4}, itime) where T = NamedTuple{(Symbol(name),)}((flip_ud(field[:,:,:,itime]),))

create_nt(name::String, data) = NamedTuple{(Symbol(name),)}((data,))
remove_field(nt::NamedTuple, k::Symbol) = Base.structdiff(nt, (; k => nothing))

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
        out = CartData(X,Y,-flip_ud(Z), fields)  
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
    out = CartData(X,Y,-flip_ud(Z), fields)  
   
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
    out = CartData(X,Y,flip_ud(Z), fields)  
   
    return out, t[timestep]
end

"""
    data, t  = read_volume(file::AbstractString; timestep::Int=1)
extracts a y,z slides at given x0 from netcdf file
"""
function read_volume(file::AbstractString; timestep::Int=1)
    @assert isnetcdf(file) "File $file is not a NetCDF file"
    ds = NCDataset(file)
    @assert haskey(ds, "z") "File $file does not contain variable 'z'"
    @assert haskey(ds, "y") "File $file does not contain variable 'y'"
    @assert haskey(ds, "x") "File $file does not contain variable 'x'"

    t  = ds["t"][:]
    fields = extract_fields_slice(ds, timestep, ("x","y","z"))
    
    X,Y,Z = xyz_grid(ds["x"][:], ds["y"][:], ds["z"][:]);
    out = CartData(X,Y,-flip_ud(Z), fields)  
   
    return out, t[timestep]
end

# extract slices and converts them in a format that GMG needs
function extract_fields_slice(ds::NCDataset, timestep::Int, types=("x","y"))
    fields = (;)
    for varname in keys(ds)
        # add all 2D fields
        if !(varname in types)
            unit = strip(get(ds[varname].var.attrib, "units", ""))
            unit = replace(unit, "0^" => "e", " "=>"")
            unit = replace(unit, "^" => "")
            save_name = varname*"_"*unit
            if length(size(ds[varname])) == 2
                if types == ("x","y")
                    fields = merge(fields, create_tuple_field_slice_xy(save_name, ds[varname]))
                elseif types == ("x","z")
                    fields = merge(fields, create_tuple_field_slice_xz(save_name, ds[varname]))
                elseif types == ("y","z")
                    fields = merge(fields, create_tuple_field_slice_yz(save_name, ds[varname]))
                else
                    error("Unsupported variable dimension for types $types")
                end
            elseif length(size(ds[varname])) == 3
                if types == ("x","y")
                    fields = merge(fields, create_tuple_field_slice_xy(save_name, ds[varname], timestep))
                elseif types == ("x","z")
                    fields = merge(fields, create_tuple_field_slice_xz(save_name, ds[varname], timestep))
                elseif types == ("y","z")
                    fields = merge(fields, create_tuple_field_slice_yz(save_name, ds[varname], timestep))
                elseif types == ("x","y","z")
                    # for 3D volume files saved as
                    fields = merge(fields, create_tuple_field(save_name, ds[varname]))
                else
                    error("Unsupported variable dimension for types $types")
                end    
            elseif length(size(ds[varname])) == 4
                if types == ("x","y","z")
                    fields = merge(fields, create_tuple_field_tstep(save_name, ds[varname], timestep))
                end

            end
        end
    end

    # create velocity vectors
    if  haskey(fields,Symbol("Vx_m/s")) && 
        haskey(fields,Symbol("Vy_m/s")) &&
        haskey(fields,Symbol("Vz_m/s"))

        Vx = fields[Symbol("Vx_m/s")]
        Vy = fields[Symbol("Vy_m/s")]
        Vz = fields[Symbol("Vz_m/s")]

        V_nt   = create_nt("V_m/s",(Vx,Vy,Vz))
        fields = merge(fields, V_nt)
        fields = remove_field(fields, Symbol("Vx_m/s"))
        fields = remove_field(fields, Symbol("Vy_m/s"))
        fields = remove_field(fields, Symbol("Vz_m/s"))
    end

    # create displacement vectors
    if  haskey(fields,Symbol("Ux_m")) && 
        haskey(fields,Symbol("Uy_m")) &&
        haskey(fields,Symbol("Uz_m"))

        Ux = fields[Symbol("Ux_m")]
        Uy = fields[Symbol("Uy_m")]
        Uz = fields[Symbol("Uz_m")]

        U_nt   = create_nt("U_m",(Ux,Uy,Uz))
        fields = merge(fields, U_nt)
        fields = remove_field(fields, Symbol("Ux_m"))
        fields = remove_field(fields, Symbol("Uy_m"))
        fields = remove_field(fields, Symbol("Uz_m"))
    end


    return fields
end


"""
    movie_slice(file::AbstractString; x0=0.0, y0=0,z0=0, slice=:xy)

reads a OpenSWPC netcdf file and saves a movie of slices at given `x0`,`y0` or `z0`
You need to indicate which slice you want (`:xz`,`:yz`, `:xy`, or `:xyz`)

"""
function movie_slice(file::AbstractString; x0=0.0, y0=0, z0=0, slice=:xy)
    ds  = NCDataset(file)
    t_v = ds["t"][:]
    if t_v[1]>1e10
        t_v[1] = 0      # sometimes, there appears to be a bug here
    end

    pvd_file = file
    pvd_file = replace(pvd_file, "." => "_")
    pvd_file = pvd_file * ".pvd"
    pvd_dir  = "slices_paraview"
    movie = movie_paraview(name=pvd_file, Initialize=true)
    for (itime,t) in enumerate(t_v)    

        if slice == :xy
            dat, _ = read_xy_slice(file,z0, timestep=itime)
        elseif slice == :xz
            dat, _ = read_xz_slice(file,y0, timestep=itime)
        elseif slice == :yz 
            dat, _ = read_yz_slice(file,x0, timestep=itime)
        elseif slice == :xyz 
            dat, _ = read_volume(file, timestep=itime)
        else
            error("Unsupported slice type: $slice. Supported are :xy, :xz, :yz, :xyz")
        end
        
        name = file*".paraview"*string(itime)*".vts"
  
        movie= write_paraview(dat, name,pvd=movie, time=t, directory=pvd_dir)
    end
    movie_paraview(pvd=movie, Finalize=true);

    println("Movie saved to ", pvd_file)
    return pvd_file
end




"""
    movie_slice(cfg::OpenSWPCConfig)

Creates movies of all slices specified in the configuration `cfg`.
"""
function movie_slice(cfg::OpenSWPCConfig)
    curdir = pwd()
    cd(cfg.odir)

    listfiles = String[]
    if cfg.xy_ps_sw
        fname = movie_slice("swpc.3d.xy.ps.nc", z0=cfg.z0_xy, slice=:xy)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.xz_ps_sw
        fname = movie_slice("swpc.3d.xz.ps.nc", y0=cfg.y0_xz, slice=:xz)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.yz_ps_sw
        fname = movie_slice("swpc.3d.yz.ps.nc", x0=cfg.x0_yz, slice=:yz)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.fs_ps_sw
        fname = movie_slice("swpc.3d.fs.ps.nc", slice=:xy)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.ob_ps_sw
        fname = movie_slice("swpc.3d.ob.ps.nc", slice=:xy)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.xy_v_sw
        fname = movie_slice("swpc.3d.xy.v.nc", z0=cfg.z0_xy, slice=:xy)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.xz_v_sw
        fname = movie_slice("swpc.3d.xz.v.nc", y0=cfg.y0_xz, slice=:xz)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.yz_v_sw
        fname = movie_slice("swpc.3d.yz.v.nc", x0=cfg.x0_yz, slice=:yz)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.xy_u_sw
        fname = movie_slice("swpc.3d.xy.u.nc", z0=cfg.z0_xy, slice=:xy)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.xz_u_sw
        fname = movie_slice("swpc.3d.xz.u.nc", y0=cfg.y0_xz, slice=:xz)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.yz_u_sw
        fname = movie_slice("swpc.3d.yz.u.nc", x0=cfg.x0_yz, slice=:yz)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.vol_v_sw
        fname = movie_slice("swpc.3d.vol.v.nc", slice=:xyz)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.vol_u_sw
        fname = movie_slice("swpc.3d.vol.u.nc", slice=:xyz)
        push!(listfiles, joinpath(cfg.odir,fname))
    end
    if cfg.vol_ps_sw
        fname = movie_slice("swpc.3d.vol.ps.nc", slice=:xyz)
        push!(listfiles, joinpath(cfg.odir,fname))
    end

    cd(curdir)
    for l in listfiles
        println("Movie-file (*.pvd) created for: ", l)
    end

    return nothing
end