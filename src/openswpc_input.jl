
using GeophysicalModelGenerator
export OpenSWPCConfig, write_input!, clean

"""
        OpenSWPCConfig

Top-level configuration for preparing an OpenSWPC input file (`input.dat`).
Fields mirror the sections found in the SWPC example input, and can be
overridden via the keyword constructor `OpenSWPCConfig(; kwargs...)`.

Sections and notable fields:

- Control: `title`, `odir`, `ntdec_r`, `strict_mode`.
- Model/Grid Size/Area: 
    - `nx`, `ny`, `nz` : number of grid cells 
    - `dx`, `dy`, `dz`: spacing
    - `nt`, `dt`: timestepping
    domain extents `xbeg`, `ybeg`, `zbeg`, start time `tbeg`, map center
    `clon`, `clat`, rotation `phi`, and Q-constant model freqs `fq_min/max/ref`.
- Snapshot Output: `snp_format`, switches like `xy_ps_sw`, decimation
    `ntdec_s`, and spatial decimations `idec/jdec/kdec`.
- Waveform Output: `sw_wav_*` switches, `ntdec_w`, `st_format`, `fn_stloc`,
    and `wav_format`.
- Earthquake Source: `stf_format`, `stftype`, `fn_stf`, `sdep_fit`, and an
    optional `source::Vector{AbstractSource}` with per-event specifications
    (from the `OpenSWPCSource` module).
- Body/Plane wave modes: `bf_mode`, `pw_mode` with plane-wave geometry
    (`pw_ztop`, `pw_zlen`, `pw_ps`, `pw_strike/dip/rake`).
- Absorbing BC: `abc_type`, `na`, `stabilize_pml`.
- Velocity model: choose among `uniform`/`grd`/`lhm`/`lgm`/`velm` via `vmodel_type` set 
    `is_ocean`, `topo_flatten`, `munk_profile`, `earth_flattening`.

    - Uniform: set when `vmodel=UniformVelocityModel(...)` is passed. See `UniformVelocityModel` for details.
    - Layered homogeneous ('lhm' or `lgm`): `fn_lhm`.
    - GMT grid: `dir_grd`, `fn_grdlst`, `node_grd`.
    - 2D/3D velocity model from file: `dir_velm`, `fn_velm`.
    - Random medium: `dir_rmed`, `fn_grdlst_rmed`, `rhomin`, `fn_rmed0`.
- Checkpoint/Restart: `is_ckp`, `ckpdir`, `ckp_interval`, `ckp_time`, `ckp_seq`.
- Reciprocity Green's Function: `green_*` fields.
- MISC: `stopwatch_mode`, `benchmark_mode`, `ipad/jpad/kpad`.

Usage:

- Start from defaults and override selectively:
    `cfg = OpenSWPCConfig(nx=128, ny=128, nz=128, dt=0.01, title="test")`
- Write an input file:
    `write_input!(cfg, "input.dat")`
- Pretty-print: showing `cfg` in the REPL lists all parameters by section.
"""
mutable struct OpenSWPCConfig
    # Control
    input_file::String
    title::String
    odir::String
    ntdec_r::Int
    strict_mode::Bool

    # Model/Grid Size and Area
    nproc_x::Int
    nproc_y::Int
    nx::Int
    ny::Int
    nz::Int
    nt::Int
    dx::Float64
    dy::Float64
    dz::Float64
    dt::Float64
    vcut::Float64
    xbeg::Float64
    ybeg::Float64
    zbeg::Float64
    tbeg::Float64
    clon::Float64
    clat::Float64
    phi::Float64
    fq_min::Float64
    fq_max::Float64
    fq_ref::Float64

    # Snapshot Output
    snp_format::String
    xy_ps_sw::Bool
    xz_ps_sw::Bool
    yz_ps_sw::Bool
    fs_ps_sw::Bool
    ob_ps_sw::Bool
    xy_v_sw::Bool
    xz_v_sw::Bool
    yz_v_sw::Bool
    fs_v_sw::Bool
    ob_v_sw::Bool
    xy_u_sw::Bool
    xz_u_sw::Bool
    yz_u_sw::Bool
    fs_u_sw::Bool
    ob_u_sw::Bool
    vol_v_sw::Bool
    vol_u_sw::Bool
    vol_ps_sw::Bool

    z0_xy::Float64
    x0_yz::Float64
    y0_xz::Float64
      
    ntdec_s::Int
    idec::Int
    jdec::Int
    kdec::Int

    # Waveform Output
    stations::Vector{StationLL} 
    sw_wav_v::Bool
    sw_wav_u::Bool
    sw_wav_stress::Bool
    sw_wav_strain::Bool
    ntdec_w::Int
    st_format::String
    fn_stloc::String
    wav_format::String
    ntdec_w_prg::Int

    # Earthquake Source
    stf_format::String
    stftype::String
    fn_stf::String
    sdep_fit::String
    source::Vector{<:AbstractSource}  # optional field to hold sources

    # Body force source mode
    bf_mode::Bool

    # Plane wave source mode
    pw_mode::Bool
    pw_ztop::Float64
    pw_zlen::Float64
    pw_ps::String
    pw_strike::Float64
    pw_dip::Float64
    pw_rake::Float64

    # Absorbing Boundary Condition
    abc_type::String
    na::Int
    stabilize_pml::Bool

    # Velocity model
    vmodel::AbstractVelocityModel
    vmodel_type::String
    is_ocean::Bool
    topo_flatten::Bool
    munk_profile::Bool
    earth_flattening::Bool

    # For uniform velocity model
    vp0::Float64
    vs0::Float64
    rho0::Float64
    qp0::Float64
    qs0::Float64
    topo0::Float64

    # For GMT grid file input
    dir_grd::String
    fn_grdlst::String
    node_grd::Int

    # For layered homogeneous medium model
    fn_lhm::String

    # For reading a 2D/3D velocity model from disk
    dir_velm::String
    fn_velm::String

    # For random medium models
    dir_rmed::String
    fn_grdlst_rmed::String
    rhomin::Float64
    fn_rmed0::String

    # Checkpoint/Restart
    is_ckp::Bool
    ckpdir::String
    ckp_interval::Int
    ckp_time::Float64
    ckp_seq::Bool

    # Reciprocity Green's Function Mode
    green_mode::Bool
    green_stnm::String
    green_cmp::String
    green_trise::Float64
    green_bforce::Bool
    green_maxdist::Float64
    green_fmt::String
    fn_glst::String

    # MISC
    stopwatch_mode::Bool
    benchmark_mode::Bool
    ipad::Int
    jpad::Int
    kpad::Int
end

_b(b::Bool) = b ? ".true." : ".false."
_qs(s::String) = "'" * s * "'"

"""
    OpenSWPCConfig(; kwargs...)

Keyword constructor with sensible defaults mirroring `example_input_file.txt`.
Override any field by passing it as a keyword.
"""
function OpenSWPCConfig(; kwargs...)
    defaults = (
        # Control
        input_file = "inputCF.inf",
        title = "swpc",
        odir = "./out",
        ntdec_r = 50,
        strict_mode = false,
        
        # Model/Grid Size and Area
        nproc_x = 2, nproc_y = 2,
        nx = 384, ny = 384, nz = 384, nt = 1000,
        dx = 0.5, dy = 0.5, dz = 0.5, dt = 0.02,
        vcut = 1.5,
        xbeg = -96.0, ybeg = -96.0, zbeg = -10.0, tbeg = 0.0,
        clon = 139.7604, clat = 35.7182, phi = 0.0,
        fq_min = 0.02, fq_max = 2.00, fq_ref = 1.0,
        
        # Snapshot Output
        snp_format = "netcdf",
        xy_ps_sw = false, xz_ps_sw = true, yz_ps_sw = false, fs_ps_sw = false, ob_ps_sw = true,
        xy_v_sw = false, xz_v_sw = true, yz_v_sw = false, fs_v_sw = false, ob_v_sw = true,
        xy_u_sw = false, xz_u_sw = true, yz_u_sw = false, fs_u_sw = false, ob_u_sw = true,
        vol_v_sw = false, vol_u_sw = false, vol_ps_sw = false,
        z0_xy = 7.0, x0_yz = 0.0, y0_xz = 0.0,
        ntdec_s = 5, idec = 2, jdec = 2, kdec = 2,
        
        # Waveform Output
        stations = StationLL[],
        sw_wav_v = true, sw_wav_u = false, sw_wav_stress = false, sw_wav_strain = false,
        ntdec_w = 5, st_format = "ll", fn_stloc = "stations.ll",
        wav_format = "sac", ntdec_w_prg = 0,
        
        # Earthquake Source
        stf_format = "xym0ij", stftype = "kupper", fn_stf = "sourceCF.dat", sdep_fit = "asis",source=AbstractSource[],
        
        # Body force mode
        bf_mode = false,
        
        # Plane wave mode
        pw_mode = false, pw_ztop = 100.0, pw_zlen = 30.0, pw_ps = "p",
        pw_strike = 0.0, pw_dip = 0.0, pw_rake = 0.0,
        
        # Absorbing BC
        abc_type = "pml", na = 20, stabilize_pml = false,
        
        # Velocity model
        vmodel = UniformVelocityModel(6.0, 3.5, 2.7; qp0=150, qs0=80, topo0=0.5),
        vmodel_type = "lhm", is_ocean = false, topo_flatten = false, munk_profile = false, earth_flattening = false,
        
        # Uniform model
        vp0 = 5.0, vs0 = 3.0, rho0 = 2.7, qp0 = 200, qs0 = 200, topo0 = 0,
        
        # GMT grid
        dir_grd = "/vmodel/ejivsm/", fn_grdlst = "grd.lst", node_grd = 0,
        
        # Layered homogeneous medium
        fn_lhm = "lhm.dat",
        
        # 2D/3D velocity model
        dir_velm = "./",   # not used if empty
        fn_velm = "velocity_model.nc",    # not used if empty

        # Random medium
        dir_rmed = "./in/", fn_grdlst_rmed = "grd.lst", rhomin = 1.0, fn_rmed0 = "dummy.ns",
        
        # Checkpoint
        is_ckp = false, ckpdir = "./out/ckp", ckp_interval = 1000000, ckp_time = 1000000.0, ckp_seq = true,
        
        # Green's function
        green_mode = false, green_stnm = "st01", green_cmp = "z", green_trise = 1.0,
        green_bforce = false, green_maxdist = 550.0, green_fmt = "llz", fn_glst = "example/green.lst",
        
        # MISC
        stopwatch_mode = false, benchmark_mode = false,
        ipad = 0, jpad = 0, kpad = 0,
    )

    # Merge kwargs over defaults
    cfg = merge(NamedTuple(defaults), kwargs)
    if isa(cfg.vmodel, UniformVelocityModel)
        vm      = cfg.vmodel
        vmodel_type = "uni"
        vp0 = vm.vp0
        vs0 = vm.vs0
        rho0 = vm.rho0
        qp0 = vm.qp0
        qs0 = vm.qs0
        topo0 = vm.topo0
        cfg = merge(cfg, (;vmodel_type,vp0,vs0,rho0,qp0,qs0,topo0))
    end

    if isa(cfg.vmodel, LHMModel)
        vmodel_type = cfg.vmodel_type
        if vmodel_type != "lhm" && vmodel_type != "lgm"
            vmodel_type = "lgm"
        end
        cfg = merge(cfg, (; vmodel_type))
    end

    if isa(cfg.vmodel, FullVelocityModel)
        vmodel_type = "velm"
        cfg = merge(cfg, (; vmodel_type))
    end

    if  eltype(cfg.source)==SourceLLMWDC
        stf_format = "llmwdc"
        cfg = merge(cfg, (; stf_format))
    elseif eltype(cfg.source)==SourceXYMWDC
        stf_format = "xymwdc"
        cfg = merge(cfg, (; stf_format))
    end


    return OpenSWPCConfig(cfg.input_file,
        cfg.title, cfg.odir, cfg.ntdec_r, cfg.strict_mode,
        cfg.nproc_x, cfg.nproc_y, cfg.nx, cfg.ny, cfg.nz, cfg.nt,
        cfg.dx, cfg.dy, cfg.dz, cfg.dt, cfg.vcut,
        cfg.xbeg, cfg.ybeg, cfg.zbeg, cfg.tbeg,
        cfg.clon, cfg.clat, cfg.phi, cfg.fq_min, cfg.fq_max, cfg.fq_ref,
        cfg.snp_format,
        cfg.xy_ps_sw, cfg.xz_ps_sw, cfg.yz_ps_sw, cfg.fs_ps_sw, cfg.ob_ps_sw,
        cfg.xy_v_sw, cfg.xz_v_sw, cfg.yz_v_sw, cfg.fs_v_sw, cfg.ob_v_sw,
        cfg.xy_u_sw, cfg.xz_u_sw, cfg.yz_u_sw, cfg.fs_u_sw, cfg.ob_u_sw,
        cfg.vol_v_sw, cfg.vol_u_sw, cfg.vol_ps_sw,
        cfg.z0_xy, cfg.x0_yz, cfg.y0_xz,
        cfg.ntdec_s, cfg.idec, cfg.jdec, cfg.kdec,
        cfg.stations,
        cfg.sw_wav_v, cfg.sw_wav_u, cfg.sw_wav_stress, cfg.sw_wav_strain,
        cfg.ntdec_w, cfg.st_format, cfg.fn_stloc, cfg.wav_format, cfg.ntdec_w_prg,
        cfg.stf_format, cfg.stftype, cfg.fn_stf, cfg.sdep_fit, cfg.source,
        cfg.bf_mode,
        cfg.pw_mode, cfg.pw_ztop, cfg.pw_zlen, cfg.pw_ps, cfg.pw_strike, cfg.pw_dip, cfg.pw_rake,
        cfg.abc_type, cfg.na, cfg.stabilize_pml,
        cfg.vmodel, cfg.vmodel_type, cfg.is_ocean, cfg.topo_flatten, cfg.munk_profile, cfg.earth_flattening,
        cfg.vp0, cfg.vs0, cfg.rho0, cfg.qp0, cfg.qs0, cfg.topo0,
        cfg.dir_grd, cfg.fn_grdlst, cfg.node_grd,
        cfg.fn_lhm,
        cfg.dir_velm, cfg.fn_velm,
        cfg.dir_rmed, cfg.fn_grdlst_rmed, cfg.rhomin, cfg.fn_rmed0,
        cfg.is_ckp, cfg.ckpdir, cfg.ckp_interval, cfg.ckp_time, cfg.ckp_seq,
        cfg.green_mode, cfg.green_stnm, cfg.green_cmp, cfg.green_trise, cfg.green_bforce,
        cfg.green_maxdist, cfg.green_fmt, cfg.fn_glst,
        cfg.stopwatch_mode, cfg.benchmark_mode, cfg.ipad, cfg.jpad, cfg.kpad,
    )
end

"""
Compact REPL printing that includes all parameters, grouped by section.
"""
function Base.show(io::IO, ::MIME"text/plain", cfg::OpenSWPCConfig)
    n = (cfg.nx, cfg.ny, cfg.nz)
    sections = [
        ("Control             ", (
            :input_file, :title, :odir, :ntdec_r, :strict_mode
        )),
        ("Grid Size           ", (
            :nx, :ny, :nz,
        )),
        ("Model Domain        ", (
            :dx, :dy, :dz, :vcut,:xbeg, :ybeg, :zbeg,
        )),
        ("Model Coord         ", (
             :clon, :clat, :phi, 
        )),
        ("Frequency           ", (
             :fq_min, :fq_max, :fq_ref
        )),
        ("Timestepping        ", (
            :nt,
            :dt,
            :tbeg,
        )),
        ("Parallelisation     ", (
            :nproc_x, :nproc_y 
        )),
        ("Output common       ", (
            :ntdec_s, :idec, :jdec, :kdec
        )),
        ("Output Free Surface ", (
            :fs_v_sw,  :fs_u_sw, :fs_ps_sw, 
        )),
        ("Output Ocean Bottom ", (
            :ob_v_sw,  :ob_u_sw, :ob_ps_sw
        )),
        ("Output XY Slice     ", (
            :z0_xy, 
            :xy_v_sw,  :xy_u_sw, :xy_ps_sw,
        )),
         ("Output XZ Slice     ", (
            :y0_xz,
            :xz_ps_sw, :xz_v_sw, :yz_v_sw
        )),
          ("Output YZ Slice     ", (
            :x0_yz,
            :yz_ps_sw, :yz_v_sw, :yz_u_sw, 
        )),
        ("Output 3D           ", (
            :vol_v_sw, :vol_u_sw, :vol_ps_sw,
        )),
        #("Waveform Output     ", (
        #    :sw_wav_v, :sw_wav_u, :sw_wav_stress, :sw_wav_strain,
        #    :ntdec_w, :st_format, :fn_stloc, :wav_format, :ntdec_w_prg
        #)),
        ("Body Force Mode     ", (
            :bf_mode
        )),
        ("Plane Wave Mode     ", (
            :pw_mode, :pw_ztop, :pw_zlen, :pw_ps, :pw_strike, :pw_dip, :pw_rake
        )),
        ("Absorbing BCs       ", (
            :abc_type, :na, :stabilize_pml
        )),
        #("Velocity Model      ", (
        #    :vmodel_type, :is_ocean, :topo_flatten, :munk_profile, :earth_flattening
        #)),
        #("Uniform Model       ", (
        #    :vp0, :vs0, :rho0, :qp0, :qs0, :topo0
        #)),
        ("GMT Grid            ", (
            :dir_grd, :fn_grdlst, :node_grd
        )),
        #("Layered Hom. Medium ", (
        #    :fn_lhm
        #)),
        ("Random Medium       ", (
            :dir_rmed, :fn_grdlst_rmed, :rhomin, :fn_rmed0
        )),
        ("Checkpoint/Restart  ", (
            :is_ckp, :ckpdir, :ckp_interval, :ckp_time, :ckp_seq
        )),
        ("Reciproc. Green's F.", (
            :green_mode, :green_stnm, :green_cmp, :green_trise, :green_bforce,
            :green_maxdist, :green_fmt, :fn_glst
        )),
        ("MISC                ", (
            :stopwatch_mode, :benchmark_mode, :ipad, :jpad, :kpad
        )),
        #("Earthquake Source   ", (
        #    :stf_format, :stftype, :fn_stf, :sdep_fit
        #)),

    ]

    println(io, "OpenSWPCConfig:")
    for (section, fields) in sections
        items = String[]
        flist = fields isa Tuple ? fields : (fields,)
        for f in flist
            val = getfield(cfg, f)
            push!(items, string(f, "=", repr(val)))
        end
        println(io, "  ", section, ": ", join(items, ", "))
    end

    # stations
    print(io, "  Waveform Output     : fn_stloc=\"$(cfg.fn_stloc)\", wav_format=\"$(cfg.wav_format)\",ntdec_w_prg=$(cfg.ntdec_w_prg), sw_wav_v/u=($(cfg.sw_wav_v), $(cfg.sw_wav_u)), sw_wav_stress/strain=($(cfg.sw_wav_stress), $(cfg.sw_wav_strain))\n")
      #Waveform Output     : sw_wav_v=true, sw_wav_u=false, sw_wav_stress=false, sw_wav_strain=false, ntdec_w=5, st_format="xy", fn_stloc="./example/stloc.xy", wav_format="sac", ntdec_w_prg=0

    for s in cfg.stations
        show(io, MIME"text/plain"(), s)
        #print(io, "                        ")
        println(io)
    end
    # EQ sources
    print(io, "  Earthquake Sources  : stf_format=\"$(cfg.stf_format)\", stftype=\"$(cfg.stftype)\", fn_stf=\"$(cfg.fn_stf)\", sdep_fit=\"$(cfg.sdep_fit)\", number of sources=$(length(cfg.source))\n")
    for s in cfg.source
        show(io, MIME"text/plain"(), s)
        #print(io, "                        ")
        println(io)
    end

    # velocity model
    println(io, "  Velocity Model      : vmodel_type=\"$(cfg.vmodel_type)\", is_ocean=$(cfg.is_ocean), topo_flatten=$(cfg.topo_flatten), munk_profile=$(cfg.munk_profile), earth_flattening=$(cfg.earth_flattening)")
    if cfg.vmodel_type == "lhm" || cfg.vmodel_type == "lgm"
        println(io, "                        Layered Homogeneous Medium, fn_lhm=\"$(cfg.fn_lhm)\" ")
    elseif cfg.vmodel_type == "velm"
        println(io, "                        Velocity model from file, dir_velm=\"$(cfg.dir_velm)\", fn_velm=\"$(cfg.fn_velm)\" ")
    end
    print(io, "                        ")
    show(io, MIME"text/plain"(), cfg.vmodel)
    println(io)
    

end

"""
    write_input!(cfg::OpenSWPCConfig)
Write the OpenSWPC input file as specified in `cfg.input_file`.
"""
function write_input!(cfg::OpenSWPCConfig)
    path=cfg.input_file
    open(path, "w") do io
        println(io, " !! ----------------------------------------------------------------------- !!")
        println(io, "  !!\n  !!  SWPC input file\n  !!\n !! ----------------------------------------------------------------------- !!\n")

        println(io, "  !! ----------------------------------------------------------------------- !!\n  !! Control\n  !!\n")
        println(io, "  title            = $( _qs(cfg.title) )           !! exe title: used for output filenames")
        println(io, "  odir             = $( _qs(cfg.odir) )          !! output directory")
        println(io, "  ntdec_r          = $(cfg.ntdec_r)               !! screen report timing (1/cycle)")
        println(io, "  strict_mode      = $( _b(cfg.strict_mode) )          !! all parameters to be explicitly definied\n")

        println(io, "  !! ----------------------------------------------------------------------- !!\n  !! Model/Grid Size and Area\n  !!\n")
        println(io, "  nproc_x          = $(cfg.nproc_x)                !! parallelization in x-dir")
        println(io, "  nproc_y          = $(cfg.nproc_y)                !! parallelization in y-dir")
        println(io, "  nx               = $(cfg.nx)              !! total grid number in x-dir")
        println(io, "  ny               = $(cfg.ny)              !! total grid number in y-dir")
        println(io, "  nz               = $(cfg.nz)              !! total grid number in z-dir")
        println(io, "  nt               = $(cfg.nt)             !! time step number\n")
        println(io, "  dx               = $(cfg.dx)              !! grid width in x-dir")
        println(io, "  dy               = $(cfg.dy)              !! grid width in y-dir")
        println(io, "  dz               = $(cfg.dz)              !! grid width in z-dir")
        println(io, "  dt               = $(cfg.dt)             !! time step width\n")
        println(io, "  vcut             = $(cfg.vcut)              !! minimum velocity\n                                      !- smaller velocities will be increased\n")
        println(io, "  xbeg             = $(cfg.xbeg)            !! minimum in x-dir")
        println(io, "  ybeg             = $(cfg.ybeg)            !! minimum in y-dir")
        println(io, "  zbeg             = $(cfg.zbeg)            !! minimum in z-dir")
        println(io, "  tbeg             = $(cfg.tbeg)              !! start time\n")
        println(io, "  clon             = $(cfg.clon)         !! center longitude")
        println(io, "  clat             = $(cfg.clat)          !! center latitude")
        println(io, "  phi              = $(cfg.phi)              !! horizontal coordinate rotation\n                                      !- measured clockwise from the north\n")
        println(io, "  fq_min           = $(cfg.fq_min)             !! minimum freq. for Q-constant model")
        println(io, "  fq_max           = $(cfg.fq_max)             !! maximum freq. for Q-constant model")
        println(io, "  fq_ref           = $(cfg.fq_ref)              !! ref. freq. for physical dispersion\n")

        println(io, "  !! ----------------------------------------------------------------------- !!\n  !! Snapshot Output\n  !!\n")
        println(io, "  snp_format       = $( _qs(cfg.snp_format) )         !! snapshot format (netcdf)\n")
        println(io, "  xy_ps%sw         = $( _b(cfg.xy_ps_sw) )          !! P&S amp. for xy section")
        println(io, "  xz_ps%sw         = $( _b(cfg.xz_ps_sw) )           !! P&S amp. for xz section")
        println(io, "  yz_ps%sw         = $( _b(cfg.yz_ps_sw) )          !! P&S amp. for yz section")
        println(io, "  fs_ps%sw         = $( _b(cfg.fs_ps_sw) )          !! P&S amp. for free surface")
        println(io, "  ob_ps%sw         = $( _b(cfg.ob_ps_sw) )           !! P&S amp. for ocean bottom\n")
        println(io, "  xy_v%sw          = $( _b(cfg.xy_v_sw) )          !! 3-comp. velocity for xy section")
        println(io, "  xz_v%sw          = $( _b(cfg.xz_v_sw) )           !! 3-comp. velocity for xz section")
        println(io, "  yz_v%sw          = $( _b(cfg.yz_v_sw) )          !! 3-comp. velocity for yz section")
        println(io, "  fs_v%sw          = $( _b(cfg.fs_v_sw) )          !! 3-comp. velocity for free surface")
        println(io, "  ob_v%sw          = $( _b(cfg.ob_v_sw) )           !! 3-comp. velocity for ocean bottom\n")
        println(io, "  xy_u%sw          = $( _b(cfg.xy_u_sw) )          !! 3-comp. disp. for xy section")
        println(io, "  xz_u%sw          = $( _b(cfg.xz_u_sw) )           !! 3-comp. disp. for xz section")
        println(io, "  yz_u%sw          = $( _b(cfg.yz_u_sw) )          !! 3-comp. disp. for yz section")
        println(io, "  fs_u%sw          = $( _b(cfg.fs_u_sw) )          !! 3-comp. disp. for free surface")
        println(io, "  ob_u%sw          = $( _b(cfg.ob_u_sw) )           !! 3-comp. disp. for ocean bottom\n")
        println(io, "  vol_v%sw         = $( _b(cfg.vol_v_sw) )          !! 3-comp. velocity for 3D volume")
        println(io, "  vol_u%sw         = $( _b(cfg.vol_u_sw) )          !! 3-comp. disp. for 3D volume")
        println(io, "  vol_ps%sw        = $( _b(cfg.vol_ps_sw) )          !! P&S amp. for 3D volume\n")

        println(io, "  z0_xy            =  $(cfg.z0_xy)             !! depth for xy cross section")
        println(io, "  x0_yz            =  $(cfg.x0_yz)             !! x-value for yz cross section")
        println(io, "  y0_xz            =  $(cfg.y0_xz)             !! y-value for xz cross section\n")
        println(io, "  ntdec_s          = $(cfg.ntdec_s)                !! time decimation of snapshot\n                                      !- (specify 1 for no decimation)\n")
        println(io, "  idec             = $(cfg.idec)                !! x-decimation for snapshot")
        println(io, "  jdec             = $(cfg.jdec)                !! y-decimation for snapshot")
        println(io, "  kdec             = $(cfg.kdec)                !! z-decimation for snapshot\n")

        println(io, "  !! ----------------------------------------------------------------------- !!\n  !! Waveform Output\n  !!\n")
        println(io, "  sw_wav_v         = $( _b(cfg.sw_wav_v) )           !! velocity trace output at stations")
        println(io, "  sw_wav_u         = $( _b(cfg.sw_wav_u) )           !! displacement trace output at stations")
        println(io, "  sw_wav_stress    = $( _b(cfg.sw_wav_stress) )           !! stress tensor trace")
        println(io, "  sw_wav_strain    = $( _b(cfg.sw_wav_strain) )           !! strain tansor trace")
        println(io, "  ntdec_w          = $(cfg.ntdec_w)                !! time decimation of waveform output")
        println(io, "  st_format        = $( _qs(cfg.st_format) )             !! station format: 'xy' or 'll'")
        println(io, "  fn_stloc         = $( _qs(cfg.fn_stloc) )  !! station location file")
        println(io, "  wav_format       = $( _qs(cfg.wav_format) )            !! 'sac' or 'csf' ('sac' recommended)")
        println(io, "  ntdec_w_prg      = $(cfg.ntdec_w_prg)                !!  waveform output during computation (0:off)\n")

        println(io, "  !! ----------------------------------------------------------------------- !!\n  !! Earthquake Source\n  !!\n")
        println(io, "  stf_format       = $( _qs(cfg.stf_format) )\n")
        println(io, "  stftype          = $( _qs(cfg.stftype) )\n")
        println(io, "  fn_stf           = $( _qs(cfg.fn_stf) )   !! Source grid file name\n")
        println(io, "  sdep_fit         = $( _qs(cfg.sdep_fit) )\n")

        println(io, "  !! --------------------------------------------------------------------- !!\n    !! Body force source mode\n    !!\n    bf_mode          = $( _b(cfg.bf_mode) )\n\n")

        println(io, "  !! --------------------------------------------------------------------- !!\n    !! Plane wave source mode\n    !!\n    pw_mode          = $( _b(cfg.pw_mode) )   !! plane wave input. neglects fn_stf")
        println(io, "  pw_ztop          = $(cfg.pw_ztop)      !! top z-coordinate of the initial plane wave")
        println(io, "  pw_zlen          = $(cfg.pw_zlen)       !! wavelength of the initial plane wave")
        println(io, "  pw_ps            = $( _qs(cfg.pw_ps) )       !! 'p' P-wave 's' S-wave")
        println(io, "  pw_strike        = $(cfg.pw_strike)       !! strike direction of plane wave (deg.)")
        println(io, "  pw_dip           = $(cfg.pw_dip)       !! dip of plane wave (deg.)")
        println(io, "  pw_rake          = $(cfg.pw_rake)       !! rake of plane S-wave polarization (deg.)\n")

        println(io, "  !! ----------------------------------------------------------------------- !!\n  !! Absorbing Boundary Condition\n  !!\n")
        println(io, "  abc_type         = $( _qs(cfg.abc_type) )            !! 'pml' or 'cerjan'")
        println(io, "  na               = $(cfg.na)               !! absorbing layer thickness")
        println(io, "  stabilize_pml    = $( _b(cfg.stabilize_pml) )           !! avoid low-v layer in PML region\n")

        println(io, "  !! ----------------------------------------------------------------------- !!\n  !! Velocity model\n  !!\n")
        println(io, "  vmodel_type      = $( _qs(cfg.vmodel_type) )             !! velocity model type 'uni'/'grd'/'lhm'")
        println(io, "  is_ocean         = $( _b(cfg.is_ocean) )                 !! topography z<0 is covered by ocean")
        println(io, "  topo_flatten     = $( _b(cfg.topo_flatten) )             !! Force topography variation to zero (formerly is_flatten)")
        println(io, "  munk_profile     = $( _b(cfg.munk_profile) )             !! velocity gradient inside the seawater column")
        println(io, "  earth_flattening = $( _b(cfg.earth_flattening) )         !! Earth-flattening tranformation\n")

        if cfg.vmodel_type == "uni"
        println(io, "   !! --------------------------------------------------------------------- !!\n    !! For uniform velocity model 'uni'\n    !!\n")
        println(io, "   vp0              = $(cfg.vp0)              !! P-wave velocity [km/s]")
        println(io, "   vs0              = $(cfg.vs0)              !! S-wave velocity [km/s]")
        println(io, "   rho0             = $(cfg.rho0)              !! mass density    [g/cm^3]")
        println(io, "   qp0              = $(cfg.qp0)              !! Qp")
        println(io, "   qs0              = $(cfg.qs0)              !! Qs")
        println(io, "   topo0            = $(cfg.topo0)                !! topography location\n")
        end
        if cfg.vmodel_type == "grd"
        println(io, "   !! --------------------------------------------------------------------- !!\n    !! For GMT grid file input 'grd' ( requires netcdf library )\n    !!\n")
        println(io, "   dir_grd          = $( _qs(cfg.dir_grd) )    !! directory for grd file")
        println(io, "   fn_grdlst        = $( _qs(cfg.fn_grdlst) )            !! grd file list")
        println(io, "   node_grd         = $(cfg.node_grd)                              !! input MPI node\n")
        end
        if cfg.vmodel_type == "lhm" || cfg.vmodel_type == "lgm"
        println(io, "   !! --------------------------------------------------------------------- !!\n    !! For layered homogeneous medium model ('lhm')\n    !!\n")
        println(io, "   fn_lhm           = $( _qs(cfg.fn_lhm) )    !! 1D velocity structure\n")
        end
        if cfg.vmodel_type == "velm" 
        println(io, "   !! --------------------------------------------------------------------- !!\n    !! For velocity model from file ('velm')\n    !!\n")
        println(io, "   topo0            = $(cfg.topo0)                !! topography location")
        println(io, "   dir_velm         = $( _qs(cfg.dir_velm) )    !! directory for velocity model file")
        println(io, "   fn_velm          = $( _qs(cfg.fn_velm) )     !! velocity model file\n")
        end

        println(io, "   !! --------------------------------------------------------------------- !!\n    !! For random medium models\n    !!\n")
        println(io, "   dir_rmed         = $( _qs(cfg.dir_rmed) )             !! location of random medium file")
        println(io, "   fn_grdlst_rmed   = $( _qs(cfg.fn_grdlst_rmed) ) !! grd file list")
        println(io, "   rhomin           = $(cfg.rhomin)                 !! minimum density threshold")
        println(io, "   fn_rmed0         = $( _qs(cfg.fn_rmed0) )          !! vel. purturb. on a uniform media\n")

        println(io, "  !! ----------------------------------------------------------------------- !!\n  !! Checkpoint/Restart\n  !!\n")
        println(io, "  is_ckp           = $( _b(cfg.is_ckp) )          !! perform checkpoint/restart")
        println(io, "  ckpdir           = $( _qs(cfg.ckpdir) )      !! output directory")
        println(io, "  ckp_interval     = $(cfg.ckp_interval)          !! interval for checkpoint check（1/cycle）")
        println(io, "  ckp_time         = $(cfg.ckp_time)          !! checkpoint time")
        println(io, "  ckp_seq          = $( _b(cfg.ckp_seq) )              !! sequential output mode\n")

        println(io, "  !! ----------------------------------------------------------------------- !!\n  !! Reciprocity Green's Function Mode\n  !!\n")
        println(io, "  green_mode       = $( _b(cfg.green_mode) )          !! reciprocity Green's function mode")
        println(io, "  green_stnm       = $( _qs(cfg.green_stnm) )           !! virtual station name from fn_stlst")
        println(io, "  green_cmp        = $( _qs(cfg.green_cmp) )              !! virtual source direction 'x', 'y', 'z'")
        println(io, "  green_trise      = $(cfg.green_trise)              !! rise time")
        println(io, "  green_bforce     = $( _b(cfg.green_bforce) )          !! also calc. body force Green's function")
        println(io, "  green_maxdist    = $(cfg.green_maxdist)              !! horizontal limit of source grid")
        println(io, "  green_fmt        = $( _qs(cfg.green_fmt) )            !! list file format: 'xyz' or 'llz'")
        println(io, "  fn_glst          = $( _qs(cfg.fn_glst) )   !! Green's function grid point list\n")

        println(io, "  !! ----------------------------------------------------------------------- !!\n  !! MISC\n  !!\n")
        println(io, "  stopwatch_mode   = $( _b(cfg.stopwatch_mode) )          !! measure computation time at routines")
        println(io, "  benchmark_mode   = $( _b(cfg.benchmark_mode) )          !! benchmark mode\n")
        println(io, "  ipad             = $(cfg.ipad)                !! memory padding size for tuning")
        println(io, "  jpad             = $(cfg.jpad)                !! memory padding size for tuning")
        println(io, "  kpad             = $(cfg.kpad)                !! memory padding size for tuning")
    end

    # write sources
    if !isempty(cfg.source)
        write_sourceCF!(cfg.fn_stf, cfg.source)
    end

    # write layered models if needed
    if cfg.vmodel_type == "lhm" || cfg.vmodel_type == "lgm"
        write_lhm!(cfg.fn_lhm, cfg.vmodel)
    end

    # write velocity model if needed 
    if cfg.vmodel_type == "velm"
        write_netcdf(cfg.vmodel.data, joinpath(cfg.dir_velm, cfg.fn_velm))
    end

    if !isempty(cfg.stations)
        write_stations_ll!(cfg.fn_stloc, cfg.stations)
    end


    return nothing
end



function clean(cfg::OpenSWPCConfig)
    rm(cfg.odir, recursive=true)
    rm(cfg.input_file, force=true)
    rm(cfg.fn_stf, force=true)
    rm(cfg.fn_lhm, force=true)
    rm(cfg.fn_stloc, force=true)   
    return nothing
end




"""
    OpenSWPCConfig(input_model::CartData; kwargs...)

Generates an input model from a GMG CartData object.
"""
function OpenSWPCConfig(input_model::CartData; kwargs...)
    nx, ny, nz = size(input_model)
    dx = input_model.x.val[2,2,2]-input_model.x.val[1,1,1]
    dy = input_model.y.val[2,2,2]-input_model.y.val[1,1,1]
    dz = input_model.z.val[2,2,2]-input_model.z.val[1,1,1]
    xbeg = minimum(input_model.x.val)   
    ybeg = minimum(input_model.y.val)
    zbeg = -maximum(input_model.z.val)      # coordinates are flipped
    vmodel = FullVelocityModel(input_model)
    @assert haskey(input_model.fields, :Qp) "CartData must have field 'Qp'"
    @assert haskey(input_model.fields, :Qs) "CartData must have field 'Qs'"
    @assert haskey(input_model.fields, :rho) "CartData must have field 'rho'"
    @assert haskey(input_model.fields, :mu) "CartData must have field 'mu'"
    @assert haskey(input_model.fields, :lambda) "CartData must have field 'lambda'"

    return OpenSWPCConfig(;
        nx=nx, ny=ny, nz=nz,
        dx=dx, dy=dy, dz=dz,
        xbeg=xbeg, ybeg=ybeg, zbeg=zbeg,
        vmodel=vmodel,
        kwargs...
    )
end