# -*- coding: utf-8 -*-
"""WP1d -- MMsurf, driven by the configuration instead of two ini files.

MMsurf turns the hourly meteorological record into the daily forcing MMsoil
consumes: rainfall, throughfall, potential transpiration per vegetation type,
potential evaporation per surface soil, open-water evaporation and LAI. It is
panel 2 of the front-end, and it runs when ``run.surface`` is true.

What this module removes
------------------------
``__inputMMsurf4MMsoil.txt``. That file was written by MMsurf and read back by
the driver, positionally, and it was AUTHORITATIVE: the driver took `Zr`,
`kTg_min`, `kTg_max`, `kT_f` and `kT_s` from it, not from the ini, so editing
them in the ini changed nothing while MMsurf was not running. Worse, the two
disagreed -- the committed handover says ``kTg_min = 0.044694325`` where the
committed ini says ``0.04440913``.

Now there is one place: the configuration. ``forcing_spec(cfg)`` returns what
the driver used to parse out of the handover, and ``run(cfg, ...)`` runs
MMsurf from the same values. Nothing is written to be read back.

The ini is still WRITTEN, into the workspace, because MMsurf's own parser
reads it -- 400 lines of positional `l += 1` that would have to be rewritten
to no benefit. It is a generated artefact of the run, not an input, and it
carries the configuration hash so it can be traced back.

A note on kt_s
--------------
The ini's nineteenth vegetation field is commented ``kT_1/s`` and holds values
above 20; the driver applied ``1/x`` before using it. So the FILE stored 1/s
while the model wanted s. The configuration holds **s**, and the inversion
happens here, on the way out -- the one place that has to know.
"""

import os


__all__ = ['MMsurfError', 'FORCING', 'forcing_spec', 'write_par_file', 'run',
           'check_forcing', 'surface_ws']


class MMsurfError(Exception):
    """MMsurf cannot be run, or its forcing cannot be used."""


# The daily forcing files, by role. Canonical names: the ones the committed
# dataset and the driver use. (The module in the repo wrote inputZONP_* for
# the first two, which is the same quantity under an older spelling and the
# reason it could not reproduce its own output -- WP1d finding B.1.)
FORCING = {
    'date': 'inputDATE.txt',
    'rf_veg': 'inputZONRF_veg_d.txt',      # rainfall
    'tf_veg': 'inputZONTF_veg_d.txt',      # throughfall -- the driver's "Pe"
    'pt_veg': 'inputZONPT_veg_d.txt',      # potential transpiration
    'lai_veg': 'inputZONLAI_veg_d.txt',
    'pe': 'inputZONPE_d.txt',              # potential soil evaporation
    'eo': 'inputZONEo_d.txt',              # open-water evaporation
    'rf_irr': 'inputZONRF_irr_d.txt',
    'tf_irr': 'inputZONTF_irr_d.txt',
    'pt_irr': 'inputZONPT_irr_d.txt',
    'crop_irr': 'inputZONcrop_irr_d.txt',
}

VEG_ORDER = ('rf_veg', 'tf_veg', 'pt_veg', 'lai_veg', 'pe', 'eo')
IRR_ORDER = ('rf_irr', 'tf_irr', 'pt_irr')

# Each forcing file is a flat column of values, blocks concatenated. How many
# blocks depends on what the quantity varies with -- and nothing in the file
# says so, which is why it is written down here.
#   per meteo zone x vegetation type : rf? no -- rainfall is per zone only
_BLOCKS = {
    'rf_veg': ('meteo',),
    'tf_veg': ('meteo', 'veg'),
    'pt_veg': ('meteo', 'veg'),
    'lai_veg': ('veg',),            # LAI does not depend on the meteo zone
    'pe': ('meteo', 'soil'),
    'eo': ('meteo',),
    'rf_irr': ('meteo', 'field'),
    'tf_irr': ('meteo', 'field'),
    'pt_irr': ('meteo', 'field'),
    'crop_irr': ('field',),
    'date': (),
}


class ForcingSpec(object):
    """What the driver used to read out of the handover file.

    Everything here comes from the configuration. ``kt_s`` is already the
    slope the model wants, so the driver's ``1.0 / float(x)`` is gone.
    """

    def __init__(self, cfg, ws):
        s = cfg.surface
        self.ws = ws
        self.nmeteo = max(1, len(s.station))
        self.nveg = len(s.vegetation)
        self.nsoil = len(s.soil)
        self.ncrop = len(s.crop)
        self.nfield = int(s.nfield) if s.irrigation else 0
        self.irrigation = bool(s.irrigation)
        self.veg_names = [v.name for v in s.vegetation]
        self.soil_names = [x.name for x in s.soil]
        self.Zr = [v.root_depth for v in s.vegetation]
        self.kTg_min = [v.ktg_min for v in s.vegetation]
        self.kTg_max = [v.ktg_max for v in s.vegetation]
        self.kT_f = [v.kt_f for v in s.vegetation]
        self.kT_s = [v.kt_s for v in s.vegetation]
        self.Zr_c = [c.root_depth for c in s.crop]
        self.kTg_min_c = [c.ktg_min for c in s.crop]
        self.kTg_max_c = [c.ktg_max for c in s.crop]
        self.kT_f_c = [c.kt_f for c in s.crop]
        self.kT_s_c = [c.kt_s for c in s.crop]
        self.files = dict(FORCING)

    def path(self, role):
        return os.path.join(self.ws, self.files[role])

    def blocks(self, role):
        """How many concatenated blocks this forcing file should hold."""
        n = 1
        for what in _BLOCKS[role]:
            n *= {'meteo': self.nmeteo, 'veg': self.nveg,
                  'soil': self.nsoil, 'field': max(self.nfield, 1)}[what]
        return n

    def roles(self):
        out = ['date'] + list(VEG_ORDER)
        if self.irrigation:
            out += list(IRR_ORDER) + ['crop_irr']
        return out

    def __repr__(self):
        return ('<ForcingSpec NMETEO=%d NVEG=%d NSOIL=%d NCROP=%d NFIELD=%d>'
                % (self.nmeteo, self.nveg, self.nsoil, self.ncrop, self.nfield))


def forcing_spec(cfg, ws):
    """Replaces the positional parsing of ``__inputMMsurf4MMsoil.txt``."""
    return ForcingSpec(cfg, ws)


def surface_ws(cfg, ws_root, case):
    """Where MMsurf writes. Run OUTPUT, so never the repository."""
    return os.path.join(str(ws_root), '%s_MMsurf' % case)


# =====================================================================
#  the parameter file MMsurf's own parser reads
# =====================================================================

def write_par_file(cfg, path, config_hash=''):
    """Render the configuration into the ini layout MMsurf parses.

    Written into the WORKSPACE, never the dataset: it is a generated
    artefact, and it carries the configuration hash so a run can be traced
    back to the settings that produced it.
    """
    s = cfg.surface
    L = ['# GENERATED by code/marmites_surface.py -- do not edit.',
         '# This file exists only because MMsurf parses it; the values live in',
         '# the run configuration. config_hash = %s' % (config_hash or 'n/a'),
         '#',
         '# METEO PARAMETERS',
         '# NMETEO', '%d' % max(1, len(s.station)),
         '# phi Lm Z Lz FC DTS z_m z_h   (phi/Lm derived from x/y, project CRS)']
    for st in (s.station or [None]):
        if st is None:
            L.append('0.0 0.0 0.0 0.0 0.0 0.0 2.0 2.0')
            continue
        phi, lm = geographic(cfg, st)
        L.append('%.6f %.6f %.3f %.3f %.3f %.3f %.3f %.3f'
                 % (phi, lm, st.z, st.tz_lon, st.fuse_shift, st.data_shift,
                    st.z_wind, st.z_hum))
    L += ['', '#VEGETATION PARAMETERS', '# NVEG', '%d' % len(s.vegetation)]
    for v in s.vegetation:
        # The file's last field is 1/s; the configuration holds s.
        L.append('%s %g %g %g %g %g %g %g %g %g %g %d %d %d %d %g %g %g %g %g'
                 % (v.name, v.h_dry, v.h_wet, v.canopy_mm, v.c_leaf,
                    v.lai_dry, v.lai_wet, v.shelter_dry, v.shelter_wet,
                    v.albedo_dry, v.albedo_wet, v.j_dry, v.j_wet,
                    v.trans_dry_wet, v.trans_wet_dry, v.root_depth,
                    v.ktg_min, v.ktg_max, v.kt_f, 1.0 / v.kt_s))
    L += ['', '#CROP PARAMETERS', '# NCRP', '%d' % len(s.crop),
          '# NFIELD', '%d' % (int(s.nfield) if s.irrigation else 0)]
    for c in s.crop:
        L.append('%s %g %g %g %g %g %g %g %g %g %g %g'
                 % (c.name, c.h, c.canopy_mm, c.c_leaf, c.lai, c.shelter,
                    c.albedo, c.root_depth, c.ktg_min, c.ktg_max, c.kt_f,
                    1.0 / c.kt_s))
    L += ['', '# SOIL PARAMETERS', '# NSOIL', '%d' % len(s.soil)]
    for x in s.soil:
        L.append('%s %g %g %g %g %d %d %d %d'
                 % (x.name, x.porosity, x.field_capacity, x.albedo_dry,
                    x.albedo_wet, x.j_dry, x.j_wet, x.trans_dry_wet,
                    x.trans_wet_dry))
    L.append('')
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(L))
    return path


def geographic(cfg, station):
    """Station latitude and longitude-WEST, from its project-CRS position.

    Penman-Monteith needs the geographic pair for its solar geometry; the
    modeller enters metres, because that is what every other coordinate in
    the project is. The old ini held the geographic pair directly and on La
    Mata it was 6.8 km out -- which never showed, because with one station
    the Thiessen polygon is the whole catchment either way.
    """
    epsg = int(getattr(cfg.grid, 'crs_epsg', 0) or 0)
    if not epsg:
        raise MMsurfError(
            'grid.crs_epsg is not set, so the station coordinates cannot be '
            'converted to the latitude/longitude Penman-Monteith needs. Set '
            'it to the project CRS (La Mata: 23029).')
    try:
        from pyproj import CRS, Transformer
    except ImportError as exc:                               # pragma: no cover
        raise MMsurfError('pyproj is needed to convert the station position '
                          'from EPSG:%d to geographic' % epsg) from exc
    t = Transformer.from_crs(CRS.from_epsg(epsg), CRS.from_epsg(4326),
                             always_xy=True)
    lon, lat = t.transform(float(station.x), float(station.y))
    return lat, -lon                       # MMsurf wants longitude WEST


# =====================================================================
#  running it
# =====================================================================

def run(cfg, dataset_dir, out_ws, config_hash='', verbose=True):
    """Run MMsurf and return the ForcingSpec pointing at what it wrote.

    Inputs (the meteo record, the irrigation series, the crop schedules) are
    read from ``<dataset_dir>/MMsurf_ws``; every output goes to ``out_ws``.
    """
    import sys

    here = os.path.dirname(os.path.abspath(__file__))
    # startMARMITESsurface imports plotPET and plotP as top-level modules, but
    # they live under MARMITESutilities/MARMITESplot -- so that folder has to
    # be on the path too, or the import fails with a name that says nothing
    # about where it came from.
    for p in (os.path.join(here, 'MARMITESsurf'),
              os.path.join(here, 'MARMITESutilities'),
              os.path.join(here, 'MARMITESutilities', 'MARMITESplot'),
              here):
        if p not in sys.path:
            sys.path.insert(0, p)
    import MARMITESutilities as MMutils
    import startMARMITESsurface as MMsurf_mod

    in_ws = os.path.join(dataset_dir, 'MMsurf_ws')
    if not os.path.isdir(in_ws):
        raise MMsurfError('MMsurf input folder not found: %s' % in_ws)
    os.makedirs(out_ws, exist_ok=True)

    s = cfg.surface
    meteo = os.path.join(in_ws, s.meteo_ts)
    if not os.path.exists(meteo):
        raise MMsurfError('the meteorological record is missing: %s' % meteo)

    # The generated parameter file lives in the WORKSPACE and is passed by
    # ABSOLUTE path: MMsurf resolves it with os.path.join(in_ws, name), which
    # returns an absolute second argument unchanged. Copying it in beside the
    # inputs would work too, and would leave a stray file in the repository
    # whenever a run crashed.
    par_name = os.path.join(out_ws, '__inputMMsurf.generated.ini')
    write_par_file(cfg, par_name, config_hash)

    irr_ts = s.irr_ts if s.irrigation else None
    if s.irrigation:
        if not os.path.exists(os.path.join(in_ws, s.irr_ts)):
            raise MMsurfError('irrigation is on but %s is missing from %s'
                              % (s.irr_ts, in_ws))
        for f in range(int(s.nfield)):
            sched = os.path.join(in_ws, s.crop_schedule % (f + 1))
            if not os.path.exists(sched):
                raise MMsurfError('irrigation is on but the crop schedule of '
                                  'field %d is missing: %s' % (f + 1, sched))

    if verbose:
        print('MMsurf: %d station(s), %d vegetation type(s), %d soil(s)%s'
              % (max(1, len(s.station)), len(s.vegetation), len(s.soil),
                 ', %d field(s)' % s.nfield if s.irrigation else ''))
        print('   input  : %s' % in_ws)
        print('   output : %s' % out_ws)

    cUTIL = MMutils.clsUTILITIES(verbose=1 if verbose else 0)
    MMsurf_mod.MMsurf(
        cUTIL, in_ws, s.meteo_ts, par_name, s.out_prefix,
        out_ws, '__inputMMsurf4MMsoil.generated.txt',
        MMsurf_plot=int(s.plot), inputFile_IRR_TS_fn=irr_ts,
        out_ws=out_ws)

    spec = forcing_spec(cfg, out_ws)
    check_forcing(spec, must_exist=True)
    return spec


def check_forcing(spec, must_exist=True, nper=None):
    """Are the forcing series there, and do they have the right shape?

    ``run.surface = 0`` means the series must already exist. Running on a
    stale or truncated file is the failure this prevents: the committed
    ``inputZONRFe_veg_d.txt`` holds 4869 values against a 1949-day record,
    which is not a whole number of blocks and went unnoticed for years.
    """
    missing, bad = [], []
    ndays = None
    for role in spec.roles():
        p = spec.path(role)
        if not os.path.exists(p):
            missing.append(os.path.basename(p))
            continue
        # One record per LINE, after a '#' header line. Counting whitespace
        # tokens instead would work for the flux files -- one number a line --
        # and silently treble the count for inputDATE.txt, whose lines are
        # "YYYY-MM-DD HH:MM, J".
        with open(p, encoding='utf-8', errors='replace') as fh:
            lines = [ln for ln in fh.read().splitlines() if ln.strip()]
        n = len(lines) - 1 if lines and lines[0].lstrip().startswith('#') \
            else len(lines)
        if role == 'date':
            ndays = n
            continue
        blocks = spec.blocks(role)
        if ndays and n != ndays * blocks:
            bad.append('%s holds %d values, expected %d (%d day(s) x %d block(s))'
                       % (os.path.basename(p), n, ndays * blocks, ndays, blocks))
    if missing and must_exist:
        raise MMsurfError(
            'run.surface is off, so the forcing series must already exist, '
            'but these are missing from %s:\n  - %s\nSet run.surface = true '
            'to generate them.' % (spec.ws, '\n  - '.join(missing)))
    if bad:
        raise MMsurfError('the forcing series do not have the expected shape:'
                          '\n  - ' + '\n  - '.join(bad))
    if nper and ndays is not None and ndays < nper:
        raise MMsurfError('the forcing covers %d day(s) but the run asks for '
                          '%d stress period(s)' % (ndays, nper))
    return ndays
