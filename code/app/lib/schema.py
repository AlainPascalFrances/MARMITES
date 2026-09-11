# -*- coding: utf-8 -*-
"""WP1d -- what the four panels show, and what every field means.

DELIBERATELY STREAMLIT-FREE, like the rest of ``code/app/lib``: the panels are
thin views over this, and what could actually be wrong -- a field in no panel,
a panel naming a field that does not exist, a label without units -- is tested
here rather than discovered in a browser.

The panel order is the modeller's:

    0  Overview        what this is, and what to fill in first
    1  Grid            the catchment polygon and the grid built inside it
    2  Surface         MMsurf: the meteorological record -> the daily forcing
    3  Model           MMsoil and MODFLOW 6, which run together
    4  Plots           the figures

Panels 2, 3 and 4 carry a MASTER SWITCH, and it is not decoration: the same
``[run]`` key the driver reads decides whether that half of the model runs.

Labels come from the ini files the panels replace, so a modeller who knows
``kTg_min`` can still find it, with a sentence saying what it does.
"""

import os
import sys

# The allowed values of an enumerated field come from marmites_config, so this
# module has to be able to import it however it was loaded -- as part of the
# app, or standalone by a test.
_CODE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _CODE not in sys.path:
    sys.path.insert(0, _CODE)

__all__ = ['PANELS', 'FIELDS', 'TABLES', 'CHOICES', 'SUBPANELS', 'panel_of',
           'describe', 'fields_of', 'choices_for', 'subpanel_for', 'is_source',
           'PanelError']


class PanelError(Exception):
    """A panel refers to something the configuration does not have."""


# ---------------------------------------------------------------- panels
# (number, title, icon, master switch or None, [section, ...], blurb)
PANELS = [
    (0, 'Overview', '💧', None, [],
     'What this case is, where its files are, and the order to fill the '
     'panels in.'),
    (1, 'Grid', '🗺️', None, ['grid'],
     'The catchment polygon comes FIRST and the grid is built inside it. '
     'Every other input is then wrapped onto whatever this produces, which '
     'is why this panel is answered before any of the others.'),
    (2, 'Surface', '🌦️', 'run.surface', ['surface'],
     'MMsurf turns the hourly meteorological record into the daily forcing '
     'MMsoil consumes. With the switch off, that forcing must already exist '
     'and is checked for presence and shape before the run starts.'),
    (3, 'Model', '🌍', 'run.model', ['soil', 'obs', 'layers', 'uzf', 'seep',
                                     'et', 'sfr', 'lak', 'crr', 'spinup'],
     'The soil water balance and MODFLOW 6. One switch, because they are no '
     'longer separable: the legacy MMsoil-only mode belonged to the Picard '
     'loop that Phase 1 removed.'),
    (4, 'Plots', '📊', 'run.plot', ['postproc'],
     'What to draw once the run finishes. Nothing here changes a flux.'),
]

# Arrays of tables: one row per zone / type, edited as a grid rather than as
# a wall of positional numbers. (dotted path, panel, singular label, the
# count it sets)
TABLES = [
    ('surface.station', 2, 'Meteorological station', 'NMETEO'),
    ('surface.vegetation', 2, 'Vegetation type', 'NVEG'),
    ('surface.crop', 2, 'Crop', 'NCRP'),
    ('surface.soil', 2, 'Surface soil', 'NSOIL'),
    ('soil.veg_class', 3, 'Vegetation class mapping', None),
]

# ---------------------------------------------------------------- fields
# dotted key -> (label, units, what it does)
_U = ''            # dimensionless / not applicable

FIELDS = {
    # ---- meta / paths -------------------------------------------------
    'meta.name': ('Run name', _U,
                  'Names the output folder, out_<timestamp>_<name>.'),
    'meta.description': ('Description', _U, 'Free text, carried into the run.'),
    'paths.case': ('Case', _U, 'Resolves to example/<case>/, the only place '
                               'the model reads input from.'),
    'paths.ws': ('Workspace', _U,
                 'Blank uses WS_ROOT from mm_paths. Everything a run produces '
                 'goes here -- never into the repository.'),
    'paths.libmf6': ('MODFLOW 6 library', _U,
                     'Blank uses mm_paths.LIBMF6; "auto" searches for it. '
                     'Without it the run stops after writing the input files.'),

    # ---- run ----------------------------------------------------------
    'run.surface': ('Run MMsurf', _U,
                    'Off means the daily forcing must already exist.'),
    'run.model': ('Run MMsoil + MODFLOW 6', _U,
                  'They run together; there is no longer a way to run one '
                  'without the other.'),
    'run.plot': ('Draw the figures', _U, 'Post-processing after the run.'),
    'run.mode': ('Coupling mode', _U,
                 'lagged: MMsoil once per stress period, from the previous '
                 'period\'s heads. iterative: re-evaluated at every MODFLOW '
                 'outer iteration, which removes the one-period lag.'),
    'run.relax': ('Under-relaxation', _U,
                  'Iterative mode only: new = relax*evaluated + '
                  '(1-relax)*previous. 0.5-0.7 is the usable range.'),
    'run.nsp': ('Stress periods', 'count',
                '0 runs the whole record. Use a small number to try a change '
                'before committing to the full run.'),
    'run.daily': ('One stress period per day', _U,
                  'Off aggregates the record by rainfall episode.'),
    'run.ats': ('Adaptive time stepping', _U,
                'On. A period MODFLOW cannot solve in one step is not a '
                'result, and without this it silently becomes one.'),
    'run.build_only': ('Write the input files and stop', _U,
                       'Useful to inspect what MODFLOW would be given.'),
    'run.max_discrepancy': ('Maximum mass-balance discrepancy', '%',
                            'The run fails above this cumulative value.'),

    # ---- panel 1: the grid --------------------------------------------
    'grid.boundary': ('Catchment polygon', _U,
                      'A PROJECTED, metric shapefile in DATA_ROOT/GIS. It '
                      'defines the active domain, the mesh boundary and the '
                      'model rectangle -- so it is the first thing to set.'),
    'grid.crs_epsg': ('Project CRS', 'EPSG',
                      'The projected, metric CRS everything is in. 0 takes it '
                      'from the layer\'s .prj. Station coordinates are '
                      'converted from here to the latitude and longitude '
                      'Penman-Monteith needs.'),
    'grid.cell_size': ('Cell size', 'm',
                       'Background cell size, for every grid kind.'),
    'grid.buffer': ('Buffer', 'm',
                    'Extend the model rectangle beyond the polygon.'),
    'grid.kind': ('Grid kind', _U,
                  'structured is the validated DIS grid and the regression '
                  'anchor. voronoi and quadtree are real meshes; disv is the '
                  'structured grid re-expressed as polygons, used to check '
                  'that the mesh path changes nothing by itself.'),
    'grid.rebuild': ('Rebuild the mesh', _U,
                     'The cache carries the mesh signature, so this is only '
                     'needed to force a rebuild.'),
    'grid.resample': ('Resampling rule', _U,
                      'auto area-weights continuous fields and majority-votes '
                      'zone maps. centre samples the cell centre -- worse per '
                      'cell, but it preserves zone proportions.'),
    'grid.voronoi.cell_far': ('Background cell size', 'm',
                              'The side of the EQUIVALENT SQUARE, so 100 aims '
                              'at 10 000 m2. Every build prints the area it '
                              'actually achieved -- read that, not this.'),
    'grid.voronoi.cell_near_stream': ('Cell size near the stream', 'm', ''),
    'grid.voronoi.stream_buffer': ('Stream corridor width', 'm', ''),
    'grid.voronoi.stream_refine': ('Refine along the streams', _U,
                                   'Off is the SFRmaker approach: map the '
                                   'network onto the background grid. CdL '
                                   'settled on off after judging a refined '
                                   'mesh too refined near streams.'),
    'grid.voronoi.seed_ponds': ('Seed a cell per pond', _U, ''),
    'grid.quadtree.refine_level': ('Refinement levels', 'count',
                                   'GRIDGEN halves a cell per level, so 2 on '
                                   'a 50 m background gives 12.5 m along the '
                                   'streams.'),
    'grid.quadtree.refine_streams': ('Refine along the streams', _U,
                                     'Needs pyshp: flopy writes the '
                                     'refinement features through a '
                                     'shapefile. Without it the producer says '
                                     'so and builds an unrefined mesh, which '
                                     'is geometrically the base grid.'),
    'grid.override.enable': ('Reproduce an existing grid', _U,
                             'Deriving the grid from the polygon means it is '
                             'not the grid the committed rasters and the '
                             'saved spin-up state were built on. Switch this '
                             'on to rebuild the old one exactly.'),
    'grid.override.xllcorner': ('Origin X', 'm', ''),
    'grid.override.yllcorner': ('Origin Y', 'm', ''),
    'grid.override.nrow': ('Rows', 'count', ''),
    'grid.override.ncol': ('Columns', 'count', ''),

    # ---- panel 2: the surface -----------------------------------------
    'surface.meteo_ts': ('Meteorological record', _U,
                         'ONE file: Date, Time, then SIX columns per station '
                         'in the order P, Ta, RHa, Pa, wind, radiation. '
                         'Hourly. The header row is compulsory and is never '
                         'parsed -- column ORDER is the contract.'),
    'surface.irrigation': ('Irrigation', _U,
                           'Adds the crop and field machinery.'),
    'surface.irr_ts': ('Irrigation series', _U,
                       'Date, Time, then one column per FIELD. Its dates are '
                       'never parsed, so it must be row-aligned with the '
                       'meteorological record.'),
    'surface.crop_schedule': ('Crop schedule pattern', _U,
                              'One file per field. Columns: start, end, '
                              'growing days, wilting days, crop index.'),
    'surface.nfield': ('Irrigated fields', 'count',
                       'Must match the columns of the irrigation series.'),
    'surface.out_prefix': ('Output prefix', _U, 'Names MMsurf\'s own files.'),
    'surface.plot': ('MMsurf figures', _U,
                     '0 writes PNGs to the workspace, 1 opens them on screen.'),
    'surface.meteo_zones': ('Meteo zones', _U,
                            'Thiessen polygons of the stations by default, so '
                            'one station means one zone over the whole '
                            'catchment. Name a layer to override.'),
    'surface.irr_zones': ('Irrigation zones', _U,
                          'A polygon layer with a field_id column: 1..NFIELD, '
                          '0 for not irrigated.'),

    # per-row fields of the surface tables
    'station.name': ('Name', _U, ''),
    'station.x': ('X', 'm', 'Project CRS. The code converts to latitude and '
                            'longitude for the solar geometry.'),
    'station.y': ('Y', 'm', 'Project CRS.'),
    'station.z': ('Altitude', 'm a.s.l.', 'Z in the ini.'),
    'station.tz_lon': ('Time-zone centre longitude', 'deg W', 'Lz.'),
    'station.fuse_shift': ('Time-zone shift', 'h',
                           'FC: the site is not in its nominal time zone.'),
    'station.data_shift': ('Data time shift', 'h',
                           'DTS: the data were not logged at standard clock '
                           'time. NOTE MMsurf applies station 1\'s value to '
                           'every station.'),
    'station.z_wind': ('Wind measurement height', 'm', 'z_m.'),
    'station.z_hum': ('Humidity measurement height', 'm', 'z_h.'),

    'vegetation.name': ('Name', _U, 'One word, no spaces.'),
    'vegetation.h_dry': ('Height, dry season', 'm', 'h_d.'),
    'vegetation.h_wet': ('Height, wet season', 'm', 'h_w.'),
    'vegetation.canopy_mm': ('Canopy storage', 'mm', 'S_w.'),
    'vegetation.c_leaf': ('Max leaf conductance', 'm/s', 'C_leaf_star.'),
    'vegetation.lai_dry': ('LAI, dry season', 'm2/m2',
                           'LAI_d. Zero means the type stops counting as '
                           'vegetation in the dry season and its area reverts '
                           'to bare soil -- which is how grass is handled.'),
    'vegetation.lai_wet': ('LAI, wet season', 'm2/m2', 'LAI_w.'),
    'vegetation.shelter_dry': ('Shelter factor, dry', _U, 'f_s_vd.'),
    'vegetation.shelter_wet': ('Shelter factor, wet', _U, 'f_s_vw.'),
    'vegetation.albedo_dry': ('Albedo, dry', _U, 'alfa_vd.'),
    'vegetation.albedo_wet': ('Albedo, wet', _U, 'alfa_vw.'),
    'vegetation.j_dry': ('Dry season starts', 'julian day', 'J_vd.'),
    'vegetation.j_wet': ('Wet season starts', 'julian day', 'J_vw.'),
    'vegetation.trans_dry_wet': ('Dry-to-wet transition', 'd', 'TRANS_vdw.'),
    'vegetation.trans_wet_dry': ('Wet-to-dry transition', 'd', 'TRANS_vwd.'),
    'vegetation.root_depth': ('Max rooting depth', 'm',
                              'Zr. Also the natural source for the UZF '
                              'extinction depth.'),
    'vegetation.ktg_min': ('Transpiration sourcing, min', _U, 'kTg_min.'),
    'vegetation.ktg_max': ('Transpiration sourcing, max', _U, 'kTg_max.'),
    'vegetation.kt_f': ('Transpiration sourcing, f', _U, 'kT_f.'),
    'vegetation.kt_s': ('Transpiration sourcing, slope s', _U,
                        'This is s, between 0 and 1. The old ini stored 1/s '
                        '-- values above 20 -- and the driver inverted it; '
                        'typing that number here is refused.'),

    'crop.name': ('Name', _U, ''),
    'crop.h': ('Max height', 'm', 'h_c.'),
    'crop.canopy_mm': ('Canopy storage', 'mm', 'S_w_c.'),
    'crop.c_leaf': ('Max leaf conductance', 'm/s', 'C_leaf_star_c.'),
    'crop.lai': ('Max LAI', 'm2/m2', 'LAI_c.'),
    'crop.shelter': ('Shelter factor', _U, 'f_s_c.'),
    'crop.albedo': ('Albedo', _U, 'alfa_c.'),
    'crop.root_depth': ('Max rooting depth', 'm', 'Zr_c.'),
    'crop.ktg_min': ('Transpiration sourcing, min', _U, 'kTg_min_c.'),
    'crop.ktg_max': ('Transpiration sourcing, max', _U, 'kTg_max_c.'),
    'crop.kt_f': ('Transpiration sourcing, f', _U, 'kT_f_c.'),
    'crop.kt_s': ('Transpiration sourcing, slope s', _U, 'As for vegetation.'),

    'soil.name': ('Name', _U, ''),
    'soil.porosity': ('Porosity, top 1 cm', 'm3/m3',
                      'por. The SURFACE soil, for bare-soil evaporation -- '
                      'not the soil column of panel 3.'),
    'soil.field_capacity': ('Field capacity, top 1 cm', 'm3/m3', 'fc.'),
    'soil.albedo_dry': ('Albedo, dry', _U, 'alfa_sd.'),
    'soil.albedo_wet': ('Albedo, wet', _U, 'alfa_sw.'),
    'soil.j_dry': ('Dry season starts', 'julian day', 'J_sd.'),
    'soil.j_wet': ('Wet season starts', 'julian day', 'J_sw.'),
    'soil.trans_dry_wet': ('Dry-to-wet transition', 'd',
                           'TRANS_sdw. MMsurf uses this for both directions.'),
    'soil.trans_wet_dry': ('Wet-to-dry transition', 'd',
                           'TRANS_swd. Read but not used by MMsurf.'),

    'veg_class.code': ('Layer value', _U,
                       'What the vegetation layer\'s class column holds.'),
    'veg_class.veg': ('Vegetation index', _U,
                      '1-based, into the vegetation table above.'),

    # ---- panel 3: soil and MODFLOW ------------------------------------
    'soil.params': ('Soil column parameters', _U,
                    'Per zone and per layer: Smax, Sfc, Sr, Si, Ks. The zone '
                    'ORDER in this file is what the zone codes refer to.'),
    'soil.zones': ('Soil zones', _U,
                   'A polygon layer whose code column matches the zone order '
                   'of the parameter file.'),
    'soil.thickness': ('Soil thickness', _U,
                       'A raster beats a polygon attribute, and a polygon '
                       'attribute beats a single value. Nothing set is an '
                       'error, not a silent zero.'),
    'soil.veg_layer': ('Vegetation layer', _U,
                       'Crown polygons. The share of each cell covered is an '
                       'exact area overlay, so a cell 37 % covered gets 37.'),
    'soil.veg_column': ('Vegetation class column', _U, ''),
    'obs.table': ('Observation points', _U,
                  'Name, x, y, layer and the initial head. "##" before a name '
                  'means the point is not drawn on the maps.'),
    'obs.layer': ('Observation layer', _U, 'Preview only; the table wins.'),
    'layers.nlay': ('Aquifer layers', 'count',
                    '2 and 6 are both authoritative parameter sets, '
                    'maintained by hand -- the 2-layer one is NOT an '
                    'aggregation of the 6-layer one.'),
    'uzf.vks_scale': ('UZF vks multiplier', _U,
                      'Offsets the EPSILON clamp MF6 forces (2.0 -> 3.5).'),
    'seep.kind': ('Seepage mechanism', _U,
                  'drn is the validated choice: a smoothed land-surface '
                  'drain. uzf uses SIMULATE_GWSEEP, which switches on and off '
                  'discontinuously.'),
    'seep.cond': ('Seepage conductance', 'm2/d',
                  'A NUMERICAL device, not a physical property: it must be '
                  'effectively free-draining. It is PER CELL, so it is '
                  'grid-dependent -- retune it when the grid changes.'),
    'et.uzf_et': ('UZF evapotranspiration', _U, 'WP2.'),
    'et.gwet_in_mf': ('Groundwater ET in MODFLOW', _U,
                      'Guarded off: ETg is computed by MARMITES and applied '
                      'as a well sink, so MODFLOW must not remove it as well.'),
    'sfr.enable': ('Stream routing (SFR)', _U,
                   'The network is the mapped hydrography, burned onto the '
                   'grid at run time.'),
    'sfr.min_slope': ('Minimum reach slope', 'm/m', ''),
    'sfr.monotonic_bed': ('Downstream-monotonic bed', _U,
                          'The SFRmaker rule: a bed that rises downstream is '
                          'DEM noise, not topography.'),
    'sfr.width': ('Channel width', 'm',
                  'value, a per-segment column, or the drainage law. The '
                  'drainage law is resolved AFTER routing, because '
                  'contributing area is only known once the reaches are '
                  'ordered.'),
    'sfr.depth': ('Channel incision', 'm',
                  'How far the bed sits below land surface.'),
    'sfr.rhk': ('Streambed conductivity', 'm/d', ''),
    'sfr.rbth': ('Streambed thickness', 'm', ''),
    'sfr.manning': ('Manning\'s n', _U, ''),
    'lak.enable': ('Lakes (LAK)', _U,
                   'One EMBEDDEDV lake per pond: every La Mata pond is '
                   'smaller than a cell, so there is nothing to excavate.'),
    'lak.depth': ('Pond depth', 'm', 'Below the rim.'),
    'lak.bedleak': ('Lakebed leakance', '1/d', 'Clay-lined charca: 1e-3.'),
    'lak.surfdep': ('Surface depression depth', 'm',
                    'Smooths the wetted area as the pond dries.'),
    'crr.enable': ('Runoff cascade (CRR)', _U,
                   'Routes runoff downslope cell to cell instead of losing '
                   'it at the cell it was generated in.'),
    'crr.beta': ('CRR beta', _U, 'Daoud et al. (2022), Eq. 23.'),
    'spinup.cycles': ('Spin-up cycles', 'count', ''),
    'spinup.strt_heads': ('Initial heads', _U,
                          'A saved, equilibrated head field. Blank starts '
                          'from the configured initial condition, which for a '
                          'cold start puts the water table above ground over '
                          'much of the catchment.'),
    'spinup.steady_means': ('Steady-state means', _U,
                            'Per-cell recharge and groundwater ET driving the '
                            'steady first period.'),

    # ---- panel 4: plots ------------------------------------------------
    'postproc.enable': ('Post-process', _U, ''),
    'postproc.preproc': ('Input maps', _U, 'Draw the inputs as maps first.'),
    'postproc.only': ('Post-process only', _U,
                      'Re-draw from a finished run without re-running it.'),
    'postproc.hydro_year_start': ('Hydrological year starts', 'month',
                                  'Drives the x-axis of every time series and '
                                  'the Sankey year index. It was in the ini '
                                  'but never reached the model, which used a '
                                  'hardcoded October.'),
    'postproc.wb_unit': ('Water-balance unit', _U, 'year or day.'),
    'postproc.obs_series': ('Observation-point series', _U, ''),
    'postproc.sankey': ('Sankey diagrams', _U,
                        'The water balance as a flow diagram, per '
                        'hydrological year and for the whole period.'),
    'postproc.sankey_min_flux': ('Sankey flux threshold', 'mm/y',
                                 'Flows below this are not drawn.'),
    'postproc.input_maps': ('Input maps', _U, ''),
    'postproc.result_maps': ('Result maps', _U, ''),
    'postproc.map_days': ('Map days', 'count',
                          '0 draws the time mean only.'),
    'postproc.tick_trimester_years': ('Trimester ticks up to', 'years', ''),
    'postproc.tick_semester_years': ('Semester ticks up to', 'years', ''),
    'postproc.sankey_full': ('Whole-period Sankey', _U,
                            'A panel for the whole run as well as one per '
                            'hydrological year.'),
    'postproc.sankey_obs_years': ('Per-point yearly Sankey', _U,
                                  'One panel per year at each observation '
                                  'point. Many figures.'),

    # ---- the rest of panel 3 ------------------------------------------
    'grid.voronoi.trans_levels': ('Transition cell sizes', 'm',
                                  'Graded steps between the stream corridor '
                                  'and the background, Daoud-style.'),
    'layers.aggregate': ('Derive 2 layers from 6', _U,
                         'Comparison only. The 2-layer parameter set is '
                         'maintained by hand and is NOT this.'),
    'obs.name_column': ('Name column', _U, 'Of the optional point layer.'),
    'obs.heads_prefix': ('Head series prefix', _U,
                         '<prefix>_<point>.txt, one file per piezometer.'),
    'obs.sm_prefix': ('Soil-moisture series prefix', _U, ''),
    'obs.ro_prefix': ('Runoff series prefix', _U, ''),
    'et.unsat_form': ('Unsaturated ET form', _U,
                      'etwc uses water content, etae capillary pressure.'),
    'et.extdp_source': ('Extinction depth from', _U,
                        'uniform, the vegetation rooting depth, or a raster.'),
    'et.extdp_default': ('Default extinction depth', 'm', ''),
    'et.extwc_source': ('Extinction water content from', _U, ''),
    'sfr.source': ('Reach table', _U,
                   'The grid-independent stream geometry the converter '
                   'writes from the hydrography layer.'),
    'lak.source': ('Pond table', _U,
                   'Centroids and DEM statistics, one row per pond.'),
    'lak.geometry': ('Pond polygons', _U,
                     'What the builder fits an embedded lake to. The pond '
                     'TABLE is not enough: it needs the footprint.'),
    'lak.polygons': ('Pond layer', _U,
                     'In DATA_ROOT/GIS, read by the converter only.'),
    'lak.maxiter': ('LAK Newton iterations', 'count', ''),
    'lak.stagechg': ('LAK stage tolerance', 'm', ''),
    'crr.sinks': ('Topographic sinks', _U,
                  'What happens to runoff that reaches a closed depression.'),
    'crr.dem': ('Sink-filled DEM', _U, 'Drives the downslope cascade.'),
    'spinup.tol': ('Spin-up tolerance', 'm',
                   'Water-table change between cycles below which the spin-up '
                   'is considered converged.'),
    'spinup.save_strt': ('Save the final heads as', _U,
                         'Blank does not save. Written to the workspace.'),
    'spinup.save_means': ('Save the steady means as', _U, ''),
    'spinup.strt_dem': ('Initial head from the DEM', _U,
                        '[a, b] gives head = a*elevation + b. Empty uses the '
                        'configured initial condition.'),
}


# A field whose value is one of a fixed set is a CHOICE, not free text.
# Typing "voroni" into a text box and finding out at run time is exactly the
# failure the panels exist to prevent, so the allowed values come from the
# schema itself wherever it declares them.
def _grid_kinds():
    import marmites_config as mcfg
    return list(mcfg.GRID_KINDS)


def _resample_modes():
    import marmites_config as mcfg
    return list(mcfg.RESAMPLE_MODES)


CHOICES = {
    'grid.kind': _grid_kinds,
    'grid.resample': _resample_modes,
    'run.mode': lambda: ['lagged', 'iterative'],
    'seep.kind': lambda: ['uzf', 'drn'],
    'et.unsat_form': lambda: ['etwc', 'etae'],
    'et.extdp_source': lambda: ['uniform', 'veg_zone', 'raster'],
    'postproc.wb_unit': lambda: ['year', 'day'],
    'crr.sinks': lambda: ['evaporate', 'route'],
    'ui.execution': lambda: ['local', 'server'],
}

# Blocks that only apply for one value of another field. The panel shows them
# as a SUB-PANEL under that choice, rather than as settings that look live and
# are not.
#   dotted prefix -> (controlling field, values that make it apply)
SUBPANELS = {
    'grid.voronoi': ('grid.kind', ('voronoi',)),
    'grid.quadtree': ('grid.kind', ('quadtree',)),
}


def choices_for(dotted):
    """The allowed values of an enumerated field, or None."""
    f = CHOICES.get(dotted)
    if f is None:
        return None
    try:
        return list(f())
    except Exception:                                    # pragma: no cover
        return None


def subpanel_for(dotted):
    """(controlling field, applicable values) if this is a conditional block."""
    for prefix, rule in SUBPANELS.items():
        if dotted == prefix or dotted.startswith(prefix + '.'):
            return rule
    return None


def panel_of(section):
    """Which panel a configuration section belongs to, or None."""
    for num, _t, _i, _s, sections, _b in PANELS:
        if section in sections:
            return num
    return None


def describe(dotted):
    """(label, units, help) for a dotted key, falling back to the key itself.

    A row field of an array of tables is looked up by its LAST two parts --
    ``surface.vegetation[2].kt_s`` and ``vegetation.kt_s`` are the same field.
    """
    import re

    if dotted in FIELDS:
        return FIELDS[dotted]
    # Strip the row index wherever it sits: surface.vegetation[2].kt_s is the
    # same field as vegetation.kt_s.
    plain = re.sub(r'\[\d+\]', '', dotted)
    if plain in FIELDS:
        return FIELDS[plain]
    parts = plain.split('.')
    if len(parts) >= 2:
        tail = '%s.%s' % (parts[-2], parts[-1])
        if tail in FIELDS:
            return FIELDS[tail]
    return (parts[-1], '', '')


def is_source(v):
    """Is this a ParamSource / VectorSource -- one value with a producer?

    Detected by behaviour rather than by name, so a third source type added
    later is handled without editing this. These are shown as ONE field: the
    producer is a choice (value | column | raster | layer | drainage), and
    splitting it into five widgets would invite setting two at once, which
    the schema then refuses.
    """
    import dataclasses

    return dataclasses.is_dataclass(v) and hasattr(v, 'producer')


def fields_of(cfg, section):
    """[(dotted, value)] for a section, nested blocks flattened.

    Arrays of tables are NOT returned -- they are edited as tables, which is
    what ``TABLES`` is for -- and a ParamSource / VectorSource stays whole.
    """
    import dataclasses

    sec = getattr(cfg, section, None)
    if sec is None:
        raise PanelError('the configuration has no section %r' % section)
    out = []
    for f in dataclasses.fields(sec):
        v = getattr(sec, f.name)
        key = '%s.%s' % (section, f.name)
        if isinstance(v, list) and v and dataclasses.is_dataclass(v[0]):
            continue
        if is_source(v):
            out.append((key, v))
            continue
        if dataclasses.is_dataclass(v):
            for g in dataclasses.fields(v):
                out.append(('%s.%s' % (key, g.name), getattr(v, g.name)))
            continue
        out.append((key, v))
    return out
