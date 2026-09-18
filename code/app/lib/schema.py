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
           'panel_name', 'HIDDEN',
           'describe', 'fields_of', 'choices_for', 'subpanel_for', 'is_source',
           'GRID_PERMANENT', 'GRID_SUBPANEL', 'GRID_DERIVED', 'GRID_GATED',
           'GRID_NEEDS_LAYER', 'GRID_NOTES', 'grid_fields',
           'SURFACE_ROWS',
           'SURFACE_FILES', 'SURFACE_PATTERNS', 'SURFACE_ON_PLOTS',
           'SURFACE_GATED', 'TIME_ROWS',
           'SURFACE_TABLE_FILES', 'INTEGER_VALUE',
           'SOIL_ROWS', 'SOIL_VEG_ROWS', 'SOIL_FILES',
           'SOIL_DATASET_FILES', 'COLUMN_OF', 'GEOMETRY_ROWS',
           'GHB_ROWS', 'DRN_ROWS', 'LAYER_LIST',
           'UZF_ROWS', 'UZF_ET_ROWS', 'UZF_LEGACY',
           'OBS_COMMON',
           'OBS_GROUPS',
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
     'Introduce here the catchment boundary: it must be a polygon shape file '
     'with projected coordinates (metric). Next select the type of grid and '
     'its parameters, test and visualize. When you are satisfied, select it '
     'as the base of the MM-MF model. All input will be wrapped onto this '
     'selected grid.'),
    (2, 'Surface and driving forces', '🌦️', 'run.surface', ['surface'],
     'MMsurf turns the hourly meteorological record into the daily forcing '
     'MMsoil consumes. With the switch off, that forcing must already exist '
     'and is checked for presence and shape before the run starts.'),
    # Panel 3 used to be one panel called Model, carrying the soil column,
    # the aquifer, the surface-water packages and the observations together.
    # They are MMsoil's inputs and MODFLOW's inputs, which is one model but
    # two different sets of questions, and the modeller answers them at
    # different times. Split three ways.
    (3, 'Soil', '🌱', 'run.model', ['soil'],
     'What MMsoil reads: the soil column zones and thickness, and the '
     'vegetation cover each cell carries.'),
    (4, 'Unsaturated zone and groundwater', '🌍', 'run.model',
     ['layers', 'ghb', 'drn', 'uzf', 'seep', 'et', 'sfr', 'lak', 'crr',
      'spinup'],
     'What MODFLOW 6 reads: the layers, the unsaturated zone, the seepage '
     'face, evapotranspiration, the surface-water packages and the initial '
     'state. The switch is the SAME one as the Soil panel\'s -- the soil water '
     'balance and MODFLOW are no longer separable, since MMsoil is '
     'stepped from inside the MODFLOW time loop.'),
    (5, 'State variables (calibration)', '📏', None, ['obs'],
     'What was MEASURED, against what the model produces: groundwater '
     'levels, soil moisture, actual evapotranspiration and surface runoff. '
     'Nothing here changes a flux -- it decides what the run is compared '
     'with.'),
    (6, 'Plots', '📊', 'run.plot', ['postproc'],
     'What to draw once the run finishes. Nothing here changes a flux.'),
]


def panel_name(number):
    """A panel's NAME, for anything the modeller reads.

    Never its number. The sidebar shows names, the numbers move whenever a
    panel is split -- three of them moved when Model became Soil,
    Unsaturated zone and State variables -- and a sentence carrying a number
    is a sentence that will one day point somewhere else. Built from PANELS,
    so a rename cannot leave prose behind.
    """
    for p in PANELS:
        if p[0] == number:
            return p[1]
    return 'panel %d' % number


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
    'run.nsp': ('Number of stress periods for test', 'count',
                '0 runs the whole record. Use a small number to try a change '
                'before committing to the full run.'),
    'run.daily': ('One stress period per day', _U,
                  'On, the run steps day by day, as the record does. Off '
                  'aggregates it by RAINFALL EPISODE: a day with rain gets a '
                  'period of its own and the dry days after it are averaged '
                  'together, up to the maximum beside this box.'),
    'run.perlen_max': ('Longest aggregated period', 'd',
                       'How many dry days may be averaged into one stress '
                       'period. Ignored while the run is daily. This was the '
                       'MODFLOW ini\'s `nper`, which is NOT a count of '
                       'periods -- the count is worked out from the rainfall '
                       'and cannot be known before the record is read.'),
    'run.ats': ('Adaptive time stepping', _U,
                'On. A period MODFLOW cannot solve in one step is not a '
                'result, and without this it silently becomes one.'),
    'run.build_only': ('Write the input files and stop', _U,
                       'Useful to inspect what MODFLOW would be given.'),
    'run.max_discrepancy': ('Maximum mass-balance discrepancy', '%',
                            'The run fails above this cumulative value.'),

    # ---- panel 1: the grid --------------------------------------------
    'grid.boundary': ('Catchment boundary (polygon)', _U,
                      'A projected, metric shapefile that defines the active '
                      'domain, the grid boundary and the model rectangle.'),
    'grid.streams': ('Hydrography layer (lines)', _U,
                     'The shape file with the stream network. It refines the '
                     'grid corridor (voronoi and quadtree), and becomes the '
                     'MF6 SFR network. Leave it blank if the catchment has '
                     'none (no refinement).'),
    'grid.ponds': ('Pond layer (polygons)', _U,
                   'The shape file with the ponds limits (charcas). Seeds a '
                   'grid cell per pond and becomes the LAK footprints. If '
                   'blank, the LAK package is not activated.'),
    'grid.dem': ('Elevation raster (DEM)', _U,
                 'The sink-filled DEM, in any format GDAL reads (asc, '
                 'GeoTIFF, or an ESRI grid folder). It is the land surface. '
                 'The DEM dataset is kept at its own resolution and the model '
                 'wraps it onto the grid cells, area-weighted. It also gives '
                 'each pond its rim and bottom.'),
    'grid.crs_epsg': ('Project CRS', 'EPSG',
                      'The projected, metric CRS of the project, extracted '
                      'from file .prj of the catchment layer.'),
    'grid.cell_size': ('Cell size', 'm',
                       'The cell itself on a structured or disv grid, and the '
                       'background GRIDGEN halves on a quadtree. A Voronoi '
                       'mesh does not use it at all -- there the size is '
                       'cell_far, and the rectangle is snapped to that.'),
    'grid.buffer': ('Buffer', 'm',
                    'Extend the model rectangle beyond the polygon.'),
    'grid.kind': ('Grid kind', _U,
                  'Grid type of the model. voronoi is the default.'),
    'grid.resample': ('Resampling rule', _U,
                      'auto area-weights: continuous fields and '
                      'majority-votes zone maps. centre: samples the cell '
                      'centre (worse per cell, but it preserves zone '
                      'proportions).'),
    'grid.voronoi.cell_pond': ('Cell size at a pond', 'm',
                               'The size a POND\'s own cell carries. 0 lets '
                               'the pond take whatever the corridor gives it '
                               '-- its footprint is refined with the stream '
                               'network, so it is meshed at the cell size '
                               'near the stream and grades out with it. Set '
                               'it and each pond gets a zone of its own at '
                               'this size, which is the ONLY way to make a '
                               'pond cell COARSER than the corridor around '
                               'it: put roughly the pond\'s own width here '
                               'for one cell per pond. It cannot exceed the '
                               'background size.'),
    'grid.quadtree.refine_ponds': ('Refine at the ponds', _U,
                                   'GRIDGEN splits every cell a feature '
                                   'touches, so the pond footprints refine '
                                   'the cells they cover whether or not a '
                                   'mapped stream runs through them. Off '
                                   'until a pond layer is given.'),
    'grid.quadtree.pond_level': ('Refinement level at the ponds', 'levels',
                                 'Each level halves the cell, so one more '
                                 'than the streams is a quarter of the cell '
                                 'area. 0 means the same level as the '
                                 'streams.'),
    'grid.voronoi.cell_far': ('Background cell size', 'm',
                              'The side of the EQUIVALENT SQUARE, so 100 aims '
                              'at 10 000 m2. Every build prints the area it '
                              'actually achieved -- read that, not this.'),
    'grid.voronoi.cell_near_stream': ('Cell size near the stream', 'm',
                                      'The size the innermost band carries. '
                                      'Cleared while the refinement is off.'),
    'grid.voronoi.stream_buffer': ('Stream corridor width (derived)', 'm',
                                   'Distance from the centreline at which the '
                                   'cells have reached the background size. '
                                   'DERIVED: each band is as wide as its own '
                                   'cells, so the corridor is their sum. '
                                   'Read-only.'),
    'grid.voronoi.stream_refine': ('Refine along the streams', _U,
                                   'Off is the SFRmaker approach: map the '
                                   'network onto the background grid. CdL '
                                   'settled on off after judging a refined '
                                   'mesh too refined near streams.'),
    'grid.voronoi.grade_ratio': ('Maximum size ratio between bands', _U,
                                 'How fast the cells may grow from the '
                                 'corridor outwards, and so how many '
                                 'transition bands there are and how wide '
                                 'the corridor becomes. 3, the default, '
                                 'grades in the fewest size changes; 1.5 to '
                                 '2 is the gentler rule of thumb and takes '
                                 'more bands. Held in (1, 3].'),
    'grid.voronoi.refine_ponds': ('Refine at the ponds', _U,
                                 'Two behaviours: (1) The footprints join '
                                 'the geometry the corridor bands are '
                                 'buffered from, so a pond is meshed at the '
                                 'cell size near the stream and grades out '
                                 'with it -- which is what reaches a charca '
                                 'that no mapped stream runs through. (2) A '
                                 'GENERATOR sits at each pond centre, so one '
                                 'cell is centred on the water instead of '
                                 'straddling it: that is the cell LAK will '
                                 'be hung on. ONE CELL PER POND is not this '
                                 'switch -- it is the pond cell size below, '
                                 'set to about the pond\'s own width.'),
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
                      'not the soil column of the %s panel.'
                      % panel_name(3)),
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

    # ---- panels 3 to 5: soil, subsurface, calibration -----------------
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
                    'flopy: `ModflowGwfdis(nlay=)`, or `ModflowGwfdisv` on a '
                    'mesh. 2 and 6 are both authoritative parameter sets, '
                    'maintained by hand -- the 2-layer one is NOT an '
                    'aggregation of the 6-layer one.'),
    'layers.hnoflo': ('No-flow / dry value', _U,
                      'The number that marks a cell as having nothing to '
                      'report. MARMITES masks on it too, so the two must be '
                      'the SAME number -- which is why it is asked once here '
                      'rather than twice in two files. It must not be 0: that '
                      'is a perfectly good head.'),
    'layers.thickness': ('Layer thickness', 'm',
                         'flopy: becomes `botm`, subtracted layer by layer '
                         'from the aquifer top. ONE RASTER PER LAYER: write '
                         '`%d` where the layer number goes and it is replaced '
                         'by 1, 2, ... up to the number of layers above -- '
                         '`thick_%d.asc` with 2 layers opens `thick_1.asc` '
                         'and `thick_2.asc`. `%d` is the whole placeholder, '
                         'so no `d` is left behind, and the count has no '
                         'leading zero. A name without `%d` is that one '
                         'raster for every layer, and a single value is '
                         'uniform over all of them. The panel prints the '
                         'expansion and whether the files are there.'),
    'layers.k': ('Hydraulic conductivity', 'm/d',
                 'flopy: `ModflowGwfnpf(k=)` -- the HORIZONTAL conductivity. '
                 'Per layer, by the same `%d` rule as the thickness.'),
    'layers.k33': ('Vertical conductivity', 'm/d or ratio',
                   'flopy: `ModflowGwfnpf(k33=)`. With **as a ratio** on '
                   '(the default, and the legacy `LAYVKA = 1`) the number is '
                   'the ANISOTROPY k/k33, so 2 means the aquifer conducts '
                   'water half as easily downwards as sideways and 1 is '
                   'isotropic; with it off the number is a conductivity in '
                   'm/d. Layering makes the ratio > 1 -- below 1 says water '
                   'moves more easily vertically than horizontally, which is '
                   'rarely what is meant. Per layer, by the same `%d` rule.'),
    'layers.k33_as_ratio': ('k33 given as a ratio', _U,
                            'On: the number above is k/k33 (the legacy '
                            '`LAYVKA = 1`). Off: it is a vertical '
                            'conductivity in m/d. It applies to every layer '
                            '-- every parameter set in this repository sets '
                            'it the same way for all of them.'),
    'layers.convertible': ('Convertible layers', _U,
                           'flopy: `ModflowGwfnpf(icelltype=)` and '
                           '`ModflowGwfsto(iconvert=)`, 1 when on. The water '
                           'table sits inside the modelled stack, so '
                           'transmissivity has to follow the saturated '
                           'thickness and storage has to switch between Sy '
                           'and Ss as a cell drains. Off means every layer is '
                           'confined. `layavg`, `laywet`, `laycbd` and `hdry` '
                           'are not asked: MODFLOW 6 has none of them.'),
    'layers.ss': ('Specific storage', '1/m',
                  'flopy: `ModflowGwfsto(ss=)`. Per layer, by the same `%d` '
                  'rule as the thickness.'),
    'layers.sy': ('Specific yield', _U,
                  'flopy: `ModflowGwfsto(sy=)`. Per layer, by the same `%d` '
                  'rule as the thickness.'),
    'ghb.enable': ('General-head boundary (GHB)', _U,
                   'Builds `ModflowGwfghb`. Off is `ghb_yn = 0` in the '
                   'parameter file, and nothing below is read.'),
    'ghb.layers': ('MODFLOW layers', 'list',
                   'Which layers carry the boundary, counted from 1 the way '
                   'MODFLOW counts them. A per-layer raster (`%d`) gives each '
                   'listed layer its own map; one value or one shapefile '
                   'applies to all of them.'),
    'ghb.head': ('Boundary head', 'm a.s.l.',
                 'flopy: `ModflowGwfghb` stress period data, `bhead`. WHERE '
                 'the boundary is comes from this source: a cell carries a '
                 'GHB where it produces a value and none where it falls back '
                 'to `fill`. That is the legacy rule -- the rasters are zero '
                 'except along the boundary -- said once instead of per '
                 'raster.'),
    'ghb.cond': ('Boundary conductance', 'm²/d',
                 'flopy: `ModflowGwfghb` stress period data, `cond`. The '
                 'conductance of the material between the cell and the head '
                 'outside it. Must be > 0.'),
    'drn.enable': ('Drains (DRN)', _U,
                   'Builds `ModflowGwfdrn`. This is the BOUNDARY drain the '
                   'parameter file called `drn` -- in La Mata six cells at '
                   'the catchment outlet. It is NOT the seepage face: that '
                   'is a second drain package, `drn_seep`, over the whole '
                   'land surface, configured lower down this same tab.'),
    'drn.layers': ('MODFLOW layers', 'list',
                   'Which layers carry the drain, counted from 1.'),
    'drn.elevation': ('Drain elevation', 'm a.s.l.',
                      'flopy: `ModflowGwfdrn` stress period data, `elev`. '
                      'Water leaves a cell once its head rises above this. '
                      'WHERE the drain is comes from this source, as for the '
                      'GHB head. Ignored when the elevation is taken at the '
                      'base of the layer below.'),
    'drn.at_layer_base': ('At the base of the layer', _U,
                          'The drain sits just above the bottom of its own '
                          'layer (`botm` + 0.01 m) instead of at the '
                          'elevation above. This is the convention the '
                          'parameter file used, where a NEGATIVE elevation '
                          'meant the layer base, '
                          'which is how every La Mata drain cell is written '
                          '(-1) -- as a switch it says what it does, as -1 '
                          'in a raster it had to be learnt from a comment.'),
    'drn.cond': ('Drain conductance', 'm²/d',
                 'flopy: `ModflowGwfdrn` stress period data, `cond`. Must be '
                 '> 0.'),
    'uzf.ntrailwaves': ('Trail waves', 'count',
                        'flopy: `ModflowGwfuzf(ntrailwaves=)`, the UZF1 '
                        '`ntrail2`. How finely a drying front is resolved. '
                        '7 is the MODFLOW 6 default and what the build has '
                        'been using; the parameter file asked for 15.'),
    'uzf.nwavesets': ('Wave sets', 'count',
                      'flopy: `ModflowGwfuzf(nwavesets=)`, the UZF1 `nsets`. '
                      'How many wetting/drying episodes a column can hold at '
                      'once. 40 is the MODFLOW 6 default; the parameter file '
                      'asked for 500.'),
    'uzf.surfdep': ('Surface depression depth', 'm',
                    'flopy: `ModflowGwfuzf` packagedata `surfdep`. The relief '
                    'within a cell, which smooths the onset of surface '
                    'leakage. Must be smaller than the cell thickness.'),
    'uzf.eps': ('Brooks-Corey exponent', _U,
                'flopy: packagedata `eps`. K(theta) = vks * Se^eps, so a '
                'bigger number throttles recharge through a deep unsaturated '
                'zone. **MODFLOW 6 enforces 3.5 to 14.0**; UZF1 accepted 2.0, '
                'which is what the La Mata NWT model used -- `vks multiplier` '
                'below exists to offset that clamp.'),
    'uzf.thtr': ('Residual water content', 'm³/m³',
                 'flopy: packagedata `thtr`. **UZF6 requires it to be > 0** '
                 'whatever the UZF1 `SPECIFYTHTR` option said, so there is no '
                 'switch for it any more -- it is always used.'),
    'uzf.thts': ('Saturated water content', 'm³/m³',
                 'flopy: packagedata `thts`. Above `thtr`, at most 1.'),
    'uzf.thti': ('Initial water content', 'm³/m³',
                 'flopy: packagedata `thti`, between `thtr` and `thts`. '
                 'Always carried, so the UZF1 `SPECIFYTHTI` option is gone '
                 'too.'),
    'uzf.vks_from': ('Unsaturated vertical K from', _U,
                     'The UZF1 `iuzfopt` with the numbers replaced by what '
                     'they meant. **layer** (iuzfopt = 2) uses each layer\'s '
                     'own `k33`, so the unsaturated column and the aquifer '
                     'cannot disagree. **raster** (iuzfopt = 1) reads the map '
                     'below instead.'),
    'uzf.vks': ('Unsaturated vertical K', 'm/d',
                'flopy: packagedata `vks`. Read only when the source above is '
                '**raster**. Per layer, by the same `%d` rule as the '
                'aquifer properties.'),
    'uzf.vks_scale': ('UZF vks multiplier', _U,
                      'NOT a MODFLOW field: it multiplies the `vks` column '
                      'before it is written. MODFLOW 6 forbids eps < 3.5 and '
                      'the NWT model ran at 2.0; a higher exponent throttles '
                      'recharge, and raising vks restores the rate the NWT '
                      'model had. It scales the UZF column ONLY, never the '
                      'aquifer `k33`.'),
    'seep.kind': ('Seepage mechanism', _U,
                  'drn is the validated choice: a smoothed land-surface '
                  'drain. uzf uses SIMULATE_GWSEEP, which switches on and off '
                  'discontinuously.'),
    'seep.cond': ('Seepage conductance', 'm2/d',
                  'A NUMERICAL device, not a physical property: it must be '
                  'effectively free-draining. It is PER CELL, so it is '
                  'grid-dependent -- retune it when the grid changes.'),
    'et.uzf_et': ('UZF evapotranspiration', _U,
                  'flopy: `ModflowGwfuzf(simulate_et=)`. Lets UZF dry the '
                  'unsaturated zone itself. GROUNDWATER ET is never done '
                  'here: MARMITES computes ETg and applies it through the '
                  'WEL package, so MODFLOW must not remove it a second '
                  'time.'),
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

    # ---- panel 6: plots ------------------------------------------------
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

    # ---- the rest of panels 3 to 5 ------------------------------------
    # DERIVED, and shown read-only: these are DISTANCES from the stream
    # centreline, not cell sizes -- the old label said sizes, which is how a
    # 70 m band ended up sitting outside a 60 m corridor.
    'grid.voronoi.trans_levels': ('Transition bands (derived)', 'm',
                                  'Buffer distances from the centreline at '
                                  'which the cell size steps up, computed '
                                  'from the corridor and the maximum size '
                                  'ratio. Read-only: editing it does '
                                  'nothing.'),
    'obs.name_column': ('Name column', _U, 'Of the optional point layer.'),
    'obs.heads_prefix': ('Head series prefix', _U,
                         '<prefix>_<point>.txt, one file per piezometer.'),
    'obs.sm_prefix': ('Soil-moisture series prefix', _U,
                      '<prefix>_<point>.txt, as for the heads.'),
    'obs.aet_prefix': ('Actual ET series prefix', _U,
                       'Measured actual evapotranspiration -- an eddy tower, '
                       'a lysimeter, or a remote-sensing product. Blank where '
                       'there is none, which is the usual case.'),
    'obs.ro_prefix': ('Runoff series prefix', _U,
                      'At the catchment outlet, so one series rather than one '
                      'per point.'),
    'et.unsat_form': ('Unsaturated ET form', _U,
                      'etwc uses water content, etae capillary pressure.'),
    'et.extdp': ('Extinction depth', 'm',
                 'flopy: `ModflowGwfuzf` perioddata `extdp`. How deep the '
                 'unsaturated zone can be dried by evapotranspiration. By the '
                 'usual rule: a raster, a COLUMN OF THE VEGETATION LAYER (as '
                 'CdL does it, so each species carries its own rooting '
                 'depth), or one value for the whole catchment. Read only '
                 'when UZF evapotranspiration is on.'),
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
    'crr.dem': ('Sink-filled DEM', _U,
                'Drives the downslope cascade. The DATASET copy, at the '
                'survey resolution, wrapped onto the model cells like every '
                'other input -- written by the converter from [grid] dem.'),
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


# A source whose single value is a ZONE, not a measurement: "the whole
# catchment is zone 1". Typed as a count, because 1.0 zones is not a thing
# and a box that offers decimals invites one. soil.thickness is the
# counter-example and is deliberately NOT here -- its value is metres.
INTEGER_VALUE = ('surface.meteo_zones', 'surface.irr_zones', 'soil.zones')


# Never OFFERED on a panel: the configuration carries it, validate() pins it,
# and a box with one legal value is not a question -- it is an invitation to
# get an error message. `et.gwet_in_mf` must stay false because MARMITES
# computes ETg and applies it as a well sink; MODFLOW removing it as well
# would take the water twice.
HIDDEN = ()

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
    # The ini's iuzfopt, with the numbers replaced by what they meant.
    'uzf.vks_from': lambda: ['layer', 'raster'],
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
    # WP1d panel 1: the settings that used to sit in the permanent block and
    # are in fact structured-only. The override reproduces an EXISTING
    # rectangle, which only means something for the two kinds the legacy
    # comparison uses -- validate() refuses it for the others rather than
    # applying a box that is no longer on screen.
    'grid.override': ('grid.kind', ('structured', 'disv')),
    'grid.cell_size': ('grid.kind', ('structured', 'disv', 'quadtree')),
}

# What each kind's sub-panel shows, ROW BY ROW. The nesting is the layout: a
# switch gets a row to itself and what it controls goes on the row below, so
# the dependency is visible without reading the help. Everything NOT listed
# here is in the permanent block; a test checks the two partition the [grid]
# block with nothing left over.
# The catchment comes first, and the two layers a GRID can depend on come
# with it -- the refinement options below cannot be answered before it is
# known whether there IS a network or a pond to refine around.


# ---- panel 4: one sub-panel per MF6 package ------------------------------
# MODFLOW AQUIFER LAYERS first. The top is NOT asked -- it is the land surface minus the
# soil column, so a third answer could only disagree with the other two --
# and neither is ibound: a cell is active when it is inside the catchment
# polygon the grid was built in.
GEOMETRY_ROWS = (
    ('layers.nlay', 'layers.hnoflo'),
    ('layers.thickness', None),
    # k and k33 sit on one row, with the switch that says what the k33
    # number MEANS directly beneath it: 2 is an anisotropy ratio or a
    # conductivity in m/d depending on that one box, and a reader who
    # cannot see them together cannot tell which.
    ('layers.k', 'layers.k33'),
    (None, 'layers.k33_as_ratio'),
    ('layers.ss', 'layers.sy'),
    ('layers.convertible', None),
)

# ---- panel 4: the boundary packages --------------------------------------
# The switch on its own row and what it controls below it, the way the grid
# sub-panels do it: a reader can see at a glance that nothing under the
# switch is read when it is off.
# Fields that are a LIST OF MODFLOW LAYERS. They get a multiselect built
# from layers.nlay rather than the "edit in the TOML tab" a bare list falls
# back to -- a boundary that applies to no layer is not buildable, so the
# panel has to be able to say which.
LAYER_LIST = ('ghb.layers', 'drn.layers')

GHB_ROWS = (
    ('ghb.enable', None),
    ('ghb.layers', None),
    ('ghb.head', 'ghb.cond'),
)
DRN_ROWS = (
    ('drn.enable', None),
    ('drn.layers', None),
    ('drn.elevation', 'drn.cond'),
    ('drn.at_layer_base', None),
)

# ---- panel 4: the unsaturated zone ---------------------------------------
# A vertical column of UZF objects per active cell. The rows group what is
# answered together: the wave machine, then the Brooks-Corey contents, then
# where the unsaturated conductivity comes from -- with the switch directly
# above the map it decides whether to read.
UZF_ROWS = (
    ('uzf.ntrailwaves', 'uzf.nwavesets'),
    ('uzf.surfdep', 'uzf.eps'),
    ('uzf.thtr', 'uzf.thts'),
    ('uzf.thti', None),
    ('uzf.vks_from', 'uzf.vks_scale'),
    ('uzf.vks', None),
)
# Evapotranspiration inside MODFLOW, which is UNSATURATED-zone ET only.
UZF_ET_ROWS = (
    ('et.uzf_et', None),
    ('et.extdp', 'et.extwc_source'),
    ('et.unsat_form', None),
)

# WHERE EVERY UZF1 NAME WENT. The parameter file carried twenty-two of
# them and MODFLOW 6 keeps eight. Written down because the question comes
# up every time someone opens the old file, and because "it is gone" is
# only a useful answer with the reason attached.
# (legacy name, what it is now)
UZF_LEGACY = (
    ('ntrail2', 'ntrailwaves, above'),
    ('nsets', 'nwavesets, above'),
    ('surfdep', 'surfdep, above'),
    ('eps', 'eps, above (now clamped to 3.5-14)'),
    ('thts, thti', 'above'),
    ('vks', 'above, when the source is a raster'),
    ('iuzfopt', '"Unsaturated vertical K from": 1 was raster, 2 the layer'),
    ('finf_user', 'perioddata `finf`; the coupled run overwrites it daily '
                  'from MMsoil'),
    ('[SPECIFYTHTR]', 'gone — UZF6 always needs thtr > 0'),
    ('[SPECIFYTHTI]', 'gone — thti is always carried'),
    ('[NOSURFLEAK]', 'gone — no equivalent'),
    ('irunflg', 'gone — rejected infiltration is routed by MVR/SFR/DRN'),
    ('ietflg', 'gone — this is the UZF evapotranspiration switch below'),
    ('iuzfcb1, iuzfcb2', 'gone — `budget_filerecord` and `save_flows`'),
    ('nuzgag, iuzrow, iuzcol, iftunit, iuzopt',
     'gone — UZF1 gage plumbing; MF6 writes observation files'),
    ('nuztop', 'gone — objects attach to explicit cell ids, so which layer '
               'is at the surface is the cell list'),
    ('uzf_iuzfbnd', 'gone — the footprint is the active catchment cells'),
)

# ---- panel 5: the measured state variables -------------------------------
# The four the model produces, each its own block. They used to be three
# boxes in a row with no heading, which said nothing about WHAT is being
# compared -- and the fourth, actual evapotranspiration, was not asked for at
# all even though the model computes it.
OBS_COMMON = ('obs.table', 'obs.name_column', 'obs.layer')
OBS_GROUPS = (
    ('Groundwater levels', 'obs.heads_prefix', 'm a.s.l., one file per '
     'piezometer'),
    ('Soil moisture', 'obs.sm_prefix', 'per soil zone and depth'),
    ('Actual evapotranspiration', 'obs.aet_prefix', 'an eddy tower, a '
     'lysimeter, or a remote-sensing product'),
    ('Surface runoff', 'obs.ro_prefix', 'at the catchment outlet'),
)


# ---- the time discretisation, on the driving-forces panel ---------------
# WHERE THE CALENDAR COMES FROM. The days themselves are MMsurf's: it reads
# the hourly record and writes one row per day into inputDATE.txt. What the
# run then does with those days -- one stress period each, or rainfall
# episodes averaged together -- is a question about the RECORD, so it is
# asked beside it rather than buried in the MODFLOW parameter file, where
# the aggregation limit sat under the misleading name `nper`.
TIME_ROWS = (
    ('run.daily', 'run.perlen_max'),
    ('run.nsp', None),
)

# ---- panel 2: the records ------------------------------------------------
# An explicit layout, like the grid sub-panels, rather than whatever order
# the dataclass happens to declare. One column per SUBJECT: the meteorology
# down the left, the irrigation down the right, each in the order it is
# filled in. Irrigation is the longer story and it is optional, so keeping it
# in one column is what makes the short one readable. ``None`` leaves a cell
# empty rather than pulling the next field up into it.
SURFACE_ROWS = (
    ('surface.meteo_ts', 'surface.irrigation'),
    ('surface.meteo_zones', 'surface.irr_zones'),
    ('surface.out_prefix', 'surface.nfield'),
    (None, 'surface.irr_ts'),
    (None, 'surface.crop_schedule'),
)


# Greyed while their switch is off, like the grid sub-panels -- but SHOWING
# what they hold, not cleared. Panel 1 clears a switched-off field because
# what it holds is derived and would be recomputed anyway (D4); these are
# FILENAMES somebody typed or picked, and emptying them on an unticked box
# would mean typing them again on the next tick.
SURFACE_GATED = {
    'surface.irrigation': ('surface.irr_zones', 'surface.nfield',
                           'surface.irr_ts', 'surface.crop_schedule'),
}

# The ones that name a FILE on this machine, so the panel offers the system
# dialog beside them instead of a name to be typed correctly.
SURFACE_FILES = ('surface.meteo_ts', 'surface.irr_ts', 'surface.crop_schedule')

# `crop_schedule` is a PATTERN, one file per field: picking
# __inputFIELD1_crop_schedule.txt has to store __inputFIELD%d_... or the
# second field would read the first one's schedule.
SURFACE_PATTERNS = ('surface.crop_schedule',)


# Which record belongs under which table. The files were listed in a block of
# their own -- "Are they there?" -- which said nothing about what they are
# for; under the table that reads them, presence and purpose are one answer.
SURFACE_TABLE_FILES = {
    'surface.station': ('surface.meteo_ts',),
    'surface.crop': ('surface.irr_ts', 'surface.crop_schedule'),
}

# MMsurf's own figures are a PLOTTING choice, so they are asked on panel 6
# with the rest of them, not here among the records that feed the run.
SURFACE_ON_PLOTS = ('surface.plot',)


# ---- panel 3: the soil column, and the vegetation beside it -------------
# [soil] carries two different subjects: where a cell gets its SOIL from, and
# where it gets its VEGETATION COVER from. They share a section because they
# share a file, not because they are one question, and shown together the
# second was read as more soil settings.
SOIL_ROWS = (
    ('soil.params', 'soil.zones'),
    ('soil.thickness', None),
)
SOIL_VEG_ROWS = (
    ('soil.veg_layer', 'soil.veg_column'),
)


# A plain field that names an ATTRIBUTE of the layer another field names.
# Not part of a VectorSource -- the vegetation cover is two strings, because
# it is read by the converter rather than wrapped like a zone map -- but the
# same question, so it gets the same list.
#   the column field -> the layer field it belongs to
COLUMN_OF = {
    'soil.veg_column': 'soil.veg_layer',
}

# Named a shapefile, so it is chosen the way every other one is.
SOIL_FILES = ('soil.veg_layer',)

# Named a file in the DATASET, not in the cartography folder -- a different
# starting point for the same dialog.
SOIL_DATASET_FILES = ('soil.params',)

GRID_PERMANENT = ('grid.boundary', 'grid.streams', 'grid.ponds', 'grid.dem',
                  'grid.crs_epsg', 'grid.kind', 'grid.resample')
_LEGACY_ROWS = (
    ('grid.cell_size', 'grid.buffer'),
    ('grid.override.enable',),
    ('grid.override.xllcorner', 'grid.override.yllcorner',
     'grid.override.nrow', 'grid.override.ncol'),
)
GRID_SUBPANEL = {
    'structured': _LEGACY_ROWS,
    'disv': _LEGACY_ROWS,
    # The modeller gives the two numbers that mean something -- the cell size
    # AT the stream and how fast it may grow -- and the corridor follows from
    # them (WP1d, panel 1 item 9). It used to be typed, and a corridor
    # narrower than the grade takes clamped the refinement silently.
    'voronoi': (
        ('grid.voronoi.cell_far', 'grid.buffer'),
        ('grid.voronoi.stream_refine',),
        ('grid.voronoi.cell_near_stream', 'grid.voronoi.grade_ratio'),
        ('grid.voronoi.stream_buffer', 'grid.voronoi.trans_levels'),
        ('grid.voronoi.refine_ponds',),
        ('grid.voronoi.cell_pond',),
    ),
    'quadtree': (
        ('grid.cell_size', 'grid.buffer'),
        ('grid.quadtree.refine_streams',),
        ('grid.quadtree.refine_level',),
        ('grid.quadtree.refine_ponds',),
        ('grid.quadtree.pond_level',),
    ),
}


def grid_fields(kind):
    """Every field on one kind's sub-panel, flattened out of its rows."""
    return tuple(d for row in GRID_SUBPANEL.get(kind, ()) for d in row)


# Shown but never typed: the geometry owns them.
GRID_DERIVED = ('grid.voronoi.stream_buffer', 'grid.voronoi.trans_levels')


def _quadtree_note(get, ponds=False):
    """What the level on screen actually gives, in metres.

    Under the BOX, because it is the sentence the refined size is read off
    and a sentence away from the number it describes is read against the
    wrong one.
    """
    bg = float(get('grid.cell_size', 0.0) or 0.0)
    lv = int(get('grid.quadtree.refine_level', 0) or 0)
    if not ponds:
        return ('GRIDGEN halves a cell per refinement level, so level %d on a '
                '%g m background gives %g m along the streams.'
                % (lv, bg, bg / (2 ** lv) if lv >= 0 else bg))
    pl = int(get('grid.quadtree.pond_level', 0) or 0)
    used = pl or lv
    got = bg / (2 ** used) if used >= 0 else bg
    if pl:
        return ('Level %d on a %g m background gives %g m over the pond '
                'footprints.' % (used, bg, got))
    return ('0 follows the streams, so level %d on a %g m background gives '
            '%g m over the pond footprints.' % (used, bg, got))


# A note drawn UNDER one field, computed from what is on screen.
GRID_NOTES = {
    'grid.quadtree.refine_level': lambda get: _quadtree_note(get),
    'grid.quadtree.pond_level': lambda get: _quadtree_note(get, ponds=True),
}


# Greyed out, and CLEARED, while their controlling switch is off (panel 1 D4).
# A switch that cannot be answered yet, because the layer it acts on has
# not been given. Keyed on the CONFIGURATION field that must be non-blank.

GRID_NEEDS_LAYER = {
    'grid.voronoi.stream_refine': 'grid.streams',
    'grid.quadtree.refine_streams': 'grid.streams',
    'grid.voronoi.refine_ponds': 'grid.ponds',
    'grid.quadtree.refine_ponds': 'grid.ponds',
}

GRID_GATED = {
    'grid.voronoi.stream_refine': ('grid.voronoi.cell_near_stream',
                                   'grid.voronoi.stream_buffer',
                                   'grid.voronoi.grade_ratio',
                                   'grid.voronoi.trans_levels'),
    'grid.quadtree.refine_streams': ('grid.quadtree.refine_level',),
    'grid.voronoi.refine_ponds': ('grid.voronoi.cell_pond',),
    'grid.quadtree.refine_ponds': ('grid.quadtree.pond_level',),
    'grid.override.enable': ('grid.override.xllcorner',
                             'grid.override.yllcorner',
                             'grid.override.nrow', 'grid.override.ncol'),
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
