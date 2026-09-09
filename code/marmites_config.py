# -*- coding: utf-8 -*-
"""MARMITES run configuration (Phase 2).

Replaces the fragile positional ``l += 1`` parsing of the MM ini file
(where a single missing line silently shifted every subsequent
parameter) with a validated TOML schema backed by a dataclass.

Two entry points:
  * ``load_mm_config(path)`` reads a TOML file into an ``MMConfig``;
  * ``legacy_ini_to_config(lines)`` / ``convert_ini_file(ini_path, toml_path)``
    read the old sequential ini and produce the same ``MMConfig`` / a TOML
    file, so existing datasets migrate without manual editing.

The MODFLOW ini is intentionally *not* modelled here: its parser
(ppMODFLOW_flopy_v3) is rewritten for the MODFLOW 6 packages in Phase 3,
so a TOML schema for the NWT-era fields would be thrown away.
"""

import dataclasses
import os
from dataclasses import dataclass, field
# (typing generics kept minimal for loader compatibility)


def _load_tomllib():
    """tomllib is stdlib on the target (Python 3.12); fall back to the
    third-party ``tomli`` on older interpreters used for testing."""
    try:
        import tomllib
        return tomllib
    except ModuleNotFoundError:
        try:
            import tomli
            return tomli
        except ModuleNotFoundError as exc:  # pragma: no cover
            raise ConfigError(
                'Reading TOML requires Python >= 3.11 (tomllib) or the '
                '`tomli` package on older interpreters.') from exc

__all__ = ['MMConfig', 'ConfigError', 'load_mm_config', 'legacy_ini_to_config',
           'convert_ini_file', 'dump_toml']


class ConfigError(Exception):
    """Invalid or inconsistent MARMITES configuration."""


# Order of the run-control fields as they appear in the legacy MM ini.
# (name, kind) where kind is one of: str, int, float, bool01
# bool01 means "1/0 in the ini -> True/False". Conditional irrigation
# lines are handled specially in legacy_ini_to_config.
_LEGACY_ORDER = [
    ('run_name', 'str'),
    ('verbose', 'int'),
    ('plt_out', 'bool01'),
    ('plt_freq', 'int'),
    ('nrangeMM', 'int'),
    ('nrangeMF', 'int'),
    ('ctrsMM', 'bool01'),
    ('ctrsMF', 'bool01'),
    ('ntick', 'int'),
    ('animation', 'int'),
    ('animation_freq', 'int'),
    ('plt_out_obs', 'bool01'),
    ('WBsankey_yn', 'bool01'),
    ('plt_WB_unit', 'str'),
    ('iniMonthHydroYear', 'int'),
    ('maxYearsTickTrimester', 'int'),
    ('maxYearsTickSemester', 'int'),
    ('plt_input', 'bool01'),
    ('MMsurf_yn', 'bool01'),
    ('MMsurf_plot', 'bool01'),
    ('MMsoil_yn', 'int'),   # -1/0/1, not a plain bool
    ('MF_yn', 'bool01'),
    ('MF_lastrun', 'bool01'),
    ('convcrit', 'float'),
    ('convcritmax', 'float'),
    ('ccnum', 'int'),
    ('MMsurf_ws', 'str'),
    ('inputFile_PAR_fn', 'str'),
    ('inputFile_TS_fn', 'str'),
    ('irr_yn', 'bool01'),   # conditional block follows
    ('outputFILE_fn', 'str'),
    ('outMMsurf_fn', 'str'),
    ('MF_ws', 'str'),
    ('MF_ini_fn', 'str'),
    ('xllcorner', 'float'),
    ('yllcorner', 'float'),
    ('gridMETEO_fn', 'str'),
    ('gridSOIL_fn', 'str'),
    ('gridSOILthick_fn', 'str'),
    ('gridSsurfhmax_fn', 'str'),
    ('gridSsurfw_fn', 'str'),
    ('SOILparam_fn', 'str'),
    ('inputObs_fn', 'str'),
    ('inputObsHEADS_fn', 'str'),
    ('inputObsSM_fn', 'str'),
    ('inputObsRo_fn', 'str'),
    ('rmseHEADSmax', 'float'),
    ('rmseSMmax', 'float'),
    ('chunks', 'int'),
]

# fields made obsolete by the Phase-1 removal of the MM-MF Picard loop;
# still accepted (for legacy datasets) but ignored by the driver.
_DEPRECATED = {'MF_yn', 'MF_lastrun', 'convcrit', 'convcritmax', 'ccnum'}


@dataclass
class MMConfig:
    """Validated MARMITES run configuration (run-control block)."""

    run_name: str
    # output / plotting
    verbose: int = 0
    plt_out: bool = True
    plt_freq: int = 30
    nrangeMM: int = 5
    nrangeMF: int = 5
    ctrsMM: bool = False
    ctrsMF: bool = False
    ntick: int = 10
    animation: int = 0
    animation_freq: int = 12
    plt_out_obs: bool = True
    WBsankey_yn: bool = True
    plt_WB_unit: str = 'year'
    iniMonthHydroYear: int = 10
    maxYearsTickTrimester: int = 5
    maxYearsTickSemester: int = 10
    plt_input: bool = True
    # run switches
    MMsurf_yn: bool = False
    MMsurf_plot: bool = False
    MMsoil_yn: int = 1          # -1 (single MMsoil calib run) / 0 / 1
    plt_input_maps: bool = True
    # surface model inputs
    MMsurf_ws: str = 'MMsurf_ws'
    inputFile_PAR_fn: str = ''
    inputFile_TS_fn: str = ''
    irr_yn: bool = False
    inputFile_TSirr_fn: str = ''
    gridIRR_fn: str = ''
    outputFILE_fn: str = ''
    outMMsurf_fn: str = ''
    # MODFLOW workspace
    MF_ws: str = 'MF_ws'
    MF_ini_fn: str = ''
    xllcorner: float = 0.0
    yllcorner: float = 0.0
    # spatial input rasters
    gridMETEO_fn: str = ''
    gridSOIL_fn: str = ''
    gridSOILthick_fn: str = ''
    gridSsurfhmax_fn: str = ''
    gridSsurfw_fn: str = ''
    SOILparam_fn: str = ''
    # observations / calibration
    inputObs_fn: str = ''
    inputObsHEADS_fn: str = ''
    inputObsSM_fn: str = ''
    inputObsRo_fn: str = ''
    rmseHEADSmax: float = 1.0
    rmseSMmax: float = 1.0
    chunks: int = 0
    # deprecated (Picard loop removed in Phase 1) -- kept for legacy files
    deprecated: dict = field(default_factory=dict)

    def validate(self) -> None:
        errs = []
        if not self.run_name:
            errs.append('run_name must not be empty')
        if not (1 <= self.iniMonthHydroYear <= 12):
            errs.append('iniMonthHydroYear must be between 1 and 12')
        if self.plt_WB_unit not in ('year', 'day'):
            errs.append("plt_WB_unit must be 'year' or 'day'")
        if self.MMsoil_yn not in (-1, 0, 1):
            errs.append('MMsoil_yn must be -1, 0 or 1')
        if self.chunks not in (0, 1):
            errs.append('chunks must be 0 or 1')
        if self.maxYearsTickTrimester > self.maxYearsTickSemester:
            errs.append('maxYearsTickTrimester must be <= maxYearsTickSemester')
        if self.irr_yn and not self.gridIRR_fn:
            errs.append('irr_yn is set but gridIRR_fn is empty')
        if errs:
            raise ConfigError('Invalid MARMITES configuration:\n  - ' + '\n  - '.join(errs))

    @classmethod
    def from_dict(cls, d: dict) -> 'MMConfig':
        known = {f.name for f in dataclasses.fields(cls)}
        deprecated = {k: d[k] for k in _DEPRECATED if k in d}
        kwargs = {k: v for k, v in d.items() if k in known and k != 'deprecated'}
        cfg = cls(**kwargs)
        cfg.deprecated.update(deprecated)
        cfg.validate()
        return cfg


def _coerce(raw: str, kind: str):
    raw = raw.strip()
    if kind == 'str':
        return raw
    if kind == 'int':
        return int(raw)
    if kind == 'float':
        return float(raw)
    if kind == 'bool01':
        return int(raw) == 1
    raise ConfigError(f'unknown field kind {kind!r}')


def legacy_ini_to_config(values):
    """Build an MMConfig from the sequential value list of a legacy MM ini.

    ``values`` is the list of content tokens as returned by
    clsUTILITIES.readFile (comment lines already stripped).
    """
    d = {}
    it = iter(range(len(values)))
    idx = 0

    def take(kind):
        nonlocal idx
        if idx >= len(values):
            raise ConfigError('legacy ini ended prematurely at field index %d' % idx)
        v = _coerce(values[idx], kind)
        idx += 1
        return v

    for name, kind in _LEGACY_ORDER:
        d[name] = take(kind)
        if name == 'irr_yn':
            # conditional irrigation block: 2 extra lines when irr_yn == 1,
            # else 3 skipped lines in the legacy layout (l += 3)
            if d['irr_yn']:
                d['inputFile_TSirr_fn'] = take('str')
                d['gridIRR_fn'] = take('str')
            else:
                idx += 3  # legacy: 3 placeholder lines skipped
    _ = it  # silence unused
    return MMConfig.from_dict(d)


def convert_ini_file(ini_path, toml_path, delimiter='#'):
    """Read a legacy MM ini file and write an equivalent TOML file."""
    values = _read_legacy_ini(ini_path, delimiter)
    cfg = legacy_ini_to_config(values)
    dump_toml(cfg, toml_path)
    return cfg


def _read_legacy_ini(path, delimiter='#'):
    """Reproduce clsUTILITIES.readFile: first char of line 1 is the
    comment delimiter; keep the content before the first delimiter on
    each subsequent non-blank line."""
    if not os.path.exists(path):
        raise ConfigError(f"ini file does not exist: {path}")
    out = []
    with open(path, encoding='utf-8-sig') as f:
        first = f.readline().split()
        dc = first[0] if first else delimiter
        for line in f:
            tok = line.split(dc)[0]
            if tok and tok != '\n' and not tok.isspace():
                out.append(tok.strip())
    return out


def load_mm_config(path):
    """Load and validate a MARMITES TOML configuration file."""
    if not os.path.exists(path):
        raise ConfigError(f"config file does not exist: {path}")
    tomllib = _load_tomllib()
    with open(path, 'rb') as f:
        data = tomllib.load(f)
    # accept either a flat table or a [marmites] section
    if 'marmites' in data and isinstance(data['marmites'], dict):
        data = data['marmites']
    return MMConfig.from_dict(data)


def dump_toml(cfg, path):
    """Write an MMConfig to a TOML file (stdlib has no writer, so emit
    a minimal, well-formed table by hand)."""
    def fmt(v):
        if isinstance(v, bool):
            return 'true' if v else 'false'
        if isinstance(v, (int, float)):
            return repr(v)
        return '"%s"' % str(v).replace('\\', '\\\\').replace('"', '\\"')

    lines = ['# MARMITES run configuration (generated from legacy ini)', '[marmites]']
    for f in dataclasses.fields(cfg):
        if f.name == 'deprecated':
            continue
        lines.append('%s = %s' % (f.name, fmt(getattr(cfg, f.name))))
    if cfg.deprecated:
        lines += ['', '# obsolete since Phase 1 (MM-MF Picard loop removed); ignored',
                  '[marmites.deprecated]']
        for k, v in cfg.deprecated.items():
            lines.append('%s = %s' % (k, fmt(v)))
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')



# ===========================================================================
# WP0 -- the RUN configuration: the ~40 argparse flags become schema keys.
#
# Design rules, all enforced below:
#   * every default equals TODAY'S flag default, so an empty file reproduces
#     current behaviour and WP0 stays a pure refactor (verified byte-identical
#     MF6 inputs).  WP1c is what flips grid.kind to "voronoi".
#   * UNKNOWN KEYS RAISE -- a typo must never fall through to a default.
#   * config_hash() stamps the run folder and every cache, so a stale grid /
#     layer set / spin-up state can never silently override the configuration
#     (the CdL grid-design cache did exactly that once).
# ===========================================================================


@dataclass
class Meta:
    config_version: int = 1
    name: str = ''                 # -> out_<stamp>_<name>; blank = <nlay>lay_<mode>
    description: str = ''


@dataclass
class Paths:
    case: str = 'LaMata'           # -> example/<case>/ (Tier A inputs)
    ws: str = ''                   # blank -> WS_ROOT/MF6_ws[_disv]
    gis_ws: str = ''               # blank -> mm_paths.GIS (figures only)
    nwt_reference: str = ''        # blank -> mm_paths.NWT_REF
    libmf6: str = ''               # blank -> mm_paths.LIBMF6


@dataclass
class Run:
    mode: str = 'lagged'           # lagged | iterative        (--mode)
    relax: float = 0.6             #                           (--relax)
    nsp: int = 0                   # 0 = all stress periods    (--nsp)
    daily: bool = True             # false = aggregated        (--aggregated)
    ats: bool = True               #                           (--no-ats)
    build_only: bool = False       #                           (--build-only)
    standalone: str = ''           # mf6.exe path              (--standalone)
    probe: bool = False            #                           (--probe)
    max_discrepancy: float = 1.0   #                           (--max-discrepancy)
    allow_bad_budget: bool = False  #                          (--allow-bad-budget)


@dataclass
class GridVoronoi:
    """WP1c. Ignored while grid.kind is 'structured' or 'disv'."""

    cell_far: float = 100.0
    cell_near_stream: float = 40.0
    stream_buffer: float = 60.0
    stream_refine: bool = True
    trans_levels: list = field(default_factory=lambda: [10.0, 20.0, 40.0, 70.0])
    seed_ponds: bool = True


@dataclass
class Grid:
    # 'structured' is today's DIS grid and stays the DEFAULT until WP1c lands;
    # 'dis' is accepted as an alias.
    kind: str = 'structured'
    rebuild: bool = False
    voronoi: GridVoronoi = field(default_factory=GridVoronoi)


@dataclass
class Layers:
    nlay: int = 6                  # 2 reads _2s1L.ini directly (--nlay)
    aggregate: bool = False        # 6->2 derivation, comparison (--aggregate)


@dataclass
class Uzf:
    vks_scale: float = 1.0         #                           (--uzf-vks-scale)


@dataclass
class Seep:
    kind: str = 'uzf'              # uzf | drn                 (--seep)
    cond: float = 10000.0          #                           (--seep-cond)


@dataclass
class Et:
    """WP2. Defaults OFF so WP0 changes no behaviour."""

    uzf_et: bool = False           # UZF SIMULATE_ET
    unsat_form: str = 'etwc'       # etwc | etae
    gwet_in_mf: bool = False       # MUST stay false: ETg comes from MM
    extdp_source: str = 'uniform'  # uniform | veg_zone | raster
    extdp_default: float = 2.0     # m
    extwc_source: str = 'thtr'


@dataclass
class ParamSource:
    """A parameter that declares WHERE its value comes from (cookbook 3a).

    Exactly one producer may be set:
        value    -- one number for the whole network
        column   -- an attribute column of the source shapefile
        raster   -- sampled from a raster
        drainage -- computed as w = a * A**b from contributing drainage area
    Resolved by the CONVERTER (WP1), never by the model.
    """

    value: float = None
    column: str = ''
    raster: str = ''
    drainage: dict = field(default_factory=dict)

    def _set(self):
        return [n for n in ('value', 'column', 'raster', 'drainage')
                if getattr(self, n) not in (None, '', {})]

    def producer(self):
        s = self._set()
        return s[0] if len(s) == 1 else None

    def validate(self, what):
        s = self._set()
        if len(s) != 1:
            raise ConfigError(
                '%s: exactly one producer must be set (value | column | raster | '
                'drainage), got %s' % (what, s or 'none'))
        if s[0] == 'drainage':
            miss = {'a', 'b'} - set(self.drainage)
            if miss:
                raise ConfigError('%s.drainage needs a and b (w = a*A**b), missing %s'
                                  % (what, sorted(miss)))

    @classmethod
    def from_value(cls, v):
        """Accept a bare number as shorthand for a value producer."""
        if isinstance(v, dict):
            return _build(cls, v, 'parameter')
        return cls(value=float(v))


@dataclass
class Sfr:
    enable: bool = False           #                           (--sfr)
    source: str = 'inputSTREAM.csv'
    min_slope: float = 1e-4
    monotonic_bed: bool = True     # SFRmaker rule (Leaf et al. 2021)
    width: ParamSource = field(
        default_factory=lambda: ParamSource(drainage={'a': 0.5, 'b': 0.35}))
    manning: ParamSource = field(default_factory=lambda: ParamSource(value=0.035))
    rhk: ParamSource = field(default_factory=lambda: ParamSource(value=0.1))
    rbth: ParamSource = field(default_factory=lambda: ParamSource(value=0.5))


@dataclass
class Lak:
    enable: bool = False           #                           (--lak)
    source: str = 'inputPONDS.csv'
    bedleak: float = 1e-3          # 1/d                       (--lak-bedleak)
    surfdep: float = 0.05          # m
    maxiter: int = 200             # LAK Newton cap (CdL)
    stagechg: float = 1e-4         # m (CdL)


@dataclass
class Crr:
    """WP5. Disabled by default; enable=false is bit-identical to WP4."""

    enable: bool = False
    beta: float = 1.0              # Daoud Eq. 23; calibrated 0.8-1.0
    sinks: str = 'evaporate'
    dem: str = 'inputDEMfill.asc'


@dataclass
class Spinup:
    cycles: int = 1                #                           (--spinup)
    tol: float = 0.05              #                           (--spinup-tol)
    strt_heads: str = ''           #                           (--strt-heads)
    steady_means: str = ''         #                           (--steady-means)
    save_strt: str = ''            #                           (--save-strt)
    save_means: str = ''           #                           (--save-means)
    strt_dem: list = field(default_factory=list)   # [a, b]     (--strt-dem)


@dataclass
class Postproc:
    enable: bool = False           #                           (--postproc)
    preproc: bool = False          #                           (--preproc)
    only: bool = False             #                           (--postproc-only)
    sankey_min_flux: float = 0.05  #                           (--sankey-min-flux)
    sankey_full: bool = True       #                           (--no-sankey-full)
    sankey_obs_years: bool = False  #                          (--sankey-obs-years)
    map_days: int = 6              #                           (--map-days)


@dataclass
class Pest:
    use_pest_params: bool = False
    params_npz: str = 'pest_optimised.npz'


@dataclass
class Ui:
    """WP1b. Read only by code/app/, never by the model."""

    execution: str = 'local'       # local | server
    runs_dir: str = ''             # blank -> WS_ROOT/runs
    port: int = 8501
    address: str = 'localhost'
    poll_secs: int = 3


_SECTIONS = {
    'meta': Meta, 'paths': Paths, 'run': Run, 'grid': Grid, 'layers': Layers,
    'uzf': Uzf, 'seep': Seep, 'et': Et, 'sfr': Sfr, 'lak': Lak, 'crr': Crr,
    'spinup': Spinup, 'postproc': Postproc, 'pest': Pest, 'ui': Ui,
}

GRID_KINDS = ('structured', 'disv', 'voronoi', 'quadtree')
_GRID_ALIAS = {'dis': 'structured'}
_NOT_YET = ('voronoi', 'quadtree')       # producers arrive with WP1c


def _build(cls, data, where):
    """Instantiate a section dataclass, RAISING on any unknown key."""
    if not isinstance(data, dict):
        raise ConfigError('%s: expected a table, got %s' % (where, type(data).__name__))
    fields_ = {f.name: f for f in dataclasses.fields(cls)}
    unknown = sorted(set(data) - set(fields_))
    if unknown:
        raise ConfigError(
            '%s: unknown key(s) %s -- valid keys are %s'
            % (where, ', '.join(repr(u) for u in unknown), ', '.join(sorted(fields_))))
    kwargs = {}
    for k, v in data.items():
        ftype = fields_[k].type
        if ftype in (GridVoronoi, 'GridVoronoi'):
            v = _build(GridVoronoi, v, '%s.%s' % (where, k))
        elif ftype in (ParamSource, 'ParamSource'):
            v = ParamSource.from_value(v)
        kwargs[k] = v
    return cls(**kwargs)


@dataclass
class RunConfig:
    """The whole run: machine-independent settings for one model configuration."""

    meta: Meta = field(default_factory=Meta)
    paths: Paths = field(default_factory=Paths)
    run: Run = field(default_factory=Run)
    grid: Grid = field(default_factory=Grid)
    layers: Layers = field(default_factory=Layers)
    uzf: Uzf = field(default_factory=Uzf)
    seep: Seep = field(default_factory=Seep)
    et: Et = field(default_factory=Et)
    sfr: Sfr = field(default_factory=Sfr)
    lak: Lak = field(default_factory=Lak)
    crr: Crr = field(default_factory=Crr)
    spinup: Spinup = field(default_factory=Spinup)
    postproc: Postproc = field(default_factory=Postproc)
    pest: Pest = field(default_factory=Pest)
    ui: Ui = field(default_factory=Ui)
    source_path: str = ''

    @classmethod
    def from_dict(cls, data, source_path=''):
        if not isinstance(data, dict):
            raise ConfigError('configuration must be a table')
        unknown = sorted(set(data) - set(_SECTIONS))
        if unknown:
            raise ConfigError(
                'unknown section(s) %s -- valid sections are %s'
                % (', '.join(repr(u) for u in unknown), ', '.join(sorted(_SECTIONS))))
        kwargs = {name: _build(_SECTIONS[name], data[name], name) for name in data}
        cfg = cls(source_path=str(source_path), **kwargs)
        cfg.validate()
        return cfg

    def validate(self):
        errs = []
        if self.meta.config_version != 1:
            errs.append('meta.config_version %r is not supported (expected 1)'
                        % self.meta.config_version)
        if self.run.mode not in ('lagged', 'iterative'):
            errs.append("run.mode must be 'lagged' or 'iterative'")
        if not (0.0 < self.run.relax <= 1.0):
            errs.append('run.relax must be in (0, 1]')
        if self.run.nsp < 0:
            errs.append('run.nsp must be >= 0 (0 = all stress periods)')
        if self.grid_kind not in GRID_KINDS:
            errs.append('grid.kind must be one of %s' % ', '.join(GRID_KINDS))
        if self.layers.nlay not in (2, 6):
            errs.append('layers.nlay must be 2 or 6')
        if self.seep.kind not in ('uzf', 'drn'):
            errs.append("seep.kind must be 'uzf' or 'drn'")
        if self.seep.kind == 'drn' and self.seep.cond <= 0:
            errs.append('seep.cond must be > 0 (a seepage face must be free-draining)')
        if self.et.unsat_form not in ('etwc', 'etae'):
            errs.append("et.unsat_form must be 'etwc' or 'etae'")
        if self.et.gwet_in_mf:
            errs.append('et.gwet_in_mf must stay false: ETg is computed by MM and '
                        'applied as a WEL sink, so MODFLOW must not remove it too '
                        '(cookbook WP2, the GWET guard)')
        if self.et.extdp_source not in ('uniform', 'veg_zone', 'raster'):
            errs.append("et.extdp_source must be 'uniform', 'veg_zone' or 'raster'")
        if self.spinup.cycles < 1:
            errs.append('spinup.cycles must be >= 1')
        if self.spinup.strt_dem and len(self.spinup.strt_dem) != 2:
            errs.append('spinup.strt_dem must be [] or [a, b]')
        if not (0.0 < self.crr.beta <= 1.0):
            errs.append('crr.beta must be in (0, 1]')
        if self.ui.execution not in ('local', 'server'):
            errs.append("ui.execution must be 'local' or 'server'")
        for name in ('width', 'manning', 'rhk', 'rbth'):
            try:
                getattr(self.sfr, name).validate('sfr.%s' % name)
            except ConfigError as exc:
                errs.append(str(exc))
        if errs:
            raise ConfigError('Invalid MARMITES run configuration'
                              + (' (%s)' % self.source_path if self.source_path else '')
                              + ':\n  - ' + '\n  - '.join(errs))

    def apply_overrides(self, assignments, echo=True):
        """Apply ``section.key=value`` strings (the CLI --set switch).

        Every override is echoed: a setting that changes a run without
        appearing anywhere is how a multi-hour run gets wasted.
        """
        applied = []
        for item in assignments or ():
            if '=' not in item:
                raise ConfigError('--set expects section.key=value, got %r' % item)
            dotted, raw = item.split('=', 1)
            parts = dotted.strip().split('.')
            if len(parts) < 2:
                raise ConfigError('--set expects section.key=value, got %r' % item)
            obj = self
            for p in parts[:-1]:
                if not hasattr(obj, p):
                    raise ConfigError('--set: no such section %r in %r' % (p, dotted))
                obj = getattr(obj, p)
            leaf = parts[-1]
            if not hasattr(obj, leaf):
                raise ConfigError('--set: no such key %r in %r' % (leaf, dotted))
            old = getattr(obj, leaf)
            setattr(obj, leaf, _coerce_like(old, raw.strip(), dotted))
            applied.append((dotted, old, getattr(obj, leaf)))
        self.validate()
        if echo and applied:
            print('config overrides (--set):')
            for dotted, old, new in applied:
                print('   %-28s %r -> %r' % (dotted, old, new))
        return applied

    def to_dict(self):
        return {name: dataclasses.asdict(getattr(self, name)) for name in _SECTIONS}

    def config_hash(self):
        """Stable digest of the RESOLVED configuration.

        Stamped into the run folder and carried by every cache (grid, layer set,
        spin-up state) so a mismatch stops the run instead of silently
        overriding what was asked for.
        """
        import hashlib
        import json
        payload = json.dumps(self.to_dict(), sort_keys=True, default=str)
        return hashlib.sha256(payload.encode('utf-8')).hexdigest()[:16]

    # Keys that decide whether saved state (spin-up heads, steady means, and in
    # WP1c a cached grid or layer set) is still valid. Deliberately NARROW: the
    # full config_hash would invalidate the state on an unrelated change such as
    # postproc.map_days, and a guard that cries wolf gets switched off.
    STATE_SCOPE = ('paths.case', 'grid.kind', 'layers.nlay', 'layers.aggregate')

    def state_scope(self):
        out = {}
        for dotted in self.STATE_SCOPE:
            sec, key = dotted.split('.')
            v = self.grid_kind if dotted == 'grid.kind' else getattr(
                getattr(self, sec), key)
            out[dotted] = v
        return out

    def state_hash(self):
        """Digest of the keys that make saved state valid or stale."""
        import hashlib
        import json
        payload = json.dumps(self.state_scope(), sort_keys=True, default=str)
        return hashlib.sha256(payload.encode('utf-8')).hexdigest()[:16]

    def write_toml(self, path):
        """Write the resolved configuration (after --set) for provenance."""
        lines = ['# MARMITES resolved run configuration',
                 '# config_hash = %s' % self.config_hash(), '']
        for name in _SECTIONS:
            sec = getattr(self, name)
            nested = []
            lines.append('[%s]' % name)
            for f in dataclasses.fields(sec):
                v = getattr(sec, f.name)
                if dataclasses.is_dataclass(v):
                    nested.append((f.name, v))
                else:
                    lines.append('%s = %s' % (f.name, _toml_value(v)))
            for sub, v in nested:
                lines.append('')
                lines.append('[%s.%s]' % (name, sub))
                for f in dataclasses.fields(v):
                    val = getattr(v, f.name)
                    if val in (None, '', {}, []):
                        continue
                    lines.append('%s = %s' % (f.name, _toml_value(val)))
            lines.append('')
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write('\n'.join(lines))
        return path

    @property
    def grid_kind(self):
        """Normalised grid kind ('dis' is accepted as 'structured')."""
        return _GRID_ALIAS.get(self.grid.kind, self.grid.kind)

    @property
    def run_tag(self):
        return self.meta.name or '%dlay_%s' % (self.layers.nlay, self.run.mode)

    def require_implemented_grid(self):
        """Fail fast and clearly on a grid producer that WP1c has not built yet."""
        if self.grid_kind in _NOT_YET:
            raise ConfigError(
                "grid.kind = %r is not implemented yet: the %s producer arrives "
                "with WP1c. Use 'structured' (today's DIS grid) or 'disv'."
                % (self.grid.kind, self.grid_kind))


def _toml_value(v):
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if v is None:
        return '""'
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, (list, tuple)):
        return '[%s]' % ', '.join(_toml_value(x) for x in v)
    if isinstance(v, dict):
        return '{%s}' % ', '.join('%s = %s' % (k, _toml_value(x)) for k, x in v.items())
    return '"%s"' % str(v).replace('\\', '\\\\').replace('"', '\\"')


def _coerce_like(old, raw, dotted):
    """Coerce a --set string to the type of the value it replaces."""
    if isinstance(old, bool):
        low = raw.lower()
        if low in ('true', '1', 'yes', 'on'):
            return True
        if low in ('false', '0', 'no', 'off'):
            return False
        raise ConfigError('--set %s: expected a boolean, got %r' % (dotted, raw))
    if isinstance(old, int) and not isinstance(old, bool):
        try:
            return int(raw)
        except ValueError:
            raise ConfigError('--set %s: expected an integer, got %r'
                              % (dotted, raw)) from None
    if isinstance(old, float):
        try:
            return float(raw)
        except ValueError:
            raise ConfigError('--set %s: expected a number, got %r'
                              % (dotted, raw)) from None
    if isinstance(old, list):
        raw = raw.strip().strip('[]')
        if not raw:
            return []
        return [float(x) for x in raw.replace(',', ' ').split()]
    return raw


def load_run_config(path):
    """Load and validate a WP0 run-configuration TOML file."""
    if not os.path.exists(path):
        raise ConfigError('config file does not exist: %s' % path)
    tomllib = _load_tomllib()
    with open(path, 'rb') as f:
        data = tomllib.load(f)
    return RunConfig.from_dict(data, source_path=path)


__all__ += ['RunConfig', 'load_run_config', 'GRID_KINDS', 'ParamSource',
            'Meta', 'Paths', 'Run', 'Grid', 'GridVoronoi', 'Layers', 'Uzf',
            'Seep', 'Et', 'Sfr', 'Lak', 'Crr', 'Spinup', 'Postproc', 'Pest', 'Ui']

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 3:
        c = convert_ini_file(sys.argv[1], sys.argv[2])
        print('Converted %s -> %s (run_name=%s)' % (sys.argv[1], sys.argv[2], c.run_name))
    else:
        print('usage: python marmites_config.py <legacy_MM.ini> <out.toml>')
