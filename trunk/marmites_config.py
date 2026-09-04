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


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 3:
        c = convert_ini_file(sys.argv[1], sys.argv[2])
        print('Converted %s -> %s (run_name=%s)' % (sys.argv[1], sys.argv[2], c.run_name))
    else:
        print('usage: python marmites_config.py <legacy_MM.ini> <out.toml>')
