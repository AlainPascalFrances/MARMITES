# -*- coding: utf-8 -*-
"""The MMsurf ini as a structured, self-describing parameter set.  WP1.7.

`__inputMMsurf.ini` is a positional file: counts followed by one whitespace-
separated line per meteo zone, vegetation type, crop and soil, with the meaning
of each column recorded only in a comment block at the top. That is the same
fragility the MM ini had before Phase 2 -- one missing line silently shifts
every parameter after it.

This module turns it into records that carry their own names, units and
descriptions, so that:

  * the Streamlit Configuration page can show every parameter ALREADY FILLED
    from the case study, with its meaning, instead of a wall of numbers;
  * the parameter list is machine-readable, which is what WP7 needs to decide
    what to calibrate (`Zr` and the kTg_* sourcing factors above all);
  * a malformed file is reported against a named field rather than shifting
    everything silently.

It is deliberately READ-ONLY and streamlit-free: the model still reads the ini
through its own path, and nothing here changes what a run does.

    from mmsurf_config import load_mmsurf_ini, enumerate_parameters
    ms = load_mmsurf_ini('example/LaMata/MMsurf_ws/__inputMMsurf.ini')
    for p in enumerate_parameters(ms):
        print(p.block, p.member, p.name, p.value, p.units, p.description)
"""

import os
from dataclasses import dataclass, field

__all__ = ['MMsurfConfig', 'Param', 'FIELDS', 'load_mmsurf_ini',
           'enumerate_parameters', 'MMsurfError']


class MMsurfError(Exception):
    """Malformed MMsurf ini."""


# (name, units, description) per column, in file order, after the leading name.
# Transcribed from the comment block of __inputMMsurf.ini -- the only place the
# meanings were previously recorded.
METEO_FIELDS = [
    ('phi', 'deg', 'latitude of the meteo station (>0 = northern hemisphere)'),
    ('Lm', 'deg', 'longitude of the station, west of Greenwich'),
    ('Z', 'm', 'altitude of the station above sea level'),
    ('Lz', 'deg', 'longitude of the centre of the local time zone, west of Greenwich'),
    ('FC', 'h', 'time shift for a site not in its nominal time zone'),
    ('DTS', 'h', 'data time shift, for data not acquired at standard clock time'),
    ('z_m', 'm', 'height of the wind-speed measurement'),
    ('z_h', 'm', 'height of the humidity measurement'),
]
VEG_FIELDS = [
    ('h_d', 'm', 'plant height, dry season'),
    ('h_w', 'm', 'plant height, wet season'),
    ('S_w', 'mm', 'canopy storage capacity'),
    ('C_leaf_star', 'm/s', 'maximum leaf conductance'),
    ('LAI_d', 'm2/m2', 'leaf area index, dry season'),
    ('LAI_w', 'm2/m2', 'leaf area index, wet season'),
    ('f_s_vd', '-', 'shelter factor, dry season'),
    ('f_s_vw', '-', 'shelter factor, wet season'),
    ('alfa_vd', '-', 'vegetation albedo, dry season'),
    ('alfa_vw', '-', 'vegetation albedo, wet season'),
    ('J_vd', 'julian day', 'start of the dry season'),
    ('J_vw', 'julian day', 'start of the wet season'),
    ('TRANS_vdw', 'd', 'dry -> wet transition length'),
    ('TRANS_vwd', 'd', 'wet -> dry transition length'),
    ('Zr', 'm', 'MAXIMUM ROOT DEPTH -- the natural source for the WP2 '
                '[et] extdp_source = "veg_zone" extinction depth'),
    ('kTg_min', '-', 'transpiration sourcing factor, minimum (1 >= k_T > 0)'),
    ('kTg_max', '-', 'transpiration sourcing factor, maximum'),
    ('kT_f', '-', 'transpiration sourcing factor f (phi > f > wilting point)'),
    ('kT_s', '-', 'transpiration sourcing factor slope (given as 1/s in the file)'),
]
CROP_FIELDS = [
    ('h_c', 'm', 'maximum plant height'),
    ('S_w_c', 'mm', 'canopy storage capacity'),
    ('C_leaf_star_c', 'm/s', 'maximum leaf conductance'),
    ('LAI_c', 'm2/m2', 'maximum leaf area index'),
    ('f_s_c', '-', 'shelter factor'),
    ('alfa_c', '-', 'albedo'),
    ('Zr_c', 'm', 'maximum root depth'),
    ('kTg_min_c', '-', 'transpiration sourcing factor, minimum'),
    ('kTg_max_c', '-', 'transpiration sourcing factor, maximum'),
    ('kT_f_c', '-', 'transpiration sourcing factor f'),
    ('kT_s_c', '-', 'transpiration sourcing factor slope'),
]
SOIL_FIELDS = [
    ('por', 'm3/m3', 'surface (1 cm) soil porosity'),
    ('fc', 'm3/m3', 'surface (1 cm) soil field capacity'),
    ('alfa_sd', '-', 'soil albedo, dry season'),
    ('alfa_sw', '-', 'soil albedo, wet season'),
    ('J_sd', 'julian day', 'start of the dry season'),
    ('J_sw', 'julian day', 'start of the wet season'),
    ('TRANS_sdw', 'd', 'dry -> wet transition length'),
    ('TRANS_swd', 'd', 'wet -> dry transition length'),
]
FIELDS = {'meteo': METEO_FIELDS, 'veg': VEG_FIELDS,
          'crop': CROP_FIELDS, 'soil': SOIL_FIELDS}

# Where each block lands in the WP0 schema once it migrates.
SCHEMA_DEST = {'meteo': '[mmsurf.meteo]', 'veg': '[mmsurf.veg]',
               'crop': '[mmsurf.crop]', 'soil': '[mmsurf.soil]',
               'openwater': '[mmsurf.openwater]'}


@dataclass
class Param:
    """One parameter, with everything needed to display or calibrate it."""

    block: str          # meteo | veg | crop | soil | openwater
    member: str         # the zone / vegetation / crop / soil name
    name: str
    value: float
    units: str
    description: str

    @property
    def dotted(self):
        return '%s.%s.%s' % (self.block, self.member, self.name)


@dataclass
class MMsurfConfig:
    """The MMsurf ini, parsed into named records."""

    path: str = ''
    meteo: list = field(default_factory=list)   # [(name, {field: value})]
    veg: list = field(default_factory=list)
    crop: list = field(default_factory=list)
    soil: list = field(default_factory=list)
    openwater: dict = field(default_factory=dict)
    nfield: int = 0

    @property
    def counts(self):
        return {'NMETEO': len(self.meteo), 'NVEG': len(self.veg),
                'NCRP': len(self.crop), 'NFIELD': self.nfield,
                'NSOIL': len(self.soil)}


def _read_tokens(path):
    """Reproduce clsUTILITIES.readFile: the first character of line 1 is the
    comment delimiter; on every later line keep what precedes the first one."""
    if not os.path.exists(path):
        raise MMsurfError('MMsurf ini does not exist: %s' % path)
    out = []
    with open(path, encoding='utf-8-sig') as fh:
        first = fh.readline().split()
        dc = first[0] if first else '#'
        for line in fh:
            tok = line.split(dc)[0]
            if tok and not tok.isspace():
                out.append(tok.strip())
    return out


def _named_row(raw, fields, block, index):
    parts = raw.split()
    if len(parts) < 1 + len(fields):
        raise MMsurfError(
            '%s entry %d ("%s"): expected a name followed by %d values '
            '(%s), got %d value(s). The ini is positional, so a short line '
            'shifts everything after it.'
            % (block, index, parts[0] if parts else '?', len(fields),
               ', '.join(f[0] for f in fields), len(parts) - 1))
    name = parts[0]
    vals = {}
    for (fname, _u, _d), tok in zip(fields, parts[1:]):
        try:
            vals[fname] = float(tok)
        except ValueError:
            raise MMsurfError('%s entry %d (%s): %s = %r is not a number'
                              % (block, index, name, fname, tok)) from None
    return name, vals


def load_mmsurf_ini(path):
    """Parse `__inputMMsurf.ini` into an MMsurfConfig."""
    tok = _read_tokens(path)
    cfg = MMsurfConfig(path=str(path))
    i = 0

    def take():
        nonlocal i
        if i >= len(tok):
            raise MMsurfError('MMsurf ini ended prematurely at entry %d' % i)
        v = tok[i]
        i += 1
        return v

    def take_int(what):
        v = take().split()[0]
        try:
            return int(v)
        except ValueError:
            raise MMsurfError('%s: expected an integer, got %r' % (what, v)) from None

    nmeteo = take_int('NMETEO')
    for k in range(nmeteo):
        cfg.meteo.append(_named_row('METEO%d %s' % (k + 1, take()),
                                    METEO_FIELDS, 'meteo', k + 1))
    nveg = take_int('NVEG')
    for k in range(nveg):
        cfg.veg.append(_named_row(take(), VEG_FIELDS, 'veg', k + 1))
    ncrop = take_int('NCRP')
    cfg.nfield = take_int('NFIELD')
    for k in range(ncrop):
        cfg.crop.append(_named_row(take(), CROP_FIELDS, 'crop', k + 1))
    nsoil = take_int('NSOIL')
    for k in range(nsoil):
        cfg.soil.append(_named_row(take(), SOIL_FIELDS, 'soil', k + 1))
    # Open water is optional: alfa_w is commented out in the La Mata file, so
    # MMsurf falls back to its own default.
    if i < len(tok):
        try:
            cfg.openwater = {'alfa_w': float(tok[i].split()[0])}
        except ValueError:
            pass
    return cfg


def enumerate_parameters(cfg):
    """Flatten an MMsurfConfig into Param records, in file order."""
    out = []
    for block, rows in (('meteo', cfg.meteo), ('veg', cfg.veg),
                        ('crop', cfg.crop), ('soil', cfg.soil)):
        for name, vals in rows:
            for fname, units, desc in FIELDS[block]:
                out.append(Param(block, name, fname, vals[fname], units, desc))
    for k, v in cfg.openwater.items():
        out.append(Param('openwater', 'water', k, v, '-', 'open-water albedo'))
    return out


if __name__ == '__main__':
    import sys
    p = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), '..', 'example', 'LaMata',
        'MMsurf_ws', '__inputMMsurf.ini')
    c = load_mmsurf_ini(p)
    print('%s\n  counts: %s' % (c.path, c.counts))
    for prm in enumerate_parameters(c):
        print('  %-26s %-12g %-10s %s'
              % (prm.dotted, prm.value, prm.units, prm.description[:58]))
