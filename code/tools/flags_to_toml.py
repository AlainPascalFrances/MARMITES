# -*- coding: utf-8 -*-
"""Translate an old run_lamata_mf6.py command line into a WP0 TOML file.

The flags are gone; this is the bridge for the commands still written down in
notebooks, in the handoff document and in old shell history.

    python code/tools/flags_to_toml.py --nlay 2 --seep drn --strt-heads hi_spinup \
        --steady-means hi_spinup --preproc --postproc

prints the equivalent configuration; add ``-o out.toml`` to write it. Flags that
became machine paths (--libmf6, --ws-root, --gis-ws) are reported separately,
because they belong in code/mm_paths.py, not in a run configuration.
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
if CODE not in sys.path:
    sys.path.insert(0, CODE)

import marmites_config as mcfg  # noqa: E402

# old flag -> (section, key) or a callable(cfg, value)
_SIMPLE = {
    '--mode': ('run', 'mode', str),
    '--relax': ('run', 'relax', float),
    '--nsp': ('run', 'nsp', int),
    '--max-discrepancy': ('run', 'max_discrepancy', float),
    '--grid': ('grid', 'kind', str),
    '--nlay': ('layers', 'nlay', int),
    '--uzf-vks-scale': ('uzf', 'vks_scale', float),
    '--seep': ('seep', 'kind', str),
    '--seep-cond': ('seep', 'cond', float),
    '--sfr-rhk': ('sfr', 'rhk', float),
    '--lak-bedleak': ('lak', 'bedleak', float),
    '--spinup': ('spinup', 'cycles', int),
    '--spinup-tol': ('spinup', 'tol', float),
    '--strt-heads': ('spinup', 'strt_heads', str),
    '--steady-means': ('spinup', 'steady_means', str),
    '--sankey-min-flux': ('postproc', 'sankey_min_flux', float),
    '--map-days': ('postproc', 'map_days', int),
    '--run-tag': ('meta', 'name', str),
    '--standalone': ('run', 'standalone', str),
}
_STORE_TRUE = {
    '--build-only': ('run', 'build_only'),
    '--aggregate': ('layers', 'aggregate'),
    '--sfr': ('sfr', 'enable'),
    '--postproc': ('postproc', 'enable'),
    '--preproc': ('postproc', 'preproc'),
    '--postproc-only': ('postproc', 'only'),
    '--sankey-obs-years': ('postproc', 'sankey_obs_years'),
    '--allow-bad-budget': ('run', 'allow_bad_budget'),
    '--probe': ('run', 'probe'),
}
_STORE_FALSE = {
    '--aggregated': ('run', 'daily'),
    '--no-ats': ('run', 'ats'),
    '--no-sankey-full': ('postproc', 'sankey_full'),
}
_OPTIONAL_ARG = {          # nargs='?' with a const
    '--save-strt': ('spinup', 'save_strt', 'hi_spinup'),
    '--save-means': ('spinup', 'save_means', 'hi_spinup'),
    '--lak': ('lak', 'source', 'GIS/lm_ponds.shp'),
}
_MACHINE = {'--libmf6': 'paths.libmf6 (or mm_paths.LIBMF6)',
            '--ws-root': 'mm_paths.WS_ROOT',
            '--ws': 'paths.ws',
            '--gis-ws': 'paths.gis_ws (or mm_paths.GIS)'}


def convert(argv):
    cfg = mcfg.RunConfig.from_dict({})
    machine, notes = [], []
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok in _SIMPLE:
            sec, key, cast = _SIMPLE[tok]
            i += 1
            val = cast(argv[i])
            # sfr.rhk and friends are ParamSource fields, not bare numbers
            if isinstance(getattr(getattr(cfg, sec), key), mcfg.ParamSource):
                val = mcfg.ParamSource(value=float(val))
            setattr(getattr(cfg, sec), key, val)
        elif tok in _STORE_TRUE:
            sec, key = _STORE_TRUE[tok]
            setattr(getattr(cfg, sec), key, True)
        elif tok in _STORE_FALSE:
            sec, key = _STORE_FALSE[tok]
            setattr(getattr(cfg, sec), key, False)
        elif tok in _OPTIONAL_ARG:
            sec, key, const = _OPTIONAL_ARG[tok]
            nxt = argv[i + 1] if i + 1 < len(argv) else ''
            val = nxt if (nxt and not nxt.startswith('--')) else const
            if val is nxt and nxt:
                i += 1
            setattr(getattr(cfg, sec), key, val)
            if tok == '--lak':
                cfg.lak.enable = True
        elif tok == '--strt-dem':
            cfg.spinup.strt_dem = [float(argv[i + 1]), float(argv[i + 2])]
            i += 2
        elif tok == '--daily':
            cfg.run.daily = True
        elif tok in _MACHINE:
            nxt = argv[i + 1] if i + 1 < len(argv) else ''
            if nxt and not nxt.startswith('--'):
                i += 1
            machine.append('%-14s -> %s   (%s)' % (tok, _MACHINE[tok], nxt or 'default'))
        elif tok.endswith('run_lamata_mf6.py') or tok in ('python', 'python.exe'):
            pass
        else:
            notes.append('unrecognised token ignored: %r' % tok)
        i += 1
    if cfg.postproc.only:
        cfg.postproc.enable = True
    cfg.validate()
    return cfg, machine, notes


def main():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument('-o', '--out', default=None)
    ap.add_argument('-h', '--help', action='store_true')
    known, rest = ap.parse_known_args()
    if known.help or not rest:
        print(__doc__)
        return 0
    cfg, machine, notes = convert(rest)
    if known.out:
        cfg.write_toml(known.out)
        print('written: %s' % known.out)
    else:
        import tempfile
        tmp = os.path.join(tempfile.gettempdir(), '_flags_to_toml.toml')
        cfg.write_toml(tmp)
        with open(tmp, encoding='utf-8') as fh:
            sys.stdout.write(fh.read())
        os.remove(tmp)
    if machine:
        print('\n# these flags became MACHINE paths -- set them in code/mm_paths.py:')
        for m in machine:
            print('#   %s' % m)
    for n in notes:
        print('# %s' % n, file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
