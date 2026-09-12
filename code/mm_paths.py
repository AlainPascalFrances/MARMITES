# -*- coding: utf-8 -*-
"""The only machine-specific paths (WP0, step 0.2).

Every other module imports its roots from here, so moving the project to a new
machine means setting them ONCE -- on PANEL 0 of the front-end, which writes
``code/configs/paths.local.toml`` -- and nothing else. This is the MARMITES
equivalent of ``CdL/code/config.py`` in the MF6models repository.

    REPO         the MARMITES checkout (self-locating; never set)
    EXAMPLE_ROOT holds one folder per case: the Tier-A inputs a run reads
    DATA_ROOT    raw geospatial sources, OUTSIDE the repo (Tier B)
    GIS          the shapefiles, read by the CONVERTER only
    WS_ROOT      root of ALL run output, OUTSIDE the repo
    NWT_REF      the legacy MODFLOW-NWT _h5_MM.h5 used by the comparison figures
    MODFLOW_DIR  folder holding the MODFLOW 6 / PEST++ executables
    PYTHON_EXE   the conda-env python the PEST++ workers run forward_run with

``DATASET`` resolves to ``<EXAMPLE_ROOT>/<case>/`` and is the only place MM and
MF read input from; nothing here writes into the repo.

Run ``python code/mm_paths.py`` to print what resolved on this machine and
which roots are missing.
"""

import os
from pathlib import Path

__all__ = ['REPO', 'DATA_ROOT', 'WS_ROOT', 'NWT_REF', 'MODFLOW_DIR', 'GIS',
           'LIBMF6', 'MF6_EXE', 'PESTPP_IES', 'PYTHON_EXE', 'TRIANGLE_EXE',
           'GRIDGEN_EXE', 'EXAMPLE_ROOT', 'SETTINGS', 'SETTABLE',
           'dataset_dir', 'resolve_input', 'report_paths', 'LEGACY_ALIASES',
           'read_settings', 'save_settings', 'reload_paths', 'source_of',
           'settings_label']

REPO = Path(__file__).resolve().parents[1]

# ==== where the machine-specific roots come from ============================
# Three sources, in this order:
#
#   1. an MM_* ENVIRONMENT VARIABLE    wins, so a PEST worker or a batch run
#                                      can override anything for one run
#   2. code/configs/paths.local.toml   what PANEL 0 writes. Machine-local and
#                                      NOT tracked: it records where one
#                                      person's data sits, which is not a
#                                      property of the project
#   3. the defaults below
#
# Panel 0 writes (2) instead of editing this file, which is tracked: a
# generated edit here would show up as a diff on every machine, and a Windows
# path written into Python source is one backslash from a broken escape.
SETTINGS = REPO / 'code' / 'configs' / 'paths.local.toml'

# key -> (environment variable, default, what it is). This IS panel 0's form.
SETTABLE = {
    'example_root': ('MM_EXAMPLE_ROOT', '',
                     'Holds one folder per case: the Tier-A inputs a run '
                     'reads, and the only place the model reads input from. '
                     'Blank = <repo>/example, which is where it belongs.'),
    'data_root':   ('MM_DATA_ROOT', r'E:\00code_ws\LAMATA_new',
                    'Raw geospatial sources, OUTSIDE the repository.'),
    'gis':         ('MM_GIS_WS', '',
                    'The shapefiles, read ONLY by the converter -- never by a '
                    'run. Blank = DATA_ROOT/GIS.'),
    'ws_root':     ('MM_WS_ROOT', r'E:\00code_ws\LaMata_MM-MF6',
                    'Root of ALL run output. Never inside the repository.'),
    'nwt_ref':     ('MM_NWT_REF',
                    r'E:\00code_ws\LaMata_new_PhD_artigo_2s3L\_h5_MM.h5',
                    'The legacy MODFLOW-NWT _h5_MM.h5 the comparison figures '
                    'read. A FILE, not a folder.'),
    'modflow_dir': ('MM_MODFLOW_DIR', r'C:\00MODFLOW',
                    'Holds the MODFLOW 6, Triangle, GRIDGEN and PEST++ '
                    'executables.'),
    'python_exe':  ('MM_PYTHON_EXE', r'C:\miniconda3\envs\flopy\python.exe',
                    'The conda-env python a PEST++ worker runs forward_run '
                    'with. A FILE, not a folder.'),
}


def read_settings():
    """What ``paths.local.toml`` holds, or ``{}``. Never raises: a missing or
    malformed file must fall back to the defaults, not stop every import."""
    try:
        import tomllib
    except ImportError:                                   # pragma: no cover
        try:
            import tomli as tomllib
        except ImportError:
            return {}
    try:
        with open(SETTINGS, 'rb') as fh:
            return {k: str(v) for k, v in
                    (tomllib.load(fh).get('paths') or {}).items()}
    except (OSError, ValueError):
        return {}


def _setting(key):
    """One key, resolved: environment, then the local file, then the default."""
    env, default, _doc = SETTABLE[key]
    return (os.environ.get(env) or _SETTINGS.get(key) or default or '')


def source_of(key):
    """WHICH of the three supplied this key, so panel 0 can say so."""
    env, _default, _doc = SETTABLE[key]
    if os.environ.get(env):
        return 'the %s environment variable' % env
    if _SETTINGS.get(key):
        return 'code/configs/paths.local.toml'
    return 'the built-in default'


def settings_label():
    """The settings file, repo-relative when it is inside the repo.

    ``relative_to`` RAISES when it is not, which is exactly the case a test
    or a relocated checkout produces -- so it is never called unguarded.
    """
    try:
        return str(SETTINGS.relative_to(REPO))
    except ValueError:
        return str(SETTINGS)


def save_settings(values):
    """Write the machine-local paths, and apply them in THIS process.

    Only non-blank values are written: a key left empty goes back to its
    default rather than being pinned to an empty string.
    """
    SETTINGS.parent.mkdir(parents=True, exist_ok=True)
    out = ["# Machine-local paths, written by panel 0. NOT tracked: this says",
           "# where one person's data sits, which is not a property of the",
           '# project. An MM_* environment variable still wins over it.',
           '[paths]']
    for key in SETTABLE:
        v = str(values.get(key, '') or '').strip()
        if v:
            out.append('%-12s = "%s"' % (key, v.replace('\\', '\\\\')))
    SETTINGS.write_text('\n'.join(out) + '\n', encoding='utf-8')
    return reload_paths()


def reload_paths():
    """Re-resolve every root after the settings changed, without a restart."""
    global _SETTINGS, DATA_ROOT, WS_ROOT, NWT_REF, MODFLOW_DIR, PYTHON_EXE
    global GIS, EXAMPLE_ROOT, LIBMF6, MF6_EXE, TRIANGLE_EXE, PESTPP_IES
    global GRIDGEN_EXE
    _SETTINGS = read_settings()
    DATA_ROOT = Path(_setting('data_root'))
    WS_ROOT = Path(_setting('ws_root'))
    NWT_REF = Path(_setting('nwt_ref'))
    MODFLOW_DIR = Path(_setting('modflow_dir'))
    PYTHON_EXE = _setting('python_exe')
    EXAMPLE_ROOT = Path(_setting('example_root') or str(REPO / 'example'))
    GIS = Path(_setting('gis') or str(DATA_ROOT / 'GIS'))
    LIBMF6 = str(MODFLOW_DIR / 'mf6.7.0_win64' / 'bin' / 'libmf6.dll')
    MF6_EXE = str(MODFLOW_DIR / 'mf6.7.0_win64' / 'bin' / 'mf6.exe')
    TRIANGLE_EXE = str(MODFLOW_DIR / 'win64' / 'triangle.exe')
    PESTPP_IES = str(MODFLOW_DIR / 'pestpp' / 'pestpp-ies.exe')
    # GRIDGEN ships two builds side by side; x64 is the one to use (WP1c.1).
    GRIDGEN_EXE = os.environ.get(
        'MM_GRIDGEN_EXE',
        str(MODFLOW_DIR / 'gridgen.1.0.02' / 'bin' / 'gridgen_x64.exe'))
    return _SETTINGS


_SETTINGS = read_settings()
DATA_ROOT   = Path(_setting('data_root'))
WS_ROOT     = Path(_setting('ws_root'))
NWT_REF     = Path(_setting('nwt_ref'))
MODFLOW_DIR = Path(_setting('modflow_dir'))
PYTHON_EXE  = _setting('python_exe')

# ==== derived (usually no need to edit) =====================================
EXAMPLE_ROOT = Path(_setting('example_root') or str(REPO / 'example'))
GIS          = Path(_setting('gis') or str(DATA_ROOT / 'GIS'))
LIBMF6       = str(MODFLOW_DIR / 'mf6.7.0_win64' / 'bin' / 'libmf6.dll')
MF6_EXE      = str(MODFLOW_DIR / 'mf6.7.0_win64' / 'bin' / 'mf6.exe')
TRIANGLE_EXE = str(MODFLOW_DIR / 'win64' / 'triangle.exe')
PESTPP_IES   = str(MODFLOW_DIR / 'pestpp' / 'pestpp-ies.exe')
# GRIDGEN ships two builds side by side; the x64 one is the one to use (WP1c.1).
GRIDGEN_EXE  = os.environ.get(
    'MM_GRIDGEN_EXE',
    str(MODFLOW_DIR / 'gridgen.1.0.02' / 'bin' / 'gridgen_x64.exe'))

# Backwards compatibility with the pre-WP0 environment variables, so an
# existing shell keeps working while the flags are being retired.
if 'MARMITES_WS_ROOT' in os.environ and 'MM_WS_ROOT' not in os.environ:
    WS_ROOT = Path(os.environ['MARMITES_WS_ROOT'])
if 'MARMITES_GIS_WS' in os.environ and 'MM_GIS_WS' not in os.environ:
    GIS = Path(os.environ['MARMITES_GIS_WS'])
if 'MARMITES_NWT_REF' in os.environ and 'MM_NWT_REF' not in os.environ:
    NWT_REF = Path(os.environ['MARMITES_NWT_REF'])


def dataset_dir(case='LaMata'):
    """Tier A inputs for a case study: ``<EXAMPLE_ROOT>/<case>/``.

    EXAMPLE_ROOT defaults to ``<repo>/example``, and that is where it belongs:
    this is the ONLY folder inside the repository the model reads from, and it
    holds strictly the files MM and MF read directly. Panel 0 can point it
    elsewhere for a case whose inputs are kept outside the checkout.
    """
    return EXAMPLE_ROOT / str(case)


# Files renamed in WP1.2. `inputPONDw.asc` and `inputPONDhmax.asc` were a legacy
# MISNOMER: they hold the STREAM network (channel width and depth), which is why
# the arrays they feed have always been called gridSsurfw / gridSsurfhmax. The
# shim keeps an un-migrated dataset working for one release.
LEGACY_ALIASES = {
    'inputSTREAMw.asc': 'inputPONDw.asc',
    'inputSTREAMhmax.asc': 'inputPONDhmax.asc',
}


def resolve_input(name, case='LaMata', warn=True):
    """Absolute path of a Tier-A input, accepting the pre-WP1.2 name.

    Prefers the current name; falls back to the legacy one with a note, so a
    dataset that has not been renamed still runs.
    """
    ds = dataset_dir(case)
    p = ds / name
    if p.exists():
        return str(p)
    legacy = LEGACY_ALIASES.get(os.path.basename(name))
    if legacy:
        alt = p.parent / legacy
        if alt.exists():
            if warn:
                print('note: %s not found; using the legacy name %s. That file '
                      'holds the STREAM network, not ponds -- rename it '
                      '(WP1.2).' % (name, legacy))
            return str(alt)
    return str(p)          # let the caller raise its own, more specific error


def report_paths(case='LaMata', stream=None):
    """Print a one-line summary of every root and flag the missing ones.

    Returns the number of roots that do not exist, so a caller can fail fast.
    """
    import sys
    out = stream or sys.stdout
    rows = [
        ('REPO',         REPO,               True),
        ('DATASET',      dataset_dir(case),  True),
        ('DATA_ROOT',    DATA_ROOT,          False),
        ('GIS',          GIS,                False),
        ('WS_ROOT',      WS_ROOT,            False),
        ('NWT_REF',      NWT_REF,            False),
        ('MODFLOW_DIR',  MODFLOW_DIR,        False),
        ('LIBMF6',       Path(LIBMF6),       False),
        ('MF6_EXE',      Path(MF6_EXE),      False),
        ('TRIANGLE_EXE', Path(TRIANGLE_EXE), False),
        ('GRIDGEN_EXE',  Path(GRIDGEN_EXE),  False),
        ('PESTPP_IES',   Path(PESTPP_IES),   False),
        ('PYTHON_EXE',   Path(PYTHON_EXE),   False),
    ]
    missing = 0
    out.write('MARMITES paths (case=%s)\n' % case)
    for name, p, required in rows:
        ok = Path(p).exists()
        if not ok:
            missing += 1
        out.write('  %-13s %-6s %s%s\n'
                  % (name, 'OK' if ok else 'MISSING', p,
                     '   <-- required' if (required and not ok) else ''))
    if missing:
        out.write('  %d path(s) missing -- set them on PANEL 0 of the '
                  'front-end, which writes %s.\n' % (missing, settings_label()))
    return missing


if __name__ == '__main__':
    import sys
    sys.exit(1 if report_paths(*(sys.argv[1:2] or ['LaMata'])) and False else 0)
