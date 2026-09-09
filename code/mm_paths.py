# -*- coding: utf-8 -*-
"""The only machine-specific paths (WP0, step 0.2).

Every other module imports its roots from here, so moving the project to a new
machine means editing the roots below -- or setting the matching ``MM_*``
environment variables -- and nothing else. This is the MARMITES equivalent of
``CdL/code/config.py`` in the MF6models repository.

    REPO         the MARMITES checkout (self-locating; never edit)
    DATA_ROOT    raw geospatial sources, OUTSIDE the repo (Tier B)
    WS_ROOT      root of ALL run output, OUTSIDE the repo
    NWT_REF      the legacy MODFLOW-NWT _h5_MM.h5 used by the comparison figures
    MODFLOW_DIR  folder holding the MODFLOW 6 / PEST++ executables
    PYTHON_EXE   the conda-env python the PEST++ workers run forward_run with

Everything else is derived. ``DATASET`` resolves to ``example/<case>/`` and is
the only place MM and MF read input from; nothing here writes into the repo.

Run ``python code/mm_paths.py`` to print what resolved on this machine and
which roots are missing.
"""

import os
from pathlib import Path

__all__ = ['REPO', 'DATA_ROOT', 'WS_ROOT', 'NWT_REF', 'MODFLOW_DIR', 'GIS',
           'LIBMF6', 'MF6_EXE', 'PESTPP_IES', 'PYTHON_EXE', 'TRIANGLE_EXE',
           'dataset_dir', 'report_paths']

# ==== the machine-specific roots (edit these, or set the MM_* env vars) ======
REPO        = Path(__file__).resolve().parents[1]
DATA_ROOT   = Path(os.environ.get('MM_DATA_ROOT',   r'E:\00code_ws\LAMATA_new'))
WS_ROOT     = Path(os.environ.get('MM_WS_ROOT',     r'E:\00code_ws\LaMata_MM-MF6'))
NWT_REF     = Path(os.environ.get('MM_NWT_REF',
                                  r'E:\00code_ws\LaMata_new_PhD_artigo_2s3L\_h5_MM.h5'))
MODFLOW_DIR = Path(os.environ.get('MM_MODFLOW_DIR', r'C:\00MODFLOW'))
PYTHON_EXE  = os.environ.get('MM_PYTHON_EXE', r'C:\miniconda3\envs\flopy\python.exe')

# ==== derived (usually no need to edit) =====================================
GIS          = Path(os.environ.get('MM_GIS_WS', str(DATA_ROOT / 'GIS')))
LIBMF6       = str(MODFLOW_DIR / 'mf6.7.0_win64' / 'bin' / 'libmf6.dll')
MF6_EXE      = str(MODFLOW_DIR / 'mf6.7.0_win64' / 'bin' / 'mf6.exe')
TRIANGLE_EXE = str(MODFLOW_DIR / 'win64' / 'triangle.exe')
PESTPP_IES   = str(MODFLOW_DIR / 'pestpp' / 'pestpp-ies.exe')

# Backwards compatibility with the pre-WP0 environment variables, so an
# existing shell keeps working while the flags are being retired.
if 'MARMITES_WS_ROOT' in os.environ and 'MM_WS_ROOT' not in os.environ:
    WS_ROOT = Path(os.environ['MARMITES_WS_ROOT'])
if 'MARMITES_GIS_WS' in os.environ and 'MM_GIS_WS' not in os.environ:
    GIS = Path(os.environ['MARMITES_GIS_WS'])
if 'MARMITES_NWT_REF' in os.environ and 'MM_NWT_REF' not in os.environ:
    NWT_REF = Path(os.environ['MARMITES_NWT_REF'])


def dataset_dir(case='LaMata'):
    """Tier A inputs for a case study: ``example/<case>/``.

    This is the ONLY folder inside the repository the model reads from, and it
    holds strictly the files MM and MF read directly.
    """
    return REPO / 'example' / str(case)


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
        out.write('  %d path(s) missing -- edit code/mm_paths.py or set the '
                  'MM_* environment variables.\n' % missing)
    return missing


if __name__ == '__main__':
    import sys
    sys.exit(1 if report_paths(*(sys.argv[1:2] or ['LaMata'])) and False else 0)
