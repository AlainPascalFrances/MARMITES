# -*- coding: utf-8 -*-
"""Start the front-end -- from a local mirror when the checkout cannot host it.

Streamlit does not run from a mapped network drive. The checkout stays the
one place the project lives -- the code, the run configurations, the machine
settings (code/configs/paths.local.toml), the dataset -- and the app runs
from a MIRROR of the code on a local disk, refreshed at every launch:

    python code/tools/launch_app.py                    # mirror when needed
    python code/tools/launch_app.py --mirror always    # mirror anyway
    python code/tools/launch_app.py --mirror-dir D:/mm_app
    python code/tools/launch_app.py -- --server.port 8502   # to streamlit

The mirror receives ``code/`` -- without ``code/configs`` and the caches --
and ``.streamlit/``, nothing else. A marker in it (``mm_paths.MIRROR_MARK``)
names the checkout, so the app, the runs it launches and the converter it
starts read and write the configurations, the settings and the dataset IN
THE CHECKOUT: a configuration saved on a panel is saved in the repository,
and there is no second copy to fall behind. A file deleted from the
checkout is deleted from the mirror at the next launch.

``--mirror auto`` (the default) mirrors when the checkout is on a network
drive (a mapped letter or a UNC path; Windows) and runs in place otherwise.
The mirror defaults to ``<LOCALAPPDATA>/MARMITES/mirror/<checkout>-<id>``
(``~/.cache/MARMITES/...`` without LOCALAPPDATA). It must be a folder of its
own -- not a git checkout, not inside this one -- because the refresh
deletes what the checkout does not have. Edit the code in the checkout,
never in the mirror: the next launch overwrites it.
"""

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

CODE = Path(os.path.abspath(__file__)).parents[1]
REPO = CODE.parent
if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))
import mm_paths  # noqa: E402

__all__ = ['TREES', 'on_network_drive', 'default_mirror', 'check_mirror_dir',
           'sync', 'mark', 'main']

# what is mirrored: (folder of the checkout, its sub-folders left out)
TREES = (('code', ('configs',)), ('.streamlit', ()))
SKIP_DIRS = {'__pycache__', '.pytest_cache', '.ruff_cache',
             '.ipynb_checkpoints', '.git'}
SKIP_SUFFIXES = ('.pyc', '.pyo')


def on_network_drive(path):
    """True for a UNC path, or a drive letter Windows maps to a share."""
    p = str(path)
    if p.startswith('\\\\') or p.startswith('//'):
        return True
    if os.name != 'nt':
        return False
    drive = os.path.splitdrive(os.path.abspath(p))[0]
    if not drive:
        return False
    import ctypes
    DRIVE_REMOTE = 4
    return ctypes.windll.kernel32.GetDriveTypeW(drive + '\\') == DRIVE_REMOTE


def default_mirror(repo):
    """A local folder of its own for this checkout's mirror."""
    base = os.environ.get('LOCALAPPDATA') or os.path.join(
        os.path.expanduser('~'), '.cache')
    tag = hashlib.sha1(str(repo).lower().encode('utf-8')).hexdigest()[:8]
    return Path(base) / 'MARMITES' / 'mirror' / ('%s-%s' % (Path(repo).name,
                                                            tag))


def check_mirror_dir(mirror, repo):
    """Refuse a folder the refresh could damage: the checkout itself, a
    folder inside it or around it, or a git checkout of its own."""
    m, r = Path(mirror).resolve(), Path(repo).resolve()
    if m == r or r in m.parents or m in r.parents:
        raise SystemExit('the mirror %s overlaps the checkout %s: give a '
                         'folder of its own (--mirror-dir)' % (mirror, repo))
    if (m / '.git').exists():
        raise SystemExit('%s is a git checkout. The mirror is refreshed by '
                         'deleting what the checkout does not have, so it '
                         'must be a folder of its own (--mirror-dir).' % mirror)


def _same(src, dst):
    try:
        a, b = src.stat(), dst.stat()
    except OSError:
        return False
    return a.st_size == b.st_size and abs(a.st_mtime - b.st_mtime) < 1.0


def sync(repo, mirror):
    """Make the mirror's code the checkout's. Returns (copied, removed).

    Copies what is new or changed (size or time), deletes what the checkout
    no longer has, and never writes to the checkout.
    """
    repo, mirror = Path(repo), Path(mirror)
    copied = removed = 0
    for top, leave_out in TREES:
        src_top, dst_top = repo / top, mirror / top
        wanted = set()
        if src_top.is_dir():
            for root, dirs, files in os.walk(src_top):
                rel = Path(root).relative_to(src_top)
                dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS
                                 and not (rel == Path('.') and d in leave_out))
                for f in files:
                    if f.endswith(SKIP_SUFFIXES) or f == mm_paths.MIRROR_MARK:
                        continue
                    r = rel / f
                    wanted.add(r)
                    if not _same(src_top / r, dst_top / r):
                        (dst_top / r).parent.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(src_top / r, dst_top / r)
                        copied += 1
        if not dst_top.is_dir():
            continue
        for root, dirs, files in os.walk(dst_top, topdown=False):
            rel = Path(root).relative_to(dst_top)
            if any(part in SKIP_DIRS for part in rel.parts):
                continue                     # caches stay, Python rebuilds them
            for f in files:
                r = rel / f
                if (r in wanted or f.endswith(SKIP_SUFFIXES)
                        or (rel == Path('.') and f == mm_paths.MIRROR_MARK)):
                    continue
                (dst_top / r).unlink()
                removed += 1
            if rel != Path('.') and not any(Path(root).iterdir()):
                Path(root).rmdir()
    return copied, removed


def mark(repo, mirror):
    """Name the checkout in the mirror, for mm_paths."""
    fn = Path(mirror) / 'code' / mm_paths.MIRROR_MARK
    fn.parent.mkdir(parents=True, exist_ok=True)
    fn.write_text(str(repo), encoding='utf-8')
    return fn


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Start the MARMITES front-end, from a local mirror of '
                    'the code when the checkout is on a network drive.')
    ap.add_argument('--mirror', choices=('auto', 'always', 'never'),
                    default='auto')
    ap.add_argument('--mirror-dir', default=None,
                    help='the local folder to mirror into')
    ap.add_argument('--sync-only', action='store_true',
                    help='refresh the mirror and stop')
    ap.add_argument('streamlit_args', nargs=argparse.REMAINDER,
                    help='after --: passed to streamlit run')
    a = ap.parse_args(argv)
    extra = list(a.streamlit_args)
    if extra and extra[0] == '--':
        extra = extra[1:]

    root = REPO
    if a.mirror == 'always' or (a.mirror == 'auto'
                                and on_network_drive(REPO)):
        mirror = Path(a.mirror_dir) if a.mirror_dir else default_mirror(REPO)
        check_mirror_dir(mirror, REPO)
        copied, removed = sync(REPO, mirror)
        mark(REPO, mirror)
        print('mirror: %s\n   refreshed from %s: %d file(s) copied, %d removed'
              % (mirror, REPO, copied, removed))
        print('configurations and settings: %s' % (REPO / 'code' / 'configs'))
        root = mirror
    elif a.sync_only:
        print('the checkout %s runs in place: nothing to mirror' % REPO)
    if a.sync_only:
        return 0

    cmd = [sys.executable, '-m', 'streamlit', 'run',
           str(root / 'code' / 'app' / 'Home.py')] + extra
    print(' '.join(cmd))
    try:
        return subprocess.call(cmd, cwd=str(root))   # .streamlit/ is found here
    except KeyboardInterrupt:
        return 130


if __name__ == '__main__':
    sys.exit(main())
