# -*- coding: utf-8 -*-
"""Run pestpp-ies on the La Mata template (WP7): the master and its local
PANTHER agents, in the PEST workspace (Home, pest_ws).

    python code/tools/pest_lamata_run.py --workers 6 --noptmax 0

``--noptmax 0`` is the base run through PEST++ itself (one forward run in
an agent: the check that the whole chain works), ``-1`` the prior Monte
Carlo, ``2``/``3`` the IES iterations. The template is copied into
``<pest_ws>/workers/worker_<k>`` (each worker its own dataset copy and
workspace) and the master runs in ``<pest_ws>/master``.

At most the server's physical cores minus two: the user keeps at least two
(CLAUDE.md rule 1, the PEST++ exception).
"""

import argparse
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
if CODE not in sys.path:
    sys.path.insert(0, CODE)

PHYSICAL_CORES = 8
LEAVE = 2


def main(argv=None):
    import pyemu
    import mm_paths
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--workers', type=int, default=PHYSICAL_CORES - LEAVE)
    ap.add_argument('--noptmax', type=int, required=True)
    ap.add_argument('--nreal', type=int, default=None,
                    help='ies_num_reals (default: as set up)')
    ap.add_argument('--port', type=int, default=4004)
    a = ap.parse_args(argv)
    if a.workers > PHYSICAL_CORES - LEAVE:
        raise SystemExit('at most %d workers: the user keeps %d cores'
                         % (PHYSICAL_CORES - LEAVE, LEAVE))
    pest_ws = str(mm_paths.PEST_WS or '')
    tpl = os.path.join(pest_ws, 'template')
    if not os.path.exists(os.path.join(tpl, 'lamata.pst')):
        raise SystemExit('no template in %s: run pest_lamata_setup.py' % tpl)
    pst = pyemu.Pst(os.path.join(tpl, 'lamata.pst'))
    pst.control_data.noptmax = a.noptmax
    if a.nreal:
        pst.pestpp_options['ies_num_reals'] = a.nreal
    pst.write(os.path.join(tpl, 'lamata.pst'), version=2)
    master = os.path.join(pest_ws, 'master')
    if os.path.isdir(master):
        shutil.rmtree(master)
    workers = os.path.join(pest_ws, 'workers')
    os.makedirs(workers, exist_ok=True)
    exe = str(mm_paths.PESTPP_IES)
    print('pestpp-ies %s, noptmax %d, %d worker(s), master %s'
          % (exe, a.noptmax, a.workers, master))
    pyemu.os_utils.start_workers(
        tpl, exe, 'lamata.pst', num_workers=a.workers,
        worker_root=workers, master_dir=master,
        port=a.port, cleanup=False)


if __name__ == '__main__':
    main()
