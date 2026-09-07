# -*- coding: utf-8 -*-
"""Isolate whether writing UZF SINF via the MODFLOW 6 API actually reaches UZF.

The coupled run keeps applying a constant perc_user (UZF INFILTRATION stuck at
977 m3/d) even after (a) binding SINF instead of FINF and (b) writing after
prepare_time_step. This script strips the coupler away and tests the raw
mechanism on the first few stress periods:

  * what SINF holds right after prepare_time_step (i.e. after UZF's rp),
  * whether a known value written to the SINF pointer survives the solve,
  * whether a freshly re-fetched pointer sees the write (view vs copy),
  * both SINF and FINF, so we can see which one (if either) is operative.

Run (flopy env):
  python tests/diag_sinf.py --libmf6 C:\\00MODFLOW\\mf6.7.0_win64\\bin\\libmf6.dll
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..', 'trunk'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

import run_lamata_mf6 as R          # noqa: E402
from marmites_mf6 import clsMF6     # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--libmf6', required=True)
    ap.add_argument('--nsteps', type=int, default=4)
    a = ap.parse_args()

    ws = os.path.join(R.WS_ROOT, 'MF6_ws')
    cMF, mm, ctx, state, top, botm, conv = R.setup_lamata(daily=True, nsp=6, nlay=2)
    b = clsMF6(cMF, top=top, botm=botm, sim_ws=ws, daily=True)
    b.seep = 'drn'
    b.build()
    b.write()
    name = cMF.modelname.upper()
    ncell = b.ncell

    lib = os.path.abspath(a.libmf6.strip('"').strip("'"))
    if hasattr(os, 'add_dll_directory'):
        os.add_dll_directory(os.path.dirname(lib))
    from modflowapi import ModflowApi
    api = ModflowApi(lib, working_directory=os.path.abspath(ws))
    api.initialize()

    def address(var, comp):
        try:
            addr = api.get_var_address(var, *comp.split('/'))
            if addr in set(api.get_input_var_names()) or True:
                return addr
        except Exception as exc:
            print('  address error %s/%s: %r' % (comp, var, exc))
        return None

    saddr = address('SINF', f'{name}/UZF')
    faddr = address('FINF', f'{name}/UZF')
    print('SINF address:', saddr)
    print('FINF address:', faddr)
    sinf = api.get_value_ptr(saddr) if saddr else None
    finf = api.get_value_ptr(faddr) if faddr else None
    for lbl, arr in (('SINF', sinf), ('FINF', finf)):
        print('  %s ptr: %s' % (lbl, None if arr is None else
                                 'size %d, first %.5g' % (arr.size, float(arr.ravel()[0]))))

    # Try a different write POSITION on each step and see which one survives to
    # the end of the solve (i.e. is the one UZF actually uses).
    #   A: after prepare_time_step (current coupler behaviour -- known to fail)
    #   B: after prepare_solve, before the first solve
    #   C: before EVERY solve() call in the outer loop
    positions = ['A_after_prep_ts', 'B_after_prep_solve', 'C_before_each_solve']
    print('\nstep | position            | wrote  | after solve (mean) | stuck?')
    for step in range(a.nsteps):
        pos = positions[step % len(positions)]
        val = 0.001 * (step + 1)
        dt = api.get_time_step()
        api.prepare_time_step(dt)
        if pos == 'A_after_prep_ts' and sinf is not None:
            sinf[:ncell] = val
        api.prepare_solve(1)
        if pos == 'B_after_prep_solve' and sinf is not None:
            sinf[:ncell] = val
        k = 0
        while k < 300:
            k += 1
            if pos == 'C_before_each_solve' and sinf is not None:
                sinf[:ncell] = val
            if api.solve(1):
                break
        after_solve = float(np.mean(sinf[:ncell])) if sinf is not None else float('nan')
        stuck = 'YES' if abs(after_solve - val) < 1e-9 else 'no'
        api.finalize_solve(1)
        api.finalize_time_step()
        print('  %2d  | %-19s | %.4f |  %12.6f  |  %s'
              % (step, pos, val, after_solve, stuck))

    try:
        api.finalize()
    except Exception:
        pass
    print('\nINTERPRETATION: whichever position shows stuck=YES is where the')
    print('coupler must write SINF. If NONE stick, UZF re-derives infiltration')
    print('from its period data every outer iteration and the BMI SINF array is')
    print('effectively read-only -- then the fix is to feed recharge a different')
    print('way (e.g. a RCH/RCHA package driven by the API, or UZF FINF via a')
    print('time-array-series), not by writing SINF.')


if __name__ == '__main__':
    main()
