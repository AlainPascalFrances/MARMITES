# -*- coding: utf-8 -*-
"""Bisect the MARMITES <-> MF6 API coupling failure.

The standalone `mf6.exe` run of the same model terminates normally, so the
groundwater model is sound and the fault lies in what the Python side does
through the API. This script isolates it by staging the exchange:

  --level 0   initialize, bind pointers, run the time loop, write NOTHING
              (pure MF6 driven through the API -- equivalent to mf6.exe)
  --level 1   + write UZF FINF only
  --level 2   + write WEL rates only (no FINF)
  --level 3   + write both, but never read heads back
  --level 4   + read heads/exfiltration back (full exchange, no MARMITES)

The first level that crashes names the culprit. Each level prints what it
did before every MF6 call, so a silent process death still tells us the
exact last operation.

Usage:
    python tests/diagnose_coupling.py --libmf6 C:\\...\\libmf6.dll --level 0 --nsp 3
    (repeat with --level 1, 2, 3, 4)
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TRUNK = os.path.abspath(os.path.join(HERE, '..'))
DS = os.path.abspath(os.path.join(HERE, '..', '..', 'example', 'LaMata'))
WS_ROOT = os.environ.get('MARMITES_WS_ROOT', os.path.join('E:' + os.sep, '00code_ws', 'LaMata_MM-MF6'))
for p in ('', 'MARMITESutilities', 'MARMITESsoil', 'ppMF_FloPy', 'ppMF6'):
    sys.path.insert(0, os.path.join(TRUNK, p))

import matplotlib  # noqa: E402
matplotlib.use('agg')


def say(msg):
    print('  [diag] %s' % msg, flush=True)     # flush: survives a hard crash


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--libmf6', required=True)
    ap.add_argument('--level', type=int, default=0, choices=[0, 1, 2, 3, 4, 5])
    ap.add_argument('--nsp', type=int, default=3)
    ap.add_argument('--ws', default=None)
    ap.add_argument('--grid', choices=['dis', 'disv'], default='dis')
    a = ap.parse_args()
    ws = os.path.abspath(a.ws or os.path.join(WS_ROOT, 'MF6_ws'))

    # --- build the model exactly as the coupled runner does ---------------
    from run_lamata_mf6 import setup_lamata
    from marmites_mf6 import clsMF6
    cMF, mm, ctx, state, top, botm, conv_fact = setup_lamata(daily=True, nsp=a.nsp,
                                                             grid=a.grid)
    b = clsMF6(cMF, top=top, botm=botm, sim_ws=ws, daily=True, grid=a.grid)
    b.build()
    b.write()
    say('model written to %s (%d SPs incl. steady)' % (ws, b.nper))

    lib = os.path.abspath(a.libmf6.strip('"').strip("'"))
    if hasattr(os, 'add_dll_directory'):
        os.add_dll_directory(os.path.dirname(lib))
    from modflowapi import ModflowApi
    api = ModflowApi(lib, working_directory=ws)

    say('initialize()')
    api.initialize()
    name = cMF.modelname.upper()

    # --- bind (read-only inspection first) --------------------------------
    say('listing variables')
    names = set(api.get_input_var_names())

    def ptr(var, comp):
        addr = api.get_var_address(var, *comp.split('/'))
        if addr not in names:
            say('  %s NOT exposed' % addr)
            return None, addr
        p = api.get_value_ptr(addr)
        say('  %-28s shape=%s dtype=%s' % (addr, getattr(p, 'shape', None),
                                           getattr(p, 'dtype', None)))
        return p, addr

    say('binding pointers:')
    p_x, _ = ptr('X', name)
    p_finf, _ = ptr('FINF', f'{name}/UZF')
    p_gwd, _ = ptr('GWD', f'{name}/UZF')
    p_q, _ = ptr('Q', f'{name}/WEL')
    p_bound, _ = ptr('BOUND', f'{name}/WEL')
    ncell = ctx.ncell
    say('MARMITES cells = %d ; NODES(X) = %s ; UZF cells = %s'
        % (ncell, getattr(p_x, 'size', '?'), getattr(p_finf, 'size', '?')))

    # node mapping (needed only at level 4)
    x_index = None
    if a.level >= 4:
        i_arr = np.array([c[1] for c in ctx.cells], dtype=int)
        j_arr = np.array([c[2] for c in ctx.cells], dtype=int)
        k_arr = np.array([b.surf_cells[n][2] for n in range(ncell)], dtype=int)
        user = k_arr * (cMF.nrow * cMF.ncol) + i_arr * cMF.ncol + j_arr
        try:
            nu = np.asarray(api.get_value(api.get_var_address('NODEUSER', name, 'DIS')))
            say('NODEUSER size=%d (X size=%d)' % (nu.size, p_x.size))
            lut = {int(u) - 1: r for r, u in enumerate(nu)}
            x_index = np.array([lut[int(u)] for u in user], dtype=int)
            say('node mapping via NODEUSER ok, max index=%d' % x_index.max())
        except Exception as exc:
            say('NODEUSER unavailable (%r) -> identity mapping' % (exc,))
            x_index = user
            say('identity max index=%d vs X size=%d' % (x_index.max(), p_x.size))

    # --- level 5: does writing Q actually change the solution? -------------
    # A silent no-op would be worse than a crash, so compare heads after the
    # same stress period with no pumping and with strong pumping.
    if a.level == 5:
        p_sim, _ = ptr('SIMVALS', f'{name}/WEL')
        i_arr = np.array([c[1] for c in ctx.cells], dtype=int)
        j_arr = np.array([c[2] for c in ctx.cells], dtype=int)
        k_arr = np.array([b.surf_cells[n][2] for n in range(ncell)], dtype=int)
        user = k_arr * (cMF.nrow * cMF.ncol) + i_arr * cMF.ncol + j_arr
        try:
            nu = np.asarray(api.get_value(api.get_var_address('NODEUSER', name, 'DIS')))
            lut = {int(u) - 1: r for r, u in enumerate(nu)}
            xi = np.array([lut[int(u)] for u in user], dtype=int)
        except Exception:
            xi = user
        heads = {}
        for label, qval in (('no pumping', 0.0), ('pumping -50 m3/d', -50.0)):
            say('--- %s ---' % label)
            if p_finf is not None:
                p_finf[:ncell] = 2.0e-4
            if p_q is not None:
                p_q[:ncell] = qval
            dt = api.get_time_step()
            api.prepare_time_step(dt)
            api.prepare_solve(1)
            k = 0
            while k < 200:
                k += 1
                if api.solve(1):
                    break
            api.finalize_solve(1)
            api.finalize_time_step()
            h = np.asarray(p_x)[xi].copy()
            heads[label] = h
            say('%s: solved in %d iters, head mean=%.4f' % (label, k, h.mean()))
            if p_sim is not None:
                s = np.asarray(p_sim).ravel()[:ncell]
                say('%s: WEL SIMVALS sum=%.4g (expected ~ %.4g)'
                    % (label, s.sum(), qval * ncell))
        d = heads['pumping -50 m3/d'] - heads['no pumping']
        say('head change from pumping: mean=%.4g m, min=%.4g m' % (d.mean(), d.min()))
        if np.allclose(d, 0.0):
            say('*** Q WRITE HAD NO EFFECT -- the rate is not reaching MF6 ***')
        else:
            say('*** Q WRITE IS EFFECTIVE (heads responded to pumping) ***')
        try:
            api.finalize()
            say('finalize() ok')
        except Exception as exc:
            say('finalize() failed: %r' % (exc,))
        return

    # --- time loop ---------------------------------------------------------
    nper_total = b.nper
    say('starting time loop over %d MF6 stress periods, level=%d' % (nper_total, a.level))
    try:
        for n in range(nper_total):
            tag = 'SP%d(%s)' % (n, 'steady' if n == 0 else 'transient')
            if a.level >= 1 and p_finf is not None:
                say('%s writing FINF[:%d] = 2e-4' % (tag, ncell))
                p_finf[:ncell] = 2.0e-4
            if a.level >= 2:
                qv = np.full(ncell, -1.0e-3)
                if p_q is not None:
                    say('%s writing Q[:%d]' % (tag, ncell))
                    p_q[:ncell] = qv
                if p_bound is not None and a.level >= 3:
                    say('%s writing BOUND (shape %s)' % (tag, p_bound.shape))
                    if p_bound.ndim == 1:
                        p_bound[:ncell] = qv
                    elif p_bound.shape[0] >= ncell:
                        p_bound[:ncell, 0] = qv
                    else:
                        p_bound[0, :ncell] = qv
            dt = api.get_time_step()
            say('%s get_time_step -> %r' % (tag, dt))
            say('%s prepare_time_step' % tag)
            api.prepare_time_step(dt)
            say('%s prepare_solve' % tag)
            api.prepare_solve(1)
            kiter = 0
            while kiter < 200:
                kiter += 1
                if api.solve(1):
                    break
            say('%s solved in %d outer iterations' % (tag, kiter))
            api.finalize_solve(1)
            api.finalize_time_step()
            if a.level >= 4 and x_index is not None:
                h = np.asarray(p_x)[x_index]
                say('%s heads min/mean/max = %.3f / %.3f / %.3f'
                    % (tag, h.min(), h.mean(), h.max()))
                if p_gwd is not None:
                    g = np.asarray(p_gwd).ravel()[:ncell]
                    say('%s GWD min/max = %.4g / %.4g' % (tag, g.min(), g.max()))
        say('TIME LOOP COMPLETED WITHOUT CRASH (level %d)' % a.level)
    finally:
        try:
            api.finalize()
            say('finalize() ok')
        except Exception as exc:
            say('finalize() failed: %r' % (exc,))


if __name__ == '__main__':
    main()
