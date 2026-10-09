# -*- coding: utf-8 -*-
"""Why La Mata does not converge at outer_dvclose 0.001: SFR evaporation
flip-flops on nearly dry reaches (2026-10-09).

Run 20261009165953 (outer_dvclose 0.001) failed 266 of 703 time steps, and
always at a layer-1 stream cell whose reach was nearly dry with the water
table just below the streambed top. MF6 6.7.0 (and develop, 2026-10-09):

  * sfr_solve takes the reach evaporation from the depth STORED by the
    previous solve: qe = EVAP * calc_surface_area_wet(n, this%depth(n)),
    a wetted area that ramps 0 -> w*L over the first 1e-5 m of depth;
  * when the reach's inflow is below E0*w*L and the head is below the
    streambed top, a stored depth > 1e-5 m makes the reach evaporate ALL
    its inflow (new depth 0), and a stored depth 0 makes it evaporate
    nothing and leak ALL its inflow (new depth ~1e-5 m);
  * sfr_fc's Picard loop (100 passes) never settles and ends in the state
    it started from; sfr_fn then perturbs the head by DEM4 = 1e-4 m, solves
    once, lands in the OTHER state, and hands GWF a derivative of
    +-inflow/1e-4 m (500 m2/d for 0.05 m3/d) where the true one is ~0;
  * so the Newton step at that cell creeps by about a millimetre per outer
    iteration, never under 0.001 m, always under 0.025 m.

La Mata's own output shows both states: on SP 43 reaches 4, 11 and 275 take
0.02-0.09 m3/d against 0.09-0.11 m3/d of potential evaporation and
evaporate all of it; on SP 42 reach 84 evaporates nothing and leaks all of
its 0.011 m3/d.

This script builds a valley with La Mata's stream-cell setup (Sy 0.01,
streambed 0.5 m below land, rhk 0.1, rbth 0.2, w 1.5 m, headwater reaches
fed 0.05 m3/d against E0 0.004 m/d x 30 m2), drives it through the API like
the coupler (ATS, EVAP written after prepare_time_step), and compares:

  e0   EVAP = E0 on every reach (what the coupler writes) at 0.001 / 0.025;
  cap  EVAP = min(E0, f * qin / (w L)), f = 0.5, at 0.001, where qin is the
       step's INFLOW (the runoff, known when the coupler writes) plus the
       upstream flow the last solve left in USFLOW (La Mata would add the
       pond outflow, QFROMMVR at the same memory path).

On 2026-10-09 (MF6 6.7.0): e0 at 0.001 failed 3 attempts (353 outer
iterations), e0 at 0.025 none (41), cap at 0.001 none (43), with 0.5 %
more stream evaporation than e0 at 0.025. Over La Mata's 43 saved days the
cap would remove 0.3 % of the stream evaporation (61 m3/d). The capped
reach can no longer evaporate more than f of what reaches it, so it never
reaches the dry/wet switch where the two states alternate.
A toy MF6 model, not La Mata.

    python code/tests/diag_sfr_evap_flipflop.py [--ws <folder>]
"""
import argparse
import os
import re
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..')))

NROW, NCOL, CELL = 7, 16, 20.0
SR = 3                               # the stream's row
NPER = 20
E0 = 0.004                           # m/d, open-water evaporation
RUNOFF = 0.05                        # m3/d into every 4th reach


def build(ws, exe):
    """The valley; outer_dvclose is set per run by set_dvclose()."""
    import flopy
    name = 'toy'
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=NPER, perioddata=[(1.0, 1, 1.0)] * NPER,
                          time_units='DAYS')
    ats = [(p, 1.0, 2e-4, 1.0, 2.0, 5.0) for p in range(NPER)]     # La Mata
    sim.get_package('tdis').ats.initialize(maxats=NPER, perioddata=ats,
                                           filename=name + '.tdis.ats')
    flopy.mf6.ModflowIms(sim, print_option='SUMMARY', complexity='COMPLEX',
                         outer_dvclose=0.001, outer_maximum=100,
                         under_relaxation='DBD', inner_dvclose=1e-4,
                         rcloserecord=[0.01, 'STRICT'],
                         linear_acceleration='BICGSTAB')
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True,
                               newtonoptions='UNDER_RELAXATION')
    x = np.arange(NCOL) * CELL
    top = np.array([100.0 - 0.01 * x + abs(i - SR) for i in range(NROW)])
    flopy.mf6.ModflowGwfdis(gwf, nlay=2, nrow=NROW, ncol=NCOL, delr=CELL,
                            delc=CELL, top=top,
                            botm=np.array([top - 20.0, top - 50.0]))
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=0.05, k33=0.05)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=[1e-4, 4e-7], sy=0.01,
                            transient={0: True})
    rtp = top[SR] - 0.5
    strt = np.array([top - 0.7, top - 0.7])
    strt[:, SR, :] = rtp - 0.002
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.0003)
    flopy.mf6.ModflowGwfdrn(gwf, auxiliary=['ddrn'], auxdepthname='ddrn',
                            pname='drn_seep', stress_period_data={0: [
                                [(0, i, j), top[i, j], 100.0, 0.125]
                                for i in range(NROW) for j in range(NCOL)
                                if i != SR]})
    slope = [(0.05, 1e-4, 0.01)[j % 3] for j in range(NCOL)]
    rd = [[j, (0, SR, j), CELL, 1.5, slope[j], rtp[j], 0.2, 0.1, 0.035,
           1 if j in (0, NCOL - 1) else 2, 1.0, 0] for j in range(NCOL)]
    cd = [[0, -1]] + [[j, j - 1, -(j + 1)] for j in range(1, NCOL - 1)] \
        + [[NCOL - 1, NCOL - 2]]
    spd = [[j, 'EVAPORATION', E0] for j in range(NCOL)] \
        + [[j, 'INFLOW', RUNOFF] for j in range(0, NCOL, 4)]
    flopy.mf6.ModflowGwfsfr(gwf, pname='sfr', nreaches=NCOL, packagedata=rd,
                            connectiondata=cd, unit_conversion=86400.0,
                            perioddata={0: spd}, save_flows=True)
    flopy.mf6.ModflowGwfoc(gwf, head_filerecord=name + '.hds',
                           saverecord=[('HEAD', 'ALL')])
    sim.write_simulation(silent=True)
    # La Mata's IMS has a plain inner_rclose; flopy wants a second word, and
    # STRICT would accept an outer iteration only when the linear solve
    # converges at its first inner iteration (ImsLinearBase testcnvg)
    p = os.path.join(ws, name + '.ims')
    txt = open(p).read()
    open(p, 'w').write(re.sub(r'(?i)(inner_rclose\s+\S+)\s+strict', r'\1', txt))


def set_dvclose(ws, dv):
    p = os.path.join(ws, 'toy.ims')
    txt = open(p).read()
    open(p, 'w').write(re.sub(r'(?i)(outer_dvclose\s+)\S+', r'\g<1>%g' % dv, txt))


def run(ws, rule, f, libmf6):
    """Drive the toy like the coupler: EVAP written after prepare_time_step,
    do_time_step (ATS retries inside). Returns (failed attempts, outer
    iterations, stream evaporation m3)."""
    from modflowapi import ModflowApi
    api = ModflowApi(libmf6, working_directory=ws)
    api.initialize(os.path.join(ws, 'mfsim.nam'))
    try:
        pre = 'TOY/SFR'
        evap = api.get_value_ptr(pre + '/EVAP')
        wl = api.get_value_ptr(pre + '/LENGTH') * 1.5
        usflow = api.get_value_ptr(pre + '/USFLOW')
        inflow = api.get_value_ptr(pre + '/INFLOW')
        simevap = api.get_value_ptr(pre + '/SIMEVAP')
        taken = 0.0
        while api.get_current_time() < api.get_end_time() - 1e-9:
            api.prepare_time_step(api.get_time_step())
            if rule == 'cap':
                # this step's INFLOW (the runoff, known) + the upstream flow
                # the last solve left in USFLOW (0 before the first step)
                qin = np.maximum(inflow + usflow, 0.0)
                evap[:] = np.minimum(E0, f * qin / wl)
            else:
                evap[:] = E0
            api.do_time_step()
            api.finalize_time_step()
            taken += simevap.sum() * api.get_time_step()
    finally:
        api.finalize()
    lst = open(os.path.join(ws, 'mfsim.lst'), errors='replace').read()
    return (len(re.findall(r'Solution 1 did not converge', lst)),
            len(re.findall(r'^\s*Model\s+\d+\s', lst, flags=re.M)), taken)


def main(argv=None):
    import mm_paths
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--ws', default=None,
                    help='a scratch folder (default: a temporary one)')
    a = ap.parse_args(argv)
    root = a.ws or tempfile.mkdtemp(prefix='mm_sfr_evap_flipflop_')
    for rule, dv, f in (('e0', 0.001, None), ('e0', 0.025, None),
                        ('cap', 0.001, 0.5)):
        ws = os.path.join(root, '%s_dv%g' % (rule, dv))
        if os.path.isdir(ws):
            shutil.rmtree(ws)
        os.makedirs(ws)
        build(ws, mm_paths.MF6_EXE)
        set_dvclose(ws, dv)
        nfail, nouter, taken = run(ws, rule, f, mm_paths.LIBMF6)
        print('EVAP %-14s outer_dvclose %-6g failed attempts %2d, '
              '%4d outer iterations, stream evaporation %.2f m3'
              % ('= E0' if rule == 'e0' else '<= %g x inflow' % f, dv,
                 nfail, nouter, taken))
    print('files in %s' % root)


if __name__ == '__main__':
    main()
