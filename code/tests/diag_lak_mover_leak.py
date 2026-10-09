# -*- coding: utf-8 -*-
"""The pond-outlet leak past the mover, on a toy chain of ponds (2026-10-09).

On La Mata's run 20261008114055 the LAK budget showed EXT-OUTFLOW at pond10
(-2,748 m3/yr) and pond13 (-535) although every outlet is moved to the
stream at FACTOR 1. EXT-OUTFLOW is the outlet's discharge MINUS what the
mover took from it (gwf-lak.f90: lak_get_external_outlet + _mover), and in
MF6 6.7:

  * the mover moves the PREVIOUS outer iteration's provider flows
    (gwf.f90: mvr_fc runs before the packages' fc; PackageMover qformvr is
    filled by lak_solve at the end of each fc);
  * those flows are ZEROED at the start of every time step
    (PackageMover%ad), so the day's water reaches a pond N links down a
    chain of ponds only after N outer iterations, from below;
  * LAK's convergence check takes the change of outlet discharge as a depth
    over the pond's area per step (lak_cc: dqout * delt / area) and the
    solver accepts it below outer_dvclose (NumericalSolution
    sln_package_convergence).

So a pond at the end of a chain can still be losing up to outer_dvclose x
area / delt when the step is accepted. On La Mata every pond stayed within
that bound; pond10 reached 96 % of it (37.2 of 38.7 m3/d) and pond13 92 %,
the two at the end of the longest chains (4 and 3 ponds upstream).

This script builds five ponds in series on a stream over an aquifer with
La Mata's K and Sy, runs it at outer_dvclose 0.025 / 0.001 / 1e-5 and
prints each pond's EXT-OUTFLOW. On 2026-10-09 (MF6 6.7.0) it gave, for the
last pond: -55 m3 (worst day -45) / -12 m3 (-1.5) / -0.2 m3 over 60 days,
at 6.7 / 8.0 / 10.8 outer iterations per step. A toy MF6 model, not La Mata.

    python code/tests/diag_lak_mover_leak.py [--ws <folder>]
"""
import argparse
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..')))

NPER = 60
NPOND = 5
SEG = 3                              # stream cells between two ponds
NROW = (NPOND + 1) * (SEG + 1)
NCOL, CELL = 5, 50.0
J = 2                                # the stream's column


def forcing(seed=1):
    rng = np.random.default_rng(seed)
    rain = np.where(rng.random(NPER) < 0.15,
                    rng.uniform(0.005, 0.03, NPER), 0.0)        # m/d
    inflow = 600.0 + 400.0 * np.sin(np.arange(NPER) / 3.0) + 3e4 * rain
    return rain, inflow


def build(ws, dvclose, exe):
    """Five ponds in series, each taking the stream through the mover and
    spilling back through it, FACTOR 1 -- La Mata's on-channel ponds."""
    import flopy
    rain, inflow = forcing()
    name = 'toy'
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=NPER, perioddata=[(1.0, 1, 1.0)] * NPER,
                          time_units='DAYS')
    flopy.mf6.ModflowIms(sim, complexity='COMPLEX', outer_dvclose=dvclose,
                         outer_maximum=500, inner_dvclose=min(dvclose, 1e-4),
                         linear_acceleration='BICGSTAB',
                         csv_outer_output_filerecord='outer.csv')
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True,
                               newtonoptions='NEWTON')
    i, j = np.arange(NROW), np.arange(NCOL)
    top = 100.0 - 0.3 * i[:, None] + 1.0 * np.abs(j[None, :] - J)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=NROW, ncol=NCOL, delr=CELL,
                            delc=CELL, top=top, botm=60.0)
    flopy.mf6.ModflowGwfic(gwf, strt=top - 1.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=0.05)            # La Mata
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.01,   # La Mata
                            transient={0: True})
    flopy.mf6.ModflowGwfrcha(gwf, recharge={p: float(rain[p])
                                            for p in range(NPER)})
    flopy.mf6.ModflowGwfdrn(gwf, stress_period_data={0: [
        [(0, a, b), float(top[a, b]), 100.0] for a in range(NROW)
        for b in range(NCOL) if b != J]})
    pond_rows = [(k + 1) * (SEG + 1) - 1 for k in range(NPOND)]
    segs, row = [], 0
    for pr in pond_rows + [NROW]:
        segs.append(list(range(row, pr)))
        row = pr + 1
    pkg, con, rno = [], [], {}
    for s, rows in enumerate(segs):
        for k, r in enumerate(rows):
            rno[(s, k)] = len(pkg)
            pkg.append([len(pkg), (0, r, J), CELL, 2.0, 0.006,
                        float(top[r, J]) - 0.5, 0.3, 0.05, 0.035, 0, 1.0, 0])
    for s, rows in enumerate(segs):
        for k in range(len(rows)):
            n = rno[(s, k)]
            c = ([rno[(s, k - 1)]] if k > 0 else []) + \
                ([-rno[(s, k + 1)]] if k < len(rows) - 1 else [])
            con.append([n] + c)
            pkg[n][9] = len(c)
    flopy.mf6.ModflowGwfsfr(gwf, pname='sfr', nreaches=len(pkg),
                            packagedata=[tuple(p) for p in pkg],
                            connectiondata=con,
                            perioddata={p: [(0, 'INFLOW', float(inflow[p]))]
                                        for p in range(NPER)},
                            mover=True, save_flows=True,
                            length_conversion=1.0, time_conversion=86400.0,
                            budget_filerecord='%s.sfr.cbc' % name)
    lakes, conns, outs, mvr = [], [], [], []
    for k, pr in enumerate(pond_rows):
        bed = float(top[pr, J])          # the pond sits on its cell
        lakes.append((k, bed + 0.5, 1))
        conns.append((k, 0, (0, pr, J), 'VERTICAL', 0.001, bed, bed, 0.0, 0.0))
        outs.append((k, k, -1, 'MANNING', bed + 0.5, 5.0, 0.035, 0.001))
        mvr.append(('sfr', rno[(k, len(segs[k]) - 1)], 'lak', k, 'FACTOR', 1.0))
        mvr.append(('lak', k, 'sfr', rno[(k + 1, 0)], 'FACTOR', 1.0))
    flopy.mf6.ModflowGwflak(gwf, pname='lak', nlakes=NPOND, noutlets=NPOND,
                            mover=True, save_flows=True, surfdep=0.05,
                            maximum_iterations=200, maximum_stage_change=1e-4,
                            length_conversion=1.0, time_conversion=86400.0,
                            packagedata=lakes, connectiondata=conns,
                            outlets=outs,
                            budget_filerecord='%s.lak.cbc' % name)
    flopy.mf6.ModflowGwfmvr(gwf, maxmvr=len(mvr), maxpackages=2,
                            packages=[('sfr',), ('lak',)], perioddata={0: mvr})
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord='%s.cbc' % name,
                           saverecord=[('BUDGET', 'ALL')])
    sim.write_simulation(silent=True)
    return sim, name


def main(argv=None):
    import flopy
    import pandas as pd
    import mm_paths
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--ws', default=None,
                    help='a scratch folder (default: a temporary one)')
    a = ap.parse_args(argv)
    root = a.ws or tempfile.mkdtemp(prefix='mm_lak_mover_leak_')
    for dv in (0.025, 0.001, 1e-5):
        ws = os.path.join(root, 'dv_%g' % dv)
        os.makedirs(ws, exist_ok=True)
        sim, name = build(ws, dv, mm_paths.MF6_EXE)
        ok, _ = sim.run_simulation(silent=True)
        if not ok:
            print('outer_dvclose %g: the run failed (see %s)' % (dv, ws))
            continue
        cb = flopy.utils.CellBudgetFile(os.path.join(ws, '%s.lak.cbc' % name),
                                        precision='double')
        ext = np.array([r['q'] for r in cb.get_data(text='EXT-OUTFLOW')])
        tom = np.array([r['q'] for r in cb.get_data(text='TO-MVR')])
        it = pd.read_csv(os.path.join(ws, 'outer.csv')).groupby('totim')[
            'nouter'].max()
        print('outer_dvclose %-7g %.1f outer iterations per step'
              % (dv, it.mean()))
        for k in range(NPOND):
            print('   pond %d (%d upstream): spill %7.0f m3, EXT-OUTFLOW %8.2f '
                  'm3, worst day %7.2f m3/d (bound %.1f)'
                  % (k + 1, k, -tom[:, k].sum(), ext[:, k].sum(),
                     ext[:, k].min(), dv * CELL * CELL))
    print('files in %s' % root)


if __name__ == '__main__':
    main()
