# -*- coding: utf-8 -*-
"""MARMITES shared utilities: error handling, input-file reading, dates,
colormap helpers.

Phase-1 port (2026): Python 3.12, typed exceptions (MarmitesError),
no module globals, `plt.register_cmap` (removed in matplotlib 3.9)
replaced by `matplotlib.colormaps.register`.
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import datetime
import os
import shutil
import sys
import traceback

import matplotlib as mpl
import matplotlib.dates  # noqa: F401  (registers mpl.dates used across MARMITES)
import numpy as np


class MarmitesError(Exception):
    """Fatal MARMITES error carrying a user-facing message."""


class clsUTILITIES:
    def __init__(self, fmt=None, verbose=1, report_fn=''):
        # default formatter built at call time (avoids import-time evaluation)
        self.fmt = fmt if fmt is not None else mpl.dates.DateFormatter('%Y-%m-%d %H:%M')
        self.verbose = verbose
        self.report_fn = report_fn

    # ------------------------------------------------------------------ #

    def ErrorExit(self, msg='Undefined error.', stdout=None, report=None):
        """Print context and abort the run by raising MarmitesError.

        Kept for API compatibility with the legacy call sites; new code
        should raise MarmitesError directly.
        """
        print('%s\nError description:' % msg)
        traceback.print_exc(file=sys.stdout)
        stamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
        print('\n##############\nWARNING!\nMARMITES terminated with ERROR!\n%s\n##############' % stamp)
        if self.verbose == 0 and stdout is not None:
            sys.stdout = stdout
            if report is not None:
                report.close()
            msg = '%s\nCheck report file:\n%s' % (msg, self.report_fn)
        raise MarmitesError(msg)

    # ------------------------------------------------------------------ #

    def readFile(self, ws, fn):
        """Read a MARMITES sequential input file.

        The first character of the first line defines the comment
        delimiter; for each subsequent line, the content before the first
        delimiter is returned (blank lines skipped).
        """
        inputFile = []
        inputFile_fn = os.path.join(ws, fn)
        if not os.path.exists(inputFile_fn):
            self.ErrorExit(msg="File [%s] doesn't exist, verify name and path!" % inputFile_fn)
        with open(inputFile_fn) as fin:
            line = fin.readline().split()
            if not line:
                self.ErrorExit('Error in file [%s]: empty first line, expected comment delimiter!' % inputFile_fn)
            delimChar = line[0]
            try:
                for line in fin:
                    line_tmp = line.split(delimChar)
                    if not line_tmp:
                        raise MarmitesError('Error in file [%s], check format!' % inputFile_fn)
                    if line_tmp[0] and line_tmp[0] != '\n' and not line_tmp[0].isspace():
                        inputFile.append(line_tmp[0])
            except MarmitesError:
                raise
            except Exception:
                self.ErrorExit('Unexpected error in file [%s]\n' % inputFile_fn)
        return inputFile

    # ------------------------------------------------------------------ #

    @staticmethod
    def which(program):
        """Locate an executable on PATH (cross-platform)."""
        return shutil.which(program)

    # ------------------------------------------------------------------ #

    @staticmethod
    def compDATE_INI(date, iniMonthHydroYear):
        year = mpl.dates.num2date(date).year
        month = mpl.dates.num2date(date).month
        if iniMonthHydroYear == 1:
            iniMonthHydroYear = 12
            year -= 1
        if month >= iniMonthHydroYear:
            date_ini = mpl.dates.date2num(datetime.datetime(year, iniMonthHydroYear, 1))
        else:
            date_ini = mpl.dates.date2num(datetime.datetime(year - 1, iniMonthHydroYear, 1))
        return date_ini, year

    @staticmethod
    def compDATE_END(date, iniMonthHydroYear):
        year = mpl.dates.num2date(date).year
        month = mpl.dates.num2date(date).month
        if iniMonthHydroYear == 12:
            iniMonthHydroYear = 1
            year += 1
        if month >= iniMonthHydroYear:
            date_end = mpl.dates.date2num(datetime.datetime(year + 1, iniMonthHydroYear, 1))
        else:
            date_end = mpl.dates.date2num(datetime.datetime(year, iniMonthHydroYear, 1))
        return date_end, year

    # ------------------------------------------------------------------ #

    @staticmethod
    def remappedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
        """Offset the midpoint of a colormap (for diverging data with
        asymmetric vmin/vmax). See matplotlib SO question 7404116."""
        cdict = {'red': [], 'green': [], 'blue': [], 'alpha': []}
        reg_index = np.hstack([
            np.linspace(start, 0.5, 128, endpoint=False),
            np.linspace(0.5, stop, 129),
        ])
        shift_index = np.hstack([
            np.linspace(0.0, midpoint, 128, endpoint=False),
            np.linspace(midpoint, 1.0, 129),
        ])
        for ri, si in zip(reg_index, shift_index):
            r, g, b, a = cmap(ri)
            cdict['red'].append((si, r, r))
            cdict['green'].append((si, g, g))
            cdict['blue'].append((si, b, b))
            cdict['alpha'].append((si, a, a))
        newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
        try:
            mpl.colormaps.register(newcmap, force=True)
        except AttributeError:  # matplotlib < 3.5 fallback
            mpl.cm.register_cmap(cmap=newcmap)
        return newcmap


if __name__ == '__main__':
    print('\nWARNING!\nStart MARMITES-MODFLOW models using the script startMARMITES_v3.py\n')

# EOF
