#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      alf
#
# Created:     01-05-2013
# Copyright:   (c) alf 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import os, sys, traceback
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

class clsUTILITIES():
    def __init__(self, fmt = mpl.dates.DateFormatter('%Y-%m-%d %H:%M'), verbose = 1, report_fn = ''):
        self.fmt = fmt
        self.verbose = verbose
        self.report_fn = report_fn

#####################################

    def ErrorExit(self, msg = 'Undefined error.', stdout = None, report = None):
        print '%s\nError description:' % msg
        traceback.print_exc(file=sys.stdout)
        print ('\n##############\nWARNING!\nMARMITES terminated with ERROR!\n%s\n##############' % (mpl.dates.DateFormatter.format_data(self.fmt, mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat()))))
        if self.verbose == 0:
            sys.stdout = stdout
            report.close()
            raise SystemExit('##############\nWARNING!\nMARMITES terminated with ERROR!\nCheck report file:\n%s.\n%s\n##############' % (self.report_fn, mpl.dates.DateFormatter.format_data(self.fmt, mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat()))))
        else:
            raise SystemExit()

#####################################

    def readFile(self, ws, fn):
        inputFile = []
        inputFile_fn = os.path.join(ws, fn)
        if os.path.exists(inputFile_fn):
            fin = open(inputFile_fn, 'r')
        else:
            self.ErrorExit(msg = "File [%s] doesn't exist, verify name and path!"%inputFile_fn)
        line = fin.readline().split()
        delimChar = line[0]
        try:
            for line in fin:
                line_tmp = line.split(delimChar)
                if not line_tmp == []:
                    if (not line_tmp[0] == '') and (not line_tmp[0] == '\n') and (not line_tmp[0].isspace()):
                        inputFile.append(line_tmp[0])
                else:
                    raise NameError('InputFileFormat')
        except NameError:
            self.ErrorExit('Error in file [%s], check format!'%inputFile_fn)
        except:
            self.ErrorExit("Unexpected error in file [%s]\n"%inputFile_fn)
        fin.close()
        del fin
        return inputFile

#####################################

    def which(self, program):
        import os
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file

        return None

#####################################

    def compDATE_INI(self, date, iniMonthHydroYear):
        year = mpl.dates.num2date(date).year
        month = mpl.dates.num2date(date).month
    #    day = mpl.dates.num2date(date).day
        if iniMonthHydroYear == 1:
            iniMonthHydroYear = 12
            year -= 1            
        if month >= iniMonthHydroYear:
            date_ini = mpl.dates.date2num(mpl.dates.datetime.datetime(year,iniMonthHydroYear,1))
        else:
            date_ini = mpl.dates.date2num(mpl.dates.datetime.datetime(year-1,iniMonthHydroYear,1))
        return date_ini, year

#####################################

    def compDATE_END(self, date, iniMonthHydroYear):
        year = mpl.dates.num2date(date).year
        month = mpl.dates.num2date(date).month
    #    day = mpl.dates.num2date(date).day
        if iniMonthHydroYear == 12:
            iniMonthHydroYear = 1
            year += 1            
        if month >= iniMonthHydroYear:
            date_end = mpl.dates.date2num(mpl.dates.datetime.datetime(year+1,iniMonthHydroYear,1))
        else:
            date_end = mpl.dates.date2num(mpl.dates.datetime.datetime(year,iniMonthHydroYear,1))
        return date_end, year

#####################################

    def remappedColorMap(self, cmap, start=0, midpoint=0.5, stop=1.0,
        name='shiftedcmap'):
        '''
        Function to offset the median value of a colormap, and scale the
        remaining color range. Useful for data with a negative minimum and
        positive maximum where you want the middle of the colormap's dynamic
        range to be at zero.
        Input
        -----
        cmap : The matplotlib colormap to be altered
        start : Offset from lowest point in the colormap's range.
        Defaults to 0.0 (no lower ofset). Should be between
        0.0 and 0.5; if your dataset mean is negative you should leave
        this at 0.0, otherwise to (vmax-abs(vmin))/(2*vmax)
        midpoint : The new center of the colormap. Defaults to
        0.5 (no shift). Should be between 0.0 and 1.0; usually the
        optimal value is abs(vmin)/(vmax+abs(vmin))
        stop : Offset from highets point in the colormap's range.
        Defaults to 1.0 (no upper ofset). Should be between
        0.5 and 1.0; if your dataset mean is positive you should leave
        this at 1.0, otherwise to (abs(vmin)-vmax)/(2*abs(vmin))
        http://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
        '''
        cdict = {
            'red': [],
            'green': [],
            'blue': [],
            'alpha': []
        }
    
        # regular index to compute the colors
        reg_index = np.hstack([
            np.linspace(start, 0.5, 128, endpoint=False),
            np.linspace(0.5, stop, 129)
        ])
    
        # shifted index to match the data
        shift_index = np.hstack([
            np.linspace(0.0, midpoint, 128, endpoint=False),
            np.linspace(midpoint, 1.0, 129)
        ])
            
        for ri, si in zip(reg_index, shift_index):
            r, g, b, a = cmap(ri)
            cdict['red'].append((si, r, r))
            cdict['green'].append((si, g, g))
            cdict['blue'].append((si, b, b))
            cdict['alpha'].append((si, a, a))
            
        newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
        plt.register_cmap(cmap=newcmap)
    
        return newcmap

#####################################
if __name__ == "__main__":
    print '\nWARNING!\nStart MARMITES-MODFLOW models using the script startMARMITES_v3.py\n'

#EOF