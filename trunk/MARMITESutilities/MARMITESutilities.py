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

class clsUTILITIES():
    def __init__(self, fmt = mpl.dates.DateFormatter('%Y-%m-%d %H:%M'), verbose = 1, s = '', report = None, report_fn = ''):
        self.fmt = fmt
        self.verbose = verbose
        self.s = s
        self.report = report
        self.report_fn = report_fn

#####################################

    def ErrorExit(self, msg = 'Undefined error.'):
        print '%s\nError description:' % msg
        traceback.print_exc(file=sys.stdout)
        print ('\n##############\nWARNING!\nMARMITES terminated with ERROR!\n%s\n##############' % (mpl.dates.DateFormatter.format_data(self.fmt, mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat()))))
        if self.verbose == 0:
            sys.stdout = self.s
            self.report.close()
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
                    if (not line_tmp[0] == '') and (not line_tmp[0] == '\n'):
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
#EOF