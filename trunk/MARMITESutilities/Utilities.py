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
#EOF