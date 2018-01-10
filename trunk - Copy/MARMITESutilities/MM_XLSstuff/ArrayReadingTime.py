#-------------------------------------------------------------------------------
# Name:        module2
# Purpose:
#
# Author:      alf
#
# Created:     14-12-2011
# Copyright:   (c) alf 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import cProfile
import numpy as np
import pstats

def s_drc():
    global s1
    print 'drc'
    s = 0
    for t in range(ndays):
        for row in range(nrow):
            for col in range(ncol):
                s += 1
                s1 = t + row + col
    print '# of cells is %d' % (nlay*s)
    print 'Stupid sum is %d\n' %s1

def s_crd():
    global s1
    print 'crd'
    s = 0
    for col in range(ncol):
        for row in range(nrow):
            for t in range(ndays):
                s += 1
                s1 = t + row + col
    print '# of cells is %d' % (nlay*s)
    print 'Stupid sum is %d\n' %s1

def s_rcd():
    global s1
    print 'rcd'
    s = 0
    for row in range(nrow):
        for col in range(ncol):
            for t in range(ndays):
                s += 1
                s1 = t + row + col
    print '# of cells is %d' % (nlay*s)
    print 'Stupid sum is %d\n' %s1

if __name__ == '__main__':

    ndays = 3650
    nrow  = 50
    ncol  = 150
    nlay  = 2

    array = np.ones((ndays,nrow,ncol, nlay), dtype = np.int)

    print '%d days x %d rows x %d cols x %d layers' % (ndays, nrow, ncol, nlay)

    print '\n#######################################'
    cProfile.run('s_drc()', 'profile.tmp')
    p = pstats.Stats('profile.tmp')
    p.sort_stats('cumulative').print_stats(5)

    print '#######################################'
    cProfile.run('s_crd()', 'profile.tmp')
    p = pstats.Stats('profile.tmp')
    p.sort_stats('cumulative').print_stats(5)

    print '#######################################'
    cProfile.run('s_rcd()', 'profile.tmp')
    p = pstats.Stats('profile.tmp')
    p.sort_stats('cumulative').print_stats(5)

    print '#######################################'
