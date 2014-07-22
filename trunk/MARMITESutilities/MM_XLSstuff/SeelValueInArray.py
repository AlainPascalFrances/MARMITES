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

def find_day():
    print "Seek per day\n"
    for t in range(ndays):
        for row,col in enumerate(array[t,:,:]):
            try:
                print 'row %d, col %d and day %d\n' % (row, list(col).index(searchValue), t)
            except:
                pass

def find_row():
    print "Seek per row\n"
    for row in range(nrow):
        for t,col in enumerate(array[:,row,:]):
            try:
                print 'row %d, col %d and day %d\n' % (row, list(col).index(searchValue), t)
            except:
                pass

def find_col():
    print "Seek per col\n"
    for col in range(ncol):
        for t,row in enumerate(array[:,:,col]):
            try:
                print 'row %d, col %d and day %d\n' % (list(row).index(searchValue), col, t)
            except:
                pass

if __name__ == '__main__':

    nrow = 50
    ncol = 100
    ndays = 500

    array = np.ones((ndays,nrow,ncol), dtype = np.int)
    searchValue = -1
    array[ndays/2,nrow/2,ncol/2] = searchValue

    print '\n#######################################'
    cProfile.run('find_day()', 'profile.tmp')
    p = pstats.Stats('profile.tmp')
    p.sort_stats('cumulative').print_stats(5)

    print '\n#######################################'
    cProfile.run('find_row()', 'profile.tmp')
    p = pstats.Stats('profile.tmp')
    p.sort_stats('cumulative').print_stats(5)

    print '\n#######################################'
    cProfile.run('find_col()', 'profile.tmp')
    p = pstats.Stats('profile.tmp')
    p.sort_stats('cumulative').print_stats(5)

