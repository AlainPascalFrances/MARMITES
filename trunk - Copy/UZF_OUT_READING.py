#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      alf
#
# Created:     14-09-2011
# Copyright:   (c) alf 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
import os
import numpy as np

def main():
    ws = r'c:\temp'
    fn = 'r13c6l2'
    ext1 = '.uzf'
    ext2 = ['58','Pz0','Pz1','Pz2','Pz3','Pz4','Pz5','Pz6','Pz9']
    inputFile_fn = []
    uzf_out_tot = []
    list_uzf_out = [['TIME', 'INF_APL',  'RUNOFF',  'INF_ACT',  'EXF',  'ET_UZ',  'ET_GW',  'deltaSTOR_UZ',  'RCH'],
    ['LAYER',  'TIME',  'HEAD',  'THICK_UZ',  'cINF_APL',  'cINF',  'cRCH',  'totSTOR_UZ',  'deltaSTOR_UZ',  'EXF',  'rINF_APL',  'rINF',  'rRCH',  'rSTOR',  'rEXF']]
    n = 0
    for e in ext2:
        inputFile_fn.append(os.path.join(ws, fn + ext1 + '_' + e))
        print '\n####\n' + inputFile_fn[n]
        if os.path.exists(inputFile_fn[n]):
            fin = open(inputFile_fn[n], 'r')
            line = fin.readline().split()
            line = fin.readline().split()
            list_uzf_out_temp = (line[1:-1])
            uzf_out_tot.append(np.loadtxt(inputFile_fn[n], skiprows = 3))
            for l in range(uzf_out_tot[n].shape[1]):
                if n == 0:
                    if len(list_uzf_out_temp) == len(list_uzf_out[0]):
                        lst = list_uzf_out[0]
                    else:
                        print "Headers don't correspond!"
                        lst = list_uzf_out_temp
                    print lst[l],uzf_out_tot[n][:,l].sum()/67/100/100
                else:
                    if len(list_uzf_out_temp) == len(list_uzf_out[1]):
                        lst = list_uzf_out[1]
                    else:
                        print "Headers don't correspond!"
                        lst = list_uzf_out_temp
                    print lst[l],uzf_out_tot[n][:,l].sum()/100/100
        else:
            print "\nThis file doesn't exist!"
        n += 1

if __name__ == '__main__':
    main()
    print '\nDone!'
