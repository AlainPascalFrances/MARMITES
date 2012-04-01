import numpy as np
import os

nrow = 65
ncol = 60
xllcorner = 739300
yllcorner = 4553050
cellsizeMF = 50
k = '5.0'

pp_in_fn = r'E:\00code_ws\LaMata_new\GIS\PilotsPoints.txt'
if os.path.exists(pp_in_fn):
    fin = open(pp_in_fn, 'r')
pp_out_fn = r'E:\00code_ws\LaMata_new\MF_ws\hk.pts'
if os.path.exists(pp_in_fn):
    fout = open(pp_out_fn, 'w')

fin.readline()
for line in fin:
    line_tmp = line.split()
    x = float(line_tmp[1])
    y = float(line_tmp[2])
    i = nrow - np.ceil((y-yllcorner)/cellsizeMF)
    j = np.ceil((x-xllcorner)/cellsizeMF) - 1
    l = 'pp%s %i %i 1 %s\n' % (line_tmp[0],i,j,k)
    fout.write(l)
fin.close()
fout.close()
print 'Done!'
