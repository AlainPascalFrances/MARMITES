# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 16:59:01 2014

@author: frances.alain@gmail.com
"""

""" Computes the water balance for a certain time span
input: ASCII file with watre fluxes wrtitten by MM
"""
import MARMITESutilities as MMutils
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.sankey import Sankey
import numpy as np
import os

#####################################

path = r'E:\00code_ws\00_TESTS\WB_plot'
fn = r'_0CATCHMENT_ts_20051001_20070930.txt'   # _0CATCHMENT_ts_20051001_20070930    _0CATCHMENT_ts_20061231_20111103    _0CATCHMENT_ts_20080531_20101108   _0CATCHMENT_ts_20080531_20090531
StartMonth = 10
DRN = 1
GHB = 0
NL = 2

fmt_DH = mpl.dates.DateFormatter('%Y-%m-%d')
cUTIL = MMutils.clsUTILITIES(fmt = fmt_DH)

inputFile_fn = os.path.join(path, fn)
if os.path.exists(inputFile_fn):
    DATA = np.loadtxt(inputFile_fn, skiprows = 1, dtype = str, delimiter = ',')
else:
    cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe input file [" + inputFile_fn + "] doesn't exist, verify name and path!")
date = DATA[:,0]
datenum = []
datenum_d = []
DATE = np.zeros(len(date), dtype = float)
for i, d in enumerate(date):
    DATE[i] = mpl.dates.datestr2num(d)

year_lst = []
index = []
if mpl.dates.num2date(DATE[0]).month<StartMonth or (mpl.dates.num2date(DATE[0]).month == StartMonth and mpl.dates.num2date(DATE[0]).day == 1):
    year_lst.append(mpl.dates.num2date(DATE[0]).year)
else:
    year_lst.append(mpl.dates.num2date(DATE[0]).year + 1)
index.append(np.argmax(DATE[:] == mpl.dates.datestr2num('%s-10-01' % year_lst[0])))

if sum(DATE[:] == mpl.dates.datestr2num('%s-10-01' % (year_lst[0]+1)))==0:
    cUTIL.ErrorExit(msg = '\nThe data file does not contain a full hydrological year starting at date 01/%d' % StartMonth)

y = 0    
while DATE[-1] >= mpl.dates.datestr2num('%d-10-01' % (year_lst[y]+1)):
    if DATE[-1] >= mpl.dates.datestr2num('%d-09-30' % (year_lst[y]+2)):
        year_lst.append(year_lst[y]+1)
        index.append(np.argmax(DATE[:] == mpl.dates.datestr2num('%d-10-01' % (year_lst[y]+1))))
        indexend = np.argmax(DATE[:] == mpl.dates.datestr2num('%d-09-30' % (year_lst[y]+2)))
        y += 1
    else:
        break
index.append(indexend)
del indexend
        
print year_lst
print index
print '\nStarting date of hydrological year(s):'
for j in index:
    print mpl.dates.DateFormatter.format_data(fmt_DH, DATE[j]), '%.2f' % float(DATA[j,1])
print 'End date of last hydrological year:'
print mpl.dates.DateFormatter.format_data(fmt_DH, DATE[index[-1]-1]), '%.2f' % float(DATA[index[-1]-1,1])

# compute fluxes for whole modelled period and hydrological years
RF=[365.0*np.sum(np.float16(DATA[:,1]))/len(DATE)]
I=[365.0*np.sum(np.float16(DATA[:,2]))/len(DATE)]
RFe=[365.0*np.sum(np.float16(DATA[:,3]))/len(DATE)]
DSsurf=[365.0*np.sum(np.float16(DATA[:,4]))/len(DATE)]    
Ro=[365.0*np.sum(np.float16(DATA[:,5]))/len(DATE)]
Esurf=[365.0*np.sum(np.float16(DATA[:,6]))/len(DATE)]
DSsoil=[365.0*np.sum(np.float16(DATA[:,7]))/len(DATE)]
EXF=[365.0*np.sum(np.float16(DATA[:,8]))/len(DATE)]
Esoil=[365.0*np.sum(np.float16(DATA[:,9]))/len(DATE)]
Tsoil=[365.0*np.sum(np.float16(DATA[:,10]))/len(DATE)]
ETsoil=[365.0*np.sum(np.float16(DATA[:,11]))/len(DATE)]
Eg=[365.0*np.sum(np.float16(DATA[:,12]))/len(DATE)]
Tg=[365.0*np.sum(np.float16(DATA[:,13]))/len(DATE)]
ETg=[365.0*np.sum(np.float16(DATA[:,14]))/len(DATE)]
Ssurf=[365.0*np.sum(np.float16(DATA[:,15]))/len(DATE)]
Rp=[365.0*np.sum(np.float16(DATA[:,18]))/len(DATE)]    
R=[365.0*np.sum(np.float16(DATA[:,20]))/len(DATE)]
DSu=[Rp[0]-R[0]]    
DSg=[365.0*np.sum(np.float16(DATA[:,25]))/len(DATE)]
DRN=[365.0*np.sum(np.float16(DATA[:,26]))/len(DATE)]
#RF 1, #I 2, #RFe 3, #DeltaSsurf 4, #Ro 5, #Esurf 6, #DeltaSsoil 7, #EXF_g 8, #Esoil 9, #Tsoil 10, #ETsoil 11, #E_g 12, #T_g 13, #ET_g 14, #Ssurf 15, #PE 16, #PT 17, #inf 18, #theta 19, #R 20, #DeltaSg1 25, #DRN1 26, #DeltaSg2 27, #DRN2 28
for k, i in enumerate(index[:-1]):
    RF.append(np.sum(np.float16(DATA[i:(index[k+1]-1),1])))
    I.append(np.sum(np.float16(DATA[i:(index[k+1]-1),2])))
    RFe.append(np.sum(np.float16(DATA[i:(index[k+1]-1),3])))
    DSsurf.append(np.sum(np.float16(DATA[i:(index[k+1]-1),4])))    
    Ro.append(np.sum(np.float16(DATA[i:(index[k+1]-1),5])))
    Esurf.append(np.sum(np.float16(DATA[i:(index[k+1]-1),6])))
    DSsoil.append(np.sum(np.float16(DATA[i:(index[k+1]-1),7])))
    EXF.append(np.sum(np.float16(DATA[i:(index[k+1]-1),8])))
    Esoil.append(np.sum(np.float16(DATA[i:(index[k+1]-1),9])))
    Tsoil.append(np.sum(np.float16(DATA[i:(index[k+1]-1),10])))
    ETsoil.append(np.sum(np.float16(DATA[i:(index[k+1]-1),11])))
    Eg.append(np.sum(np.float16(DATA[i:(index[k+1]-1),12])))
    Tg.append(np.sum(np.float16(DATA[i:(index[k+1]-1),13])))
    ETg.append(np.sum(np.float16(DATA[i:(index[k+1]-1),14])))
    Ssurf.append(np.sum(np.float16(DATA[i:(index[k+1]-1),15])))
    Rp.append(np.sum(np.float16(DATA[i:(index[k+1]-1),18])))    
    R.append(np.sum(np.float16(DATA[i:(index[k+1]-1),20])))
    DSu.append(Rp[k+1]-R[k+1])    
    DSg.append(np.sum(np.float16(DATA[i:(index[k+1]-1),25])))
    DRN.append(np.sum(np.float16(DATA[i:(index[k+1]-1),26])))

# Sankey plots    
fig = []
ax = []
for k in range(len(RF)):
    print '\n########################'
    if k == 0:
        print 'Water balance of the whole modelled period:'
        title = "La Mata catchment - Water balance (units in $mm.y^{-1}$)\nWhole modelled period"
    else:        
        print 'Water balance of hydrological year %d/%d:' % (year_lst[k-1], year_lst[k-1] + 1)
        title = "La Mata catchment - Water balance (units in $mm.y^{-1}$)\nHydrological year %d/%d" % (year_lst[k-1], year_lst[k-1] + 1)
    print '\nMMsurf: RF %s, I %s, RFe %s, Esurf %s, DSsurf %s\nMMsoil: Esoil %s, Tsoil %s, ETsoil %s, Ro %s, Rp %s, DSsoil  %s\nUZF: R %s, EXF  %s, DSu %s\nMF: Eg %s, Tg %s, ETg %s,  DRN %s, DSg %s' % (RF[k], I[k], RFe[k], Esurf[k], DSsurf[k], Esoil[k], Tsoil[k], ETsoil[k], Ro[k], Rp[k], DSsoil[k], R[k], EXF[k], DSu[k], Eg[k], Tg[k], ETg[k], DRN[k], DSg[k])
    fig.append(plt.figure()) # figsize=(8.27, 11.7), dpi = 72) 
    ax.append(fig[k].add_subplot(1, 1, 1, xticks=[], yticks=[], title=title))
    pltsankey = Sankey(ax=ax[k], format='%.1F', scale=1.0/1000.0, offset = 0.25)
    pltsankey.add(label='MMsurf', facecolor='lightblue', trunklength = 1,
               flows=[RF[k], -I[k], -RFe[k], DSsurf[k], -Esurf[k]],
               labels=['$RF$', '$I$', '$RFe$',r'$\Delta S_{surf}$', '$E_{surf}$'],
               orientations=[1,1,-1,0,1],
               pathlengths = [0.15, 0.15, 0.15, 0.15, 0.15])
    pltsankey.add(label='MMsoil', facecolor='khaki', trunklength = 1,
               flows=[RFe[k], -Rp[k], -Esoil[k], -Tsoil[k], DSsoil[k], EXF[k], -Ro[k]],
               labels=[None,'$Rp$','$E_{soil}$','$T_{soil}$',r'$\Delta S_{soil}$', None,'$Ro$'],
               orientations=[1,-1,1,1,0,-1,0],
               pathlengths = [0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15],
               prior=0, connect=(2,0))
    pltsankey.add(label='MF_UZF', facecolor='aliceblue', trunklength = 1,
               flows=[Rp[k], -R[k], DSu[k]],
               labels=[None, '$R$', '$\Delta S_u$'],
               orientations=[1, -1, 0],
               pathlengths = [0.15, 0.15, 0.15],
               prior=1, connect=(1, 0))
    pltsankey.add(label='MF_L1', facecolor='MediumSlateBlue', trunklength = 1,
               flows=[R[k], DSg[k], DRN[k], -EXF[k], ],
               labels=[None, '$\Delta S_g$', '$DRN$', '$EXF$'],
               orientations=[1, 0, 1, -1],
               pathlengths = [0.15, 0.15, 0.15, 0.15],
               prior=2, connect=(1, 0))
    diagrams = pltsankey.finish()
    diagrams[-1].patch.set_hatch('/')
    for d in diagrams:
        for t in d.texts:
            t.set_fontsize(8)
    plt.legend(loc='best')
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=10)
#    # Notice that the explicit connections are handled automatically, but the
#    # implicit ones currently are not.  The lengths of the paths and the trunks
#    # must be adjusted manually, and that is a bit tricky.
    if k == 0:
        print '\nPlot of the whole modelled period done!\n########################'
    else:
        print '\nPlot of hydrological year %d/%d done!\n########################' % (year_lst[k-1], year_lst[k-1] + 1)

print "\nWater fluxes imported from file %s!\n" % inputFile_fn



# EOF