#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      frances08512
#
# Created:     30-03-2011
# Copyright:   (c) frances08512 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
#import CreateColors
import itertools

def Eg(PE, y0, b, dtwt, dll):
    '''
    Groundwater evaporation, equation 17 of Shah et al 2007
    '''
    if dtwt<=dll:
        Eg_tmp = PE
    else:
        Eg_tmp = PE*(y0 + np.exp(-b*(dtwt-dll)))
    if Eg_tmp/PE < 0.5/100:
        Eg_tmp = 0.0
    return Eg_tmp

paramEg = {'sand'        : {'dll':160.0,'y0':0.000,'b':0.0171, 'ext_d':500.0},
       'loamy sand'      : {'dll':210.0,'y0':0.002,'b':0.0130, 'ext_d':700.0},
       'sandy loam'      : {'dll':300.0,'y0':0.004,'b':0.0065, 'ext_d':1300.0},
       'sandy clay loam' : {'dll':300.0,'y0':0.006,'b':0.0046, 'ext_d':2000.0},
       'sandy clay'      : {'dll':200.0,'y0':0.005,'b':0.0042, 'ext_d':2100.0},
       'loam'            : {'dll':330.0,'y0':0.004,'b':0.0028, 'ext_d':2600.0},
       'silty clay'      : {'dll':370.0,'y0':0.007,'b':0.0046, 'ext_d':3300.0},
       'clay loam'       : {'dll':330.0,'y0':0.008,'b':0.0027, 'ext_d':4000.0},
       'silt loam'       : {'dll':380.0,'y0':0.006,'b':0.0019, 'ext_d':4200.0},
       'silt'            : {'dll':310.0,'y0':0.007,'b':0.0021, 'ext_d':4300.0},
       'silty clay loam' : {'dll':400.0,'y0':0.007,'b':0.0021, 'ext_d':4500.0},
       'clay'            : {'dll':450.0,'y0':0.006,'b':0.0019, 'ext_d':6200.0},
       }

lblspc = 0.05
mkscale = 0.5
#colors_ = CreateColors.main(hi=10, hf=150, numbcolors = (len(paramEg.keys())))
lines = itertools.cycle(['-','--','-.',':','.',',','o','v','^','<','>','1','2','3','4','s','p','*','h','H','+','x','D','d','|','_'])
fig = plt.figure(figsize=(11.7, 8.27))

dtwt = np.asarray(range(-1000,6500,10))
#dtwt_norm = []
GWET = []
lbls = ['sand', 'loam', 'silt', 'clay']
for s, lbl in enumerate(lbls):
    GWET.append([])
    #dtwt_norm.append([])
    y0  = paramEg[lbl]['y0']
    b   = paramEg[lbl]['b']
    dll = paramEg[lbl]['dll']
    for d in range(len(dtwt)):
        #dtwt_norm[s].append(dtwt[d]/paramEg[paramEg[lbl]]['ext_d'])
        GWET[s].append(Eg(1.0, y0, b,dtwt[d],dll))

ax1=fig.add_subplot(2,2,1)
plt.setp(ax1.get_xticklabels(), fontsize=8)
plt.setp(ax1.get_yticklabels(), fontsize=8)
for l, (x, lbl) in enumerate(zip(GWET, lbls)) :
    ax1.plot(x, dtwt/1000.0, lines.next(), color = 'black', label=lbl)
plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)
plt.grid(True)
#plt.xlim((-0.01, 1.01))
plt.ylim(0.0, 3.0)
plt.xlabel(r'$E_g/PE$')
plt.ylabel('$d\ (m)$')
ax1.set_ylim(ax1.get_ylim()[::-1])

#ax5=fig.add_subplot(1,2,2)
#plt.setp(ax5.get_xticklabels(), fontsize=8)
#plt.setp(ax5.get_yticklabels(), fontsize=8)
#for l, (x, y, color, lbl) in enumerate(zip(GWET, dtwt_norm, colors_, paramEg.keys())) :
#    ax5.plot(x, y, '-', color=color, label=lbl)
#plt.legend(paramEg.keys(), loc=0, labelspacing=lblspc, markerscale=mkscale)
#leg = plt.gca().get_legend()
#ltext  = leg.get_texts()
#plt.setp(ltext, fontsize=8 )
#plt.grid(True)
#plt.xlim((-0.05, 1.05))
#plt.ylim((-0.05, 1.05))
#plt.xlabel('Eg/PE')
#plt.ylabel('dtwt/ext_d')
#ax5.set_ylim(ax5.get_ylim()[::-1])

plt.show()
img_path = 'Eg_function'
plt.savefig(img_path,dpi=300)
print 'Plot printed:\n%s' % img_path