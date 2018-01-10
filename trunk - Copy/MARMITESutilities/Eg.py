#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      frances.alain@gmail.com
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
    return Eg_tmp

paramEg = {'sand'            : {'dll':16.0,'y0':0.000,'b':0.171, 'ext_d':50.0},
           'loamy sand'      : {'dll':21.0,'y0':0.002,'b':0.130, 'ext_d':70.0},
           'sandy loam'      : {'dll':30.0,'y0':0.004,'b':0.065, 'ext_d':130.0},
           'sandy clay loam' : {'dll':30.0,'y0':0.006,'b':0.046, 'ext_d':200.0},
           'sandy clay'      : {'dll':20.0,'y0':0.005,'b':0.042, 'ext_d':210.0},
           'loam'            : {'dll':33.0,'y0':0.004,'b':0.028, 'ext_d':265.0},
           'silty clay'      : {'dll':37.0,'y0':0.007,'b':0.046, 'ext_d':335.0},
           'clay loam'       : {'dll':33.0,'y0':0.008,'b':0.027, 'ext_d':405.0},
           'silt loam'       : {'dll':38.0,'y0':0.006,'b':0.019, 'ext_d':420.0},
           'silt'            : {'dll':31.0,'y0':0.007,'b':0.021, 'ext_d':430.0},
           'silty clay loam' : {'dll':40.0,'y0':0.007,'b':0.021, 'ext_d':450.0},
           'clay'            : {'dll':45.0,'y0':0.006,'b':0.019, 'ext_d':620.0},
           'sandy loam field_enrico': {'dll':100.0,'y0':0.000,'b':0.013, 'ext_d':475.0},
           'sandy loam field': {'dll':115.6,'y0':0.004,'b':0.013, 'ext_d':1000.0}
           }
d_pts = [0,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.50,3.00,3.50,4.00,4.50,5.00]
GWET_pts = [1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00,9.01E-01,6.59E-01,4.66E-01,3.37E-01,1.84E-01,1.12E-01,7.59E-02,5.54E-02,4.33E-02,3.56E-02]

lblspc = 0.05
mkscale = 0.5
#colors_ = CreateColors.main(hi=10, hf=150, numbcolors = (len(paramEg.keys())))
lines = itertools.cycle(['-','--','-.',':','--','-.'])
fig = plt.figure(figsize=(11.7, 8.27))

dtwt = np.asarray(range(-100,5000,1))
#dtwt_norm = []
GWET = []
lbls = ['sand', 'loam', 'silt', 'clay', 'sandy loam field'] #, 'sandy loam field_enrico']
for s, lbl in enumerate(lbls):
    GWET.append([])
    #dtwt_norm.append([])
    y0  = paramEg[lbl]['y0']
    b   = paramEg[lbl]['b']
    dll = paramEg[lbl]['dll']
    for d in range(len(dtwt)):
        #dtwt_norm[s].append(dtwt[d]/paramEg[paramEg[lbl]]['ext_d'])
        GWET[s].append(Eg(1.0, y0, b,dtwt[d],dll))
    GWET[s] = np.ma.where(np.asarray(GWET[s])>0.001,GWET[s],-50.0)

ax1=fig.add_subplot(2,2,1)
plt.setp(ax1.get_xticklabels(), fontsize=8)
plt.setp(ax1.get_yticklabels(), fontsize=8)
for l, (x, lbl, col) in enumerate(zip(GWET, lbls, ['black','black','black','black','grey','grey'])) :
    ax1.plot(x, dtwt/100.0, lines.next(), color = col, label=lbl)
ax1.plot(GWET_pts, d_pts, 'o', color = 'grey', label = 'field/HYDRUS data')
plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)
plt.grid(True)
plt.xlim((0.0, 1.0))
plt.ylim(0.0, 5.00)
plt.xlabel(r'$E_g/PE_g$')
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

#plt.show()
img_path = 'Eg_function_201603'
plt.savefig(img_path,dpi=300)
print 'Plot printed:\n%s' % img_path