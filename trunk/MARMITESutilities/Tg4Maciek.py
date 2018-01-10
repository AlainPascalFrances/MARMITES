# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:25:53 2014

@author: apf
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
#import CreateColors
import itertools

# 0 is Tg/PTg increase with theta, 1 is decrease
version = 0

# INPUT
phi = 0.35
theta_wp = 0.05
ThickSoil = 1000.0
PT = 3.5
kT_min =  0.1 #[0.1, 0.8, 0.1]  # [0.1, 0.8, 0.1]
kT_max =  0.9 #[0.9, 0.9, 0.2]  # [0.9, 0.9, 0.2]
daynum= 350

# plot parameters (legend)
lblspc = 0.05
mkscale = 0.15
handletextpad = 0.15
borderpad = 0.05
borderaxespad = 0.5
n_x = 0.97
n_y = 0.500

incr = 0.001

days = np.arange(1,daynum + 1,1)
Tsoil = np.zeros([2,daynum], float)
Tg    = np.zeros([2,daynum], float)
T     = np.zeros([2,daynum], float)
PTg   = np.zeros([2,daynum], float)
theta = np.zeros([2,daynum], float)
PTarray = np.ones([daynum], float)*PT
for i, c in enumerate([2.0]):
    for d in range(len(days)):
        if d == 0:
            theta_tmp = phi*ThickSoil
        else:
            theta_tmp = theta[i,d-1]
        Tsoil[i,d] = PT * (theta_tmp/ThickSoil - theta_wp)/(phi - theta_wp)
        theta[i,d] = theta_tmp - Tsoil[i,d]
        PTg[i,d] = PT-Tsoil[i,d]
        if version == 0:
            kT = kT_min+(kT_max-kT_min)*np.power(1.0-np.power((phi-theta[i,d]/ThickSoil)/(phi - theta_wp),c),1.0/c)
        else:
            kT = kT_max - (kT_max-kT_min)*np.power(1.0-np.power((phi-theta[i,d]/ThickSoil)/(phi - theta_wp),c),1.0/c)
        Tg[i,d] = PTg[i,d]*kT
        T[i,d] = Tsoil[i,d] + Tg[i,d]

lines = itertools.cycle(['-','--','-.','-'])

fig = plt.figure(figsize=(11.7, 8.27))
ax2=fig.add_subplot(2,2,1)
plt.setp(ax2.get_xticklabels(), fontsize=8)
plt.setp(ax2.get_yticklabels(), fontsize=8)
plt2_0 = ax2.plot(days, T[0], lines.next(), color = 'black', label = r'$T$')
plt2_1 = ax2.plot(days, Tsoil[0], lines.next(), color = 'black', label = r'$T_{soil}$')
plt2_2 = ax2.plot(days, Tg[0], lines.next(), color = 'black', label = r'$T_g$')
plt2_3 = ax2.plot(days, PTarray, '-', color = 'red', label = r'$PT$')
plt2_4 = ax2.plot(days, PTg[0], '-.', color = 'red', label = r'$PT_g$')
plt.grid(True)
plt.ylabel(r'$mm$')
plt.xlabel(r'$day$')
plt.ylim((0,PT+0.5))
ax2a = ax2.twinx()
plt.setp(ax2a.get_yticklabels(), fontsize=8, color = 'darkgray')
plt2a = ax2a.plot(days, theta[0]/ThickSoil, lines.next(), color = 'darkgray', label = r'$\theta$')
plt.ylim((theta_wp-.05,phi+.05))
plt.ylabel(r'$(\%)$', color = 'darkgray')
plts =  plt2_3 + plt2_4 + plt2_0 + plt2_1 + plt2_2 + plt2a
labs = [l.get_label() for l in plts]
plt.legend(plts, labs, loc=5, labelspacing=lblspc, markerscale=mkscale, handletextpad = handletextpad, borderaxespad = borderaxespad, borderpad = borderpad, ncol = 3, columnspacing = 0.05)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)

plt.subplots_adjust(left=0.075, bottom=0.05, right=0.95, top=0.95, wspace=0.25, hspace=0.3)    
plt.show()
img_path = 'Tg_function%d'%version
plt.savefig(img_path,dpi=300)
print 'Plot printed:\n%s' % img_path