# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:25:53 2014

@author: apf
"""

import matplotlib.pyplot as plt
import numpy as np
#import CreateColors
import itertools

# INPUT
phi = 0.40
theta_wp = 0.1
PT = 6.0
kT_min = 0.2
kT_max = 0.8
kT_n = [0.5,1.0,2.0]
daynum= 120

# plot parameters (legend)
lblspc = 0.05
mkscale = 0.5
handletextpad = 0.25
borderaxespad = 0.5
fig = plt.figure(figsize=(11.7, 8.27))

incr = 0.001
theta = np.arange(theta_wp,phi, incr)
ax1=fig.add_subplot(2,2,1)
plt.setp(ax1.get_xticklabels(), fontsize=8)
plt.setp(ax1.get_yticklabels(), fontsize=8)
for l, n in enumerate(kT_n):
    kT = kT_min+(kT_max-kT_min)*np.power(1-np.power((phi-theta)/(phi - theta_wp),n),1/n)
    kT[-1] = kT_max
    ax1.plot(theta, kT, color = 'black')
    x_lbl = theta_wp + (phi-theta_wp)/2.0
    y_lbl = kT_min+(kT_max-kT_min)*np.power(1-np.power((phi-x_lbl)/(phi - theta_wp),n),1/n)
    ax1.annotate('$n = %.1f$' % n, (x_lbl, y_lbl), horizontalalignment='right', verticalalignment='bottom', fontsize = 10, color = 'black')
plt.grid(True)
plt.xlim((theta_wp, phi))
plt.ylim((kT_min,kT_max))
plt.xlabel(r'$\theta\ (\%)$')
plt.ylabel(r'$k_T$')

days = np.arange(1,daynum + 1,1)
Tsoil = np.zeros([2,daynum], float)
Tg    = np.zeros([2,daynum], float)
T     = np.zeros([2,daynum], float)
theta = np.zeros([2,daynum], float)
for i, n in enumerate([0.5, 2.0]):
    for d in range(len(days)):
        if d == 0:
            theta_tmp = phi*1000.0
        else:
            theta_tmp = theta[i,d-1]
        Tsoil[i,d] = PT * (theta_tmp/1000.0 - theta_wp)/(phi - theta_wp)
        theta[i,d] = theta_tmp - Tsoil[i,d]
        kT = kT_min+(kT_max-kT_min)*np.power(1-np.power((phi-theta[i,d]/1000.0)/(phi - theta_wp),n),1/n)
        Tg[i,d] = Tsoil[i,d]*(1.0/kT-1.0)
        T[i,d] = Tsoil[i,d] + Tg[i,d]
        if T[i,d] > PT:
            Tg[i,d] = PT - Tsoil[i,d]
            T[i,d]  = PT                

lines = itertools.cycle(['-','--','-.','-'])
ax2=fig.add_subplot(4,2,2)
plt.setp(ax2.get_xticklabels(), fontsize=8)
plt.setp(ax2.get_yticklabels(), fontsize=8)
ax2.plot(days, T[0], lines.next(), color = 'black', label = r'$T$')
ax2.plot(days, Tsoil[0], lines.next(), color = 'black', label = r'$T_{soil}$')
ax2.plot(days, Tg[0], lines.next(), color = 'black', label = r'$T_g$')
plt.grid(True)
plt.ylabel(r'$(mm)$')
plt.ylim((0,PT+0.5))
ax2a = ax2.twinx()
plt.setp(ax2a.get_yticklabels(), fontsize=8, color = 'darkgray')
ax2a.plot(days, theta[0]/1000.0, lines.next(), color = 'darkgray', label = r'$\theta$')
ax2a.text(.15, .025, '$n = 0.5$', horizontalalignment='left', verticalalignment='bottom', transform=ax2a.transAxes, bbox=dict(facecolor='white'), fontsize = 10)
plt.ylim((theta_wp,phi+.025))
plt.ylabel(r'$(\%)$', color = 'darkgray')

lines = itertools.cycle(['-','--','-.','-'])
ax3=fig.add_subplot(4,2,4)
plt.setp(ax3.get_xticklabels(), fontsize=8)
plt.setp(ax3.get_yticklabels(), fontsize=8)
plt3_0 = ax3.plot(days, T[1], lines.next(), color = 'black', label = r'$T$')
plt3_1 = ax3.plot(days, Tsoil[1], lines.next(), color = 'black', label = r'$T_{soil}$')
plt3_2 = ax3.plot(days, Tg[1], lines.next(), color = 'black', label = r'$T_g$')
plt.grid(True)
plt.ylabel(r'$(mm)$')
plt.xlabel(r'$days$')
plt.ylim((0,PT+0.5))
ax3a = ax3.twinx()
plt.setp(ax3a.get_yticklabels(), fontsize=8, color = 'darkgray')
plt3a = ax3a.plot(days, theta[1]/1000.0, lines.next(), color = 'darkgray', label = r'$\theta$')
plt.ylim((theta_wp,phi+.025))
plt.ylabel(r'$(\%)$', color = 'darkgray')
plts =  plt3_0 + plt3_1 + plt3_2 + plt3a
labs = [l.get_label() for l in plts]
plt.legend(plts, labs, loc=1, labelspacing=lblspc, markerscale=mkscale, handletextpad = handletextpad, borderaxespad = borderaxespad, ncol = 2, columnspacing = 0.05)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)
ax3a.text(.15, .025, '$n = 2.0$', horizontalalignment='left', verticalalignment='bottom', transform=ax3a.transAxes, bbox=dict(facecolor='white'), fontsize = 10)

plt.show()
img_path = 'kT_function'
plt.savefig(img_path,dpi=300)
print 'Plot printed:\n%s' % img_path