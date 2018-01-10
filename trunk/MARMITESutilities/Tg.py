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
kT_min_lst =  [0.1, 0.8, 0.1]  # [0.1, 0.8, 0.1]
kT_max_lst =  [0.9, 0.9, 0.2]  # [0.9, 0.9, 0.2]
titlesuffix_lst = ['intermediate phreatophytic behavior', 'major phreatophytic behavior', 'minor phreatophytic behavior']   # ['intermediate'] #
plt_numb_lst = ['b)', 'c)', 'd)']
pos1 = [2,5,6]
pos2 = [4,7,8]
Tg_a_lst = [0.5,1.0,2.0]
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


kT_min = kT_min_lst[0]
kT_max = kT_max_lst[0]
fig = plt.figure(figsize=(11.7, 8.27))
ax1=fig.add_subplot(2,2,1)
plt.setp(ax1.get_xticklabels(), fontsize=8)
plt.setp(ax1.get_yticklabels(), fontsize=8) 
plt.title(r'a) - $T_g/PT_g$ function of $\theta $', fontsize = 12)
theta = np.arange(theta_wp,phi+incr, incr)
for l, a in enumerate(Tg_a_lst):
    if version == 0:
        kT =    kT_min+(kT_max-kT_min)*np.power(np.abs(np.power((phi-theta)/(phi - theta_wp),a)-1.0),1.0/a)
    else:
    # inverse relationship with soil moisture
        kT =    kT_max - (kT_max-kT_min)*np.power(np.abs(np.power((phi-theta)/(phi - theta_wp),a)-1.0),1.0/a)
    kT[-1] = kT_max
    ax1.plot(theta, kT, color = 'black')
    x_lbl = theta_wp + (phi-theta_wp)/2.0
    if version == 0:
        y_lbl = kT_min+(kT_max-kT_min)*np.power(np.abs(np.power((phi-x_lbl)/(phi - theta_wp),a)-1.0),1.0/a)
    else:
        # inverse relationship with soil moisture
        y_lbl = kT_max - (kT_max-kT_min)*np.power(np.abs(np.power((phi-x_lbl)/(phi - theta_wp),a)-1.0),1.0/a)
    ax1.annotate('$a = %.1f$' % a, (x_lbl, y_lbl), horizontalalignment='right', verticalalignment='bottom', fontsize = 10, color = 'black')
plt.grid(True)
plt.xlim((theta_wp, phi))
plt.ylim((kT_min,kT_max))
plt.xlabel(r'$\theta$ (%)')
plt.ylabel(r'$T_g/PT_g$')
ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
ax1.text(theta_wp-0.03, kT_min, '${k_T}_{min}$', horizontalalignment='center', verticalalignment='center', fontsize = 10)
ax1.text(theta_wp-0.03, kT_max, '${k_T}_{max}$', horizontalalignment='center', verticalalignment='center', fontsize = 10)

for j, (kT_min, kT_max, titlesuffix, plt_numb) in enumerate(zip(kT_min_lst, kT_max_lst, titlesuffix_lst, plt_numb_lst)):    
    days = np.arange(1,daynum + 1,1)
    Tsoil = np.zeros([2,daynum], float)
    Tg    = np.zeros([2,daynum], float)
    T     = np.zeros([2,daynum], float)
    PTg   = np.zeros([2,daynum], float)
    theta = np.zeros([2,daynum], float)
    PTarray = np.ones([daynum], float)*PT
    for i, a in enumerate([0.5, 2.0]):
        for d in range(len(days)):
            if d == 0:
                theta_tmp = phi*ThickSoil
            else:
                theta_tmp = theta[i,d-1]
            Tsoil[i,d] = PT * (theta_tmp/ThickSoil - theta_wp)/(phi - theta_wp)
            theta[i,d] = theta_tmp - Tsoil[i,d]
            PTg[i,d] = PT-Tsoil[i,d]
            if version == 0:
                kT = kT_min+(kT_max-kT_min)*np.power(1.0-np.power((phi-theta[i,d]/ThickSoil)/(phi - theta_wp),a),1.0/a)
            else:
                kT = kT_max - (kT_max-kT_min)*np.power(1.0-np.power((phi-theta[i,d]/ThickSoil)/(phi - theta_wp),a),1.0/a)
            Tg[i,d] = PTg[i,d]*kT
            T[i,d] = Tsoil[i,d] + Tg[i,d]
    
    lines = itertools.cycle(['-','--','-.','-'])
    ax2=fig.add_subplot(4,2,pos1[j])
    plt.title('%s - %s' % (plt_numb,titlesuffix), fontsize = 12)
    plt.setp(ax2.get_xticklabels(), fontsize=8)
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    plt2_0 = ax2.plot(days, T[0], lines.next(), color = 'black', label = r'$T$')
    plt2_1 = ax2.plot(days, Tsoil[0], lines.next(), color = 'black', label = r'$T_{soil}$')
    plt2_2 = ax2.plot(days, Tg[0], lines.next(), color = 'black', label = r'$T_g$')
    plt2_3 = ax2.plot(days, PTarray, '-', color = 'red', label = r'$PT$')
    plt2_4 = ax2.plot(days, PTg[0], '-.', color = 'red', label = r'$PT_g$')
    plt.grid(True)
    plt.ylabel('mm')
    plt.xlabel('day')
    plt.ylim((0,PT+0.5))
    ax2a = ax2.twinx()
    plt.setp(ax2a.get_yticklabels(), fontsize=8, color = 'darkgray')
    plt2a = ax2a.plot(days, theta[0]/ThickSoil, lines.next(), color = 'darkgray', label = r'$\theta$')
    plt.ylim((theta_wp-.05,phi+.05))
    plt.ylabel(r'$(\%)$', color = 'darkgray')
    if j == 0:
        plts =  plt2_3 + plt2_4 + plt2_0 + plt2_1 + plt2_2 + plt2a
        labs = [l.get_label() for l in plts]
        plt.legend(plts, labs, loc=0, labelspacing=lblspc, markerscale=mkscale, handletextpad = handletextpad, borderaxespad = borderaxespad, borderpad = borderpad, ncol = 3, columnspacing = 0.05)
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize=10)
    ax2a.text(n_x, n_y, r'$a$ = %.1f, ${k_T}_{min}$ = %.1f, ${k_T}_{max}$ = %.1f' % (0.5, kT_min, kT_max), horizontalalignment='right', verticalalignment='center', transform=ax2a.transAxes, bbox=dict(facecolor='white'), fontsize = 10)
    
    lines = itertools.cycle(['-','--','-.','-'])
    ax3=fig.add_subplot(4,2,pos2[j])
    plt.setp(ax3.get_xticklabels(), fontsize=8)
    plt.setp(ax3.get_yticklabels(), fontsize=8)
    plt3_0 = ax3.plot(days, T[1], lines.next(), color = 'black', label = r'$T$')
    plt3_1 = ax3.plot(days, Tsoil[1], lines.next(), color = 'black', label = r'$T_{soil}$')
    plt3_2 = ax3.plot(days, Tg[1], lines.next(), color = 'black', label = r'$T_g$')
    plt3_3 = ax3.plot(days, PTarray, '-', color = 'red', label = r'$PT$')
    plt3_4 = ax3.plot(days, PTg[1], '-.', color = 'red', label = r'$PT_g$')    
    plt.grid(True)
    plt.ylabel('mm')
    plt.ylim((0,PT+0.5))
    ax3a = ax3.twinx()
    plt.setp(ax3a.get_yticklabels(), fontsize=8, color = 'darkgray')
    plt3a = ax3a.plot(days, theta[1]/ThickSoil, lines.next(), color = 'darkgray', label = r'$\theta$')
    plt.ylim((theta_wp-.05,phi+.05))
    plt.ylabel('(%)', color = 'darkgray')
    #plts =  plt3_0 + plt3_1 + plt3_2 + plt3a
    #labs = [l.get_label() for l in plts]
    #plt.legend(plts, labs, loc=1, labelspacing=lblspc, markerscale=mkscale, handletextpad = handletextpad, borderaxespad = borderaxespad, ncol = 2, columnspacing = 0.05)
    #leg = plt.gca().get_legend()
    #ltext  = leg.get_texts()
    #plt.setp(ltext, fontsize=10)
    ax3a.text(n_x, n_y, r'$a$ = %.1f, ${k_T}_{min}$ = %.1f, ${k_T}_{max}$ = %.1f' % (2.0, kT_min, kT_max), horizontalalignment='right', verticalalignment='center', transform=ax3a.transAxes, bbox=dict(facecolor='white'), fontsize = 10)

plt.subplots_adjust(left=0.075, bottom=0.05, right=0.95, top=0.95, wspace=0.25, hspace=0.3)    
plt.show()
img_path = 'Tg_function%d'%version
plt.savefig(img_path,dpi=300)
print 'Plot printed:\n%s' % img_path