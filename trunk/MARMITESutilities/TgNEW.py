# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:25:53 2014

@author: apf
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
#import CreateColors
import itertools, datetime
mpl.rcParams['mathtext.fontset'] = 'stix'
#matplotlib.rcParams['font.family'] = 'STIXGeneral'
mpl.pyplot.legend(r'ABC123 vs $\mathrm{ABC123}^{123}$')

# INPUT
phi = 0.45
theta_fc = 0.35
theta_wp = 0.05
ThickSoil = 1000.0
PT = 4.0
titlesuffix_lst = ['major phreatophytic behavior', 'intermediate phreatophytic behavior', 'minor phreatophytic behavior']
plt_numb_lst = ['b)', 'c)', 'd)']
pos = [2,4,6]
daynum= 350
#veg param
f_lst = [0.10, 0.217336864608619, 0.30]
s_inv_lst = [20, 28.26674282, 35]
kTmax_lst = [0.9, 0.712676596, 0.35]
kTmin_lst = [0.6, 0.19772923, 0.05]
# kT obs
months = ['June', 'August', 'September']
theta_obs = [0.348,0.125,0.060]
kT_obs = [0.7,0.233,0.205]

cT = 'royalblue'
cTs = 'blue'
cTg = 'deepskyblue'
cPTa = 'red'
cPT = 'red'
cSM = 'green'

# plot parameters (legend)
lblspc = 0.05
mkscale = 0.15
handletextpad = 0.15
borderpad = 0.05
borderaxespad = 0.5
n_x = 0.97
n_y = 0.500

incr = 0.001

s_lst = []
for s_inv in s_inv_lst:
    s_lst.append(1/s_inv)

fig = plt.figure(figsize=(11.7, 8.27))
ax1=fig.add_subplot(1,2,1)
plt.setp(ax1.get_xticklabels(), fontsize=8)
plt.setp(ax1.get_yticklabels(), fontsize=8) 
plt.title(r'a) - $T_g/PT_g$ function of $\theta $', fontsize = 12)
theta = np.arange(theta_wp,phi+incr, incr)
linesTg = itertools.cycle(['--','-','-.'])
for l, (s, s_inv, f, kT_max, kT_min, plt_numb) in enumerate(zip(s_lst, s_inv_lst, f_lst, kTmax_lst, kTmin_lst, plt_numb_lst)):
    kT = kT_max - (kT_max-kT_min)/(1 + np.exp((theta - f)/s))
    ax1.plot(theta, kT, linesTg.__next__(), color = cTg)
    ax1.annotate('$1/s = %.1f,\ f = %.2f$\n$k_{T_{max}} = %.1f,\ k_{T_{min}} = %.1f$\n$(T_g\ parameters\ of\ figure\ %s$' % (s_inv, f, kT_max, kT_min, plt_numb), (f, kT_max - (kT_max-kT_min)/(1 + np.exp((f - f)/s))), horizontalalignment='left', verticalalignment='bottom', bbox=dict(facecolor='white'), fontsize = 10, color = 'black')
plt_leg = ax1.plot(theta_obs,kT_obs, 'o', label = r'Observed $k_T$ (Reyes-Acosta, 2015)')
plt.grid(True)
plt.xlim((theta_wp, phi))
plt.ylim((0.0,1.0))
plt.xlabel(r'$\theta$ (%)')
plt.ylabel(r'$T_g/PT_g$')  #, color = cTg)
#"plt.setp(ax1.get_yticklabels(), color = cTg)
ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
#ax1.text(theta_wp-0.015, kT_min, '${k_T}_{min}$', horizontalalignment='center', verticalalignment='center', fontsize = 10)
#ax1.text(phi+0.015, kT_max, '${k_T}_{max}$', horizontalalignment='center', verticalalignment='center', fontsize = 10)
ax1.text(phi, -0.0375, r'$\phi$', horizontalalignment='center', verticalalignment='center', fontsize = 10)
ax1.text(theta_fc, -0.0375, r'$\theta_{fc}$', horizontalalignment='center', verticalalignment='center', fontsize = 10)
ax1.text(theta_wp, -0.0375, r'$\theta_{wp}$', horizontalalignment='center', verticalalignment='center', fontsize = 10)
labs = [l.get_label() for l in plt_leg]
plt.legend(plt_leg, labs, loc=1, labelspacing=lblspc, markerscale=1, handletextpad=handletextpad, borderaxespad=borderaxespad, borderpad=borderpad)
leg = plt.gca().get_legend()
ltext = leg.get_texts()
#plt.setp(ltext, fontsize=10)

for j, (s, f, kT_max, kT_min, titlesuffix, plt_numb) in enumerate(zip(s_lst, f_lst, kTmax_lst, kTmin_lst, titlesuffix_lst, plt_numb_lst)):
    days = np.arange(1,daynum + 1,1)
    Tsoil = np.zeros(daynum, float)
    Tg    = np.zeros(daynum, float)
    T     = np.zeros(daynum, float)
    PTg   = np.zeros(daynum, float)
    theta = np.zeros(daynum, float)
    PTarray = np.ones([daynum], float)*PT
    for d in range(len(days)):
        if d == 0:
            theta_tmp = phi*ThickSoil
        else:
            theta_tmp = theta[d-1]
        Tsoil[d] = PT * (theta_tmp/ThickSoil - theta_wp)/(phi - theta_wp)
        theta[d] = theta_tmp - Tsoil[d]
        PTg[d] = PT-Tsoil[d]
        kT = kT_max - (kT_max-kT_min)/(1 + np.exp((theta[d]/ThickSoil - f)/s))
        Tg[d] = PTg[d]*kT
        T[d] = Tsoil[d] + Tg[d]
    lines = itertools.cycle(['-','--','-.','-'])
    ax2=fig.add_subplot(3,2,pos[j])
    plt.title('%s - %s' % (plt_numb,titlesuffix), fontsize = 12)
    plt.setp(ax2.get_xticklabels(), fontsize=8)
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    plt2_3 = ax2.plot(days, PTarray, '-', color = cPTa, label = r'$PT$')
    plt2_0 = ax2.plot(days, T, lines.__next__(), color = cT, label = r'$T$')
    plt2_1 = ax2.plot(days, Tsoil, lines.__next__(), color = cTs, label = r'$T_{soil}$')
    plt2_2 = ax2.plot(days, Tg, linesTg.__next__(), color = cTg, label = r'$T_g\ (fig.\ %s$'% plt_numb)
    plt2_4 = ax2.plot(days, PTg, '-.', color = cPT, label = r'$PT_g$')
    plt.grid(True)
    plt.ylabel('mm')
    plt.xlabel('day')
    plt.ylim((0,PT+0.1))
    ax2a = ax2.twinx()
    plt.setp(ax2a.get_yticklabels(), fontsize=8, color = cSM)
    plt2a = ax2a.plot(days, theta/ThickSoil, lines.__next__(), color = cSM, label = r'$\theta$')
    plt.ylim((theta_wp-.05,phi+.01))
    plt.ylabel(r'$(\%)$', color = cSM)
    plt.xlim((0, daynum))
    # ax2a.text(n_x, n_y, '$1/s = %.1f,\ f = %.2f$\n$k_{T_{max}} = %.1f,\ k_{T_{min}} = %.1f$' % (s, f, kT_max, kT_min), horizontalalignment='right', verticalalignment='center', transform=ax2a.transAxes, bbox=dict(facecolor='white'), fontsize = 10)
    if j == 0:
        plts =  plt2_3 + plt2_4 + plt2_0 + plt2_1 + plt2_2
    elif j==1:
        plts += (plt2_2)
    elif j==2:
        plts += (plt2_2+ plt2a)
labs = [l.get_label() for l in plts]
plt.legend(plts, labs, loc='center left', bbox_to_anchor=(-0.25, 2.5), labelspacing=lblspc, markerscale=mkscale, handletextpad = handletextpad, borderaxespad = borderaxespad, borderpad = borderpad, ncol = 1) #, columnspacing = 0.05)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)


plt.subplots_adjust(left=0.075, bottom=0.05, right=0.95, top=0.95, wspace=0.25, hspace=0.3)    
#plt.show()
img_path = 'Tg_functionNEW_%s'% (datetime.date.today().strftime("%Y%m%d"))
plt.savefig(img_path,dpi=300)
print('Plot printed:\n%s' % img_path)