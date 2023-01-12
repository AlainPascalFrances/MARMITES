# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:25:53 2014

@author: apf
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
#import CreateColors
import itertools, datetime, adjustText
mpl.rcParams['mathtext.fontset'] = 'stix'
#matplotlib.rcParams['font.family'] = 'STIXGeneral'
mpl.pyplot.legend(r'ABC123 vs $\mathrm{ABC123}^{123}$')

# INPUT
phi = 0.45
theta_fc = 0.35
theta_wp = 0.05
ThickSoil = 1000.0
PT = 4.0
plt_numb_lst = ['b) - Transpiration simulation']
pos = [2,4,6,8]
daynum= 350
#veg param
f_lst = [0.102]#[0.10, 0.217336864608619, 0.30]
s_lst = [22.15]#[20, 28.26674282, 35]
kTs_max_lst = [0.955]#[0.9, 0.712676596, 0.35]
kTs_min_lst = [0.136]#[0.6, 0.19772923, 0.05]
# kTs obs
kTs_source = ['$Reyes\ (2015),\ \it{Q.i.}$', '$Reyes\ (2015),\ \it{Q.p.}$', '$David\ \it{et\ al.},\ (2013),\ \it{Q.s.}$', '$Pinto\ \it{et\ al.},\ (2014),\ \it{Q.s.}$']
kTs_color = ['red','orange','green','limegreen']
theta_obs = [[0.580,0.108],
             [0.437,0.157],
             [0.119047619,0.19047619,0.071428571,0.428571429,0.547619048,0.666666667,0.666666667,0.428571429,0.726190476,0.19047619,0.071428571,0],
             [0.202020202,0.141414141,0.393939394,0.606060606,0.555555556,0.727272727,0.424242424,0.535353535,0.282828283,0.141414141,0.080808081,0.04040404]]
kTs_obs = [[0.040,0.573],
          [0.7,0.205],
          [0.709367946,0.854499043,0.73393802,1,0.940585009,0.955264969,1,0.954905545,0.971951918,0.622202057,0.234186979,0.213772708],
          [0.829109881,1,0.764239532,1,1,1,1,1,1,0.786471594,0.013845341,0.016042992]]
month_lst = [['06', '09'],
         ['06', '09'],
         ['09','10','11','12','01','02','03','04','05','06','07','08'],
         ['09','10','11','12','01','02','03','04','05','06','07','08']]

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
n_y = 0.7500

incr = 0.001

s_inv_lst = []
for s in s_lst:
    s_inv_lst.append(1/s)

fig = plt.figure(figsize=(11.7, 8.27))
ax1=fig.add_subplot(1,2,1)
plt.setp(ax1.get_xticklabels(), fontsize=8)
plt.setp(ax1.get_yticklabels(), fontsize=8)
title = plt.title(r'a) - $k_{T_{soil}}$ function of $\theta_{norm}$', fontsize = 12)
theta = np.arange(theta_wp,phi+incr, incr)
thetanorm = (theta - theta_wp) / (phi - theta_wp)
linesTg = itertools.cycle(['-','--','-.',':'])
plts = []
x = []
y = []
m = []
texts = []
for l, (s, s_inv, f, kTs_max, kTs_min, plt_numb) in enumerate(zip(s_inv_lst, s_lst, f_lst, kTs_max_lst, kTs_min_lst, plt_numb_lst)):
    kTs = kTs_max - (kTs_max-kTs_min)/(1 + np.exp((thetanorm - f)/s))
    ax1.plot(thetanorm, kTs, linesTg.__next__(), color = cTg)
    #ax1.annotate('$1/s = %.1f,\ f = %.2f$\n$k_{T_{max}} = %.1f,\ k_{T_{min}} = %.1f$\n$(T_g\ parameters\ of\ figure\ %s$' % (s_inv, f, kTs_max, kTs_min, plt_numb), (f, kTs_max - (kTs_max-kTs_min)/(1 + np.exp((f - f)/s))), horizontalalignment='left', verticalalignment='bottom', bbox=dict(facecolor='white'), fontsize = 10, color = 'black')
for l, (lbl,c,month) in enumerate(zip(kTs_source,kTs_color,month_lst)):
        plts += ax1.plot(theta_obs[l],kTs_obs[l], 'o', label = lbl, color = c, markersize = 5)
        for e, (x_, y_, m_) in enumerate(zip(theta_obs[l], kTs_obs[l], month)):
            x.append(x_)
            y.append(y_)
            m.append(m_)
for e, (x_, y_, m_) in enumerate(zip(x, y, m)):
    texts.append(plt.text(x_, y_, m_, fontsize = 8, color = 'black'))
adjustText.adjust_text(texts, add_objects = [title]) #, expand_objects = (1.25, 1.5), force_objects  = (0.15,0.3))
plt.grid(True)
plt.xlim((0.0,1.0))
plt.ylim((0.0,1.0))
plt.xlabel(r'$\theta_{norm}$')
plt.ylabel(r'$k_{T_{soil}}$')  #, color = cTg)
#"plt.setp(ax1.get_yticklabels(), color = cTg)
ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
#ax1.text(theta_wp-0.015, kTs_min, '${k_{T_s}}_{min}$', horizontalalignment='center', verticalalignment='center', fontsize = 10)
#ax1.text(phi+0.015, kTs_max, '${k_{T_s}}_{max}$', horizontalalignment='center', verticalalignment='center', fontsize = 10)
#ax1.text(phi, -0.0375, r'$\phi$', horizontalalignment='center', verticalalignment='center', fontsize = 10)
#ax1.text(theta_fc, -0.0375, r'$\theta_{fc}$', horizontalalignment='center', verticalalignment='center', fontsize = 10)
#ax1.text(theta_wp, -0.0375, r'$\theta_{wp}$', horizontalalignment='center', verticalalignment='center', fontsize = 10)
#labs = [l.get_label() for l in ax1.get_legend_handles_labels()]
plt.legend(plts, ax1.get_legend_handles_labels()[1], loc=7, labelspacing=lblspc, markerscale=1, handletextpad=handletextpad,
           borderaxespad=borderaxespad, borderpad=borderpad)
leg = plt.gca().get_legend()
ltext = leg.get_texts()
plt.setp(ltext, fontsize=10)

linesTg = itertools.cycle(['-','--','-.',':'])
for j, (s, f, kTs_max, kTs_min, plt_numb) in enumerate(zip(s_inv_lst, f_lst, kTs_max_lst, kTs_min_lst, plt_numb_lst)):
    days = np.arange(1,daynum + 1,1)
    Tsoil = np.zeros(daynum, float)
    Tg    = np.zeros(daynum, float)
    T     = np.zeros(daynum, float)
    PTg   = np.zeros(daynum, float)
    theta = np.zeros(daynum, float)
    PTarray = np.ones([daynum], float)*PT
    for d in range(len(days)):
        if d == 0:
            theta_tmp = phi
        else:
            theta_tmp = theta[d-1]
        thetanorm = (theta_tmp - theta_wp) / (phi - theta_wp)
        Tsoil[d] = PT * thetanorm
        theta[d] = (theta_tmp*ThickSoil - Tsoil[d])/ThickSoil
        PTg[d] = PT-Tsoil[d]
        kTs = kTs_max - (kTs_max-kTs_min)/(1 + np.exp((thetanorm - f)/s))
        Tg[d] = PTg[d]*kTs
        T[d] = Tsoil[d] + Tg[d]
    lines = itertools.cycle(['-','--','-.','-'])
    ax2=fig.add_subplot(1,2,2)
    plt.title('%s' % (plt_numb), fontsize = 12)
    plt.setp(ax2.get_xticklabels(), fontsize=8)
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    plt2_3 = ax2.plot(days, PTarray, '-', color = cPTa, label = r'$PT$')
    plt2_0 = ax2.plot(days, T, lines.__next__(), color = cT, label = r'$T$')
    plt2_1 = ax2.plot(days, Tsoil, lines.__next__(), color = cTs, label = r'$T_{soil}$')
    plt2_2 = ax2.plot(days, Tg, linesTg.__next__(), color = cTg, label = r'$T_g$')
    plt2_4 = ax2.plot(days, PTg, '-.', color = cPT, label = r'$PT_g$')
    plt.grid(True)
    plt.ylabel('mm')
    if j>2:
        plt.xlabel('day')
    plt.ylim((0,PT+0.1))
    ax2a = ax2.twinx()
    plt.setp(ax2a.get_yticklabels(), fontsize=8, color = cSM)
    plt2a = ax2a.plot(days, theta, lines.__next__(), color = cSM, label = r'$\theta$')
    plt.ylim((theta_wp-.05,phi+.01))
    plt.ylabel(r'$(\%)$', color = cSM)
    plt.xlim((0, daynum))
    ax2a.text(n_x, n_y, '$1/s = %.1f,\ f = %.2f$\n$k_{{T_s}_{min}} = %.2f,\ k_{{T_s}_{max}} = %.2f$' % (1/s, f, kTs_min, kTs_max), horizontalalignment='right', verticalalignment='center', transform=ax2a.transAxes, bbox=dict(facecolor='white'), fontsize = 10)
    if j == 0:
        plts =  plt2_3 + plt2_4 + plt2_0 + plt2_1 + plt2_2
    elif j==1:
        plts += (plt2_2)
    elif j==2:
        plts += (plt2_2+ plt2a)
labs = [l.get_label() for l in plts]
plt.legend(plts, labs, loc=7, labelspacing=lblspc, markerscale=mkscale, handletextpad = handletextpad, borderaxespad = borderaxespad, borderpad = borderpad, ncol = 1) #, columnspacing = 0.05)  #bbox_to_anchor=(-0.25, 2.5),
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=10)

plt.subplots_adjust(left=0.075, bottom=0.05, right=0.95, top=0.95, wspace=0.25, hspace=0.3)
#plt.show()
img_path = 'Tg_functionNEW1_%s'% (datetime.date.today().strftime("%Y%m%d"))
plt.savefig(img_path,dpi=300)
print('Plot printed:\n%s' % img_path)