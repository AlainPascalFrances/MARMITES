#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Alain Francés
#
# Created:     22-11-2012
# Copyright:   (c) alf 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import numpy as np
import os, sys
import matplotlib as mpl
if mpl.get_backend!='agg':
    mpl.use('agg')
import matplotlib.pyplot as plt

#####################
#input section
ws_fn = r'E:\00WIP\00MyPapers_ws\MMMFalgorithm\FrontiersWater' #E:\00code_ws\LaMata_new_PhD_artigo_2s2L\out_202110100955_2s2L_2013_SySsdim_DRN0p9_Kh1_full_1stSS_3runs'
# read input files
TXT_in_fnn = 'Eddy_MM_202110.csv'

# read XLS file
inputFile_TS_fn = os.path.join(ws_fn, TXT_in_fnn)
if os.path.exists(inputFile_TS_fn):
    data = np.loadtxt(inputFile_TS_fn, skiprows = 1, dtype = str, delimiter = ';')
else:
    sys.exit("\nFATAL ERROR!\nThe input file [" + inputFile_TS_fn + "] doesn't exist, verify name and path!")
date   = data[:,0]
date   = mpl.dates.datestr2num(date)
ET_MM  = data[:,2].astype(float)
RF     = data[:,3].astype(float)
ET_ECT  = data[:,4].astype(float)
Ro_MM  = data[:,5].astype(float)
Ro_obs = data[:,6].astype(float)

# end of reading data from TXT file

#####################################
def compRMSE (sim, obs) :
    def sqre_diff (v, w) :
        return (v - w) ** 2
    s = len(sim)
    v = sum(map(sqre_diff, sim, obs))
    return np.sqrt(v / s)

#####################################

def compE (sim, obs, hnoflo) :
    def sqre_diff (v, w) :
        return (v - w) ** 2
    #s = len(sim)
    if sum(sqre_diff(obs, compAVGE(obs)))>0.0:
        v = 1 - sum(map(sqre_diff, sim, obs))/sum(sqre_diff(obs, compAVGE(obs)))
    else:
        v = hnoflo
    return v

#####################################

def compAVGE(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

#####################################

def compR(x, y, hnoflo):
    # pearson correlation coefficient r
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = compAVGE(x)
    avg_y = compAVGE(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    if np.sqrt(xdiff2 * ydiff2)> 0.0 :
        return diffprod / np.sqrt(xdiff2 * ydiff2)
    else:
        return hnoflo

#####################################

def compCalibCrit(sim, obs, hnoflo): 
    rmse = None
    rsr  = None
    nse  = None
    r    = None
    aa = np.array([sim,obs])
    aa = np.transpose(aa)
    if hnoflo > 0:
        bb = aa[~(aa > hnoflo - 1000.0).any(1)]
    else:
        bb = aa[~(aa < hnoflo + 1000.0).any(1)]
    if bb[:,0].any():
        rmse = (compRMSE(bb[:,0], bb[:,1]))
        if np.std(bb[:,1]) > 0:
            rsr = (rmse/(np.std(bb[:,1])))
        nse = (compE(bb[:,0], bb[:,1], hnoflo))
        r = (compR(bb[:,0], bb[:,1], hnoflo))
    return rmse, rsr, nse, r

#####################################
hnoflo = 99999.99
ECT = np.where(~np.isnan(ET_ECT),ET_ECT,hnoflo)
MM = np.where(~np.isnan(ET_MM),ET_MM,hnoflo)
print('ET_MM vs ET_ECT\nrmse, rsr, nse, r')
print('%.2f, %.2f, %.2f, %.2f' % compCalibCrit(ECT, MM, hnoflo))

dateFmt=mpl.dates.DateFormatter('%Y-%b-%d')
dateminorFmt=mpl.dates.DateFormatter('%b')
lblspc = 0.05
mkscale = 0.5
bdpd = 0.1
hdltxtpd = 0.05
colspc = 0.1

plt_suptitle = 'ET ECT vs ET MM + RF' 
iniMonthHydroYear = 10
date_ini = mpl.dates.datestr2num('2008-10-01')-15
date_end = mpl.dates.datestr2num('2010-09-30')+15

fig = plt.figure(num=None, figsize=(8.27, 11.7), dpi = 60)    #(8.5,15), dpi=30)
fig.suptitle(plt_suptitle)

# # Ro
# ax0=fig.add_subplot(8,1,1) #8,1,1
# plt.setp(ax0.get_xticklabels(), visible = False)
# plt.setp(ax0.get_yticklabels(), fontsize=8)
# plt.plot_date(date,Ro_MM,fmt = '-', c='blue', linewidth=1.0, label = 'MM')
# plt.ylabel('$Ro$ (mm)', fontsize=10)
# plt.plot_date(date,Ro_obs, marker='o', ls = 'None', color = 'lightblue', markeredgecolor = 'green', markerfacecolor = 'lightgreen', markersize = 2, label = 'flume') # ls='--', color = 'blue'
# ymax = np.ma.max(np.ma.masked_invalid(Ro_obs))
# plt.ylim(0, 0.6)
# plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 5)
# leg = plt.gca().get_legend()
# ltext  = leg.get_texts()
# plt.setp(ltext, fontsize=8)
# ax0.grid(b=True, which='major', axis = 'both')
# ax0.xaxis.grid(b=True, which='minor', color='0.65')
# ax0.xaxis.set_major_formatter(dateFmt)
# ax0.xaxis.set_major_locator(mpl.dates.YearLocator(1, month = iniMonthHydroYear, day = 1))
# bymonth = []
# month_tmp = 3
# while len(bymonth)<3:
#    if (iniMonthHydroYear+month_tmp) <13:
#        bymonth.append(iniMonthHydroYear+month_tmp)
#    else:
#        bymonth.append(iniMonthHydroYear+month_tmp - 12)
#    month_tmp += 3
# del month_tmp
# ax0.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth = bymonth))
# ax0.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
# plt.setp(ax0.get_xticklabels(minor=True), visible=False)

# ET
ax1=fig.add_subplot(4,1,1) #8,1,2, sharex=ax0)
plt.setp(ax1.get_xticklabels(), fontsize=8)
plt.setp(ax1.get_yticklabels(), fontsize=8)
plt1a = plt.plot_date(date,ET_ECT, fmt = '-', c='red', linewidth=1.0, label='ECT')
plt1b = plt.plot_date(date,ET_MM, fmt = '-', c='orange', linewidth=1.0, label='MM')
plt.xlabel('Date', fontsize=10)
plt.ylabel('$ET$ (mm)', fontsize=10)
ax1.set_ylim(0.0,6.0)

ax2 = ax1.twinx()
plt2a = ax2.bar(date,RF,color='darkblue', linewidth=0, align = 'center', label='RF')
for tl in ax2.get_yticklabels():
    tl.set_color('darkblue')
    tl.set_fontsize(8)
plt.ylabel('$RF$ (mm)', fontsize=10, color = 'darkblue')
ax2.set_ylim(0,60)
plt.gca().invert_yaxis()
plts =  plt1a + plt1b + [plt2a]
labs = [l.get_label() for l in plts]
plt.legend(plts, labs, loc=6, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 3, columnspacing = colspc, numpoints = 5)
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=8)

ax1.grid(b=True, which='major', axis = 'both')
ax1.xaxis.grid(b=True, which='minor', color='0.65')
ax1.xaxis.set_major_formatter(dateFmt)
ax1.xaxis.set_major_locator(mpl.dates.YearLocator(1, month = iniMonthHydroYear, day = 1))
bymonth = []
month_tmp = 3
while len(bymonth)<3:
    if (iniMonthHydroYear+month_tmp) <13:
        bymonth.append(iniMonthHydroYear+month_tmp)
    else:
        bymonth.append(iniMonthHydroYear+month_tmp - 12)
    month_tmp += 3
del month_tmp
ax1.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth = bymonth))
ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
plt.setp(ax1.get_xticklabels(minor=True), visible=False)


ax1.grid(b=True, which='major', axis = 'both')
ax1.xaxis.grid(b=True, which='minor', color='0.65')
plt.setp(ax1.get_xticklabels(minor=True), visible=True)
plt.xlabel('Date', fontsize=10)
labels=ax1.get_xticklabels()
plt.setp(labels, 'rotation', 90)
ax1.xaxis.set_minor_formatter(dateminorFmt)
labels=ax1.get_xminorticklabels()
plt.setp(labels, fontsize=8)
plt.setp(labels, 'rotation', 90)
del labels

ax1.set_xlim(date_ini,date_end)

#plt.show()
plt_export_fn = '%s\%s.png' % (ws_fn, TXT_in_fnn.split('.')[0])
plt.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
plt.savefig(plt_export_fn,dpi=150)
print('Plot printed:\n%s' % plt_export_fn)