# -*- coding: utf-8 -*-

__author__ = "Alain P. Francés <frances08512@itc.nl>"
__version__ = "0.3"
__date__ = "2012"

import matplotlib as mpl
import matplotlib.pyplot as plt
import subprocess as sp
import numpy as np
import os, itertools, shutil
from matplotlib.sankey import Sankey

#####################################

def plotTIMESERIES(cMF, i, j, flx, flxLbl, flxIndex_lst, Sm, Sr, plt_export_fn, plt_suptitle, plt_title, clr_lst, hmax, hmin, obs_name, l_obs, nsl, iniMonthHydroYear, date_ini, date_end):
    """
    Plot the time series of the fluxes observed at one point of the catchment
    Use Matplotlib
    _______________________________________________________________________________

    INPUTS
            STATE VARIABLES
                SP              Stress period
                P               Daily rainfall
                PT              Daily potential transpiration
                PE              Daily potential evaporation
                Pe              Daily Excess rainfall
                Esoil           Daily evaporation (bare soil)
                Tsoil           Daily transpiration
                S               Daily soil moisture
                Rp              Daily percolation
                POND            Daily ponding
                Ro              Daily runoff
                Rg               Daily gross recharge
                h               Daily water level
                hobs           Daily obsserved water level
    ______________________________________________________________________________
    ______________________________________________________________________________
    """

    #fmt_DH = mpl.dates.DateFormatter('%Y-%m-%d %H:%M')
#    cUTIL = MMutils.clsUTILITIES(fmt = fmt_DH)
#    date_ini, year_ini = cUTIL.compDATE_INI(cMF.inputDate[0], iniMonthHydroYear)
#    date_end, year_end = cUTIL.compDATE_END(cMF.inputDate[-1], iniMonthHydroYear)
    date_ini -= 15
    date_end += 15

    dateFmt=mpl.dates.DateFormatter('%Y-%b-%d')
    dateminorFmt=mpl.dates.DateFormatter('%b')
    lblspc = 0.05
    mkscale = 1.0
    bdpd = 0.1
    hdltxtpd = 0.05
    colspc = 0.1

    fig = plt.figure(num=None, figsize=(8.27, 11.7), dpi = 60)    #(8.5,15), dpi=30)
    fig.suptitle(plt_suptitle)
    fig.text(x = 0.5, y = 0.05, s = plt_title, horizontalalignment = 'center', verticalalignment = 'bottom', fontsize = 9)

    ax1=fig.add_subplot(8,1,1)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), fontsize=8)
    ax1.bar(cMF.inputDate,flx[flxIndex_lst['iRF']],color='darkblue', linewidth=0, align = 'center', label=flxLbl[flxIndex_lst['iRF']])
    ax1.bar(cMF.inputDate,flx[flxIndex_lst['iRFe']],color='deepskyblue', linewidth=0, align = 'center', label=flxLbl[flxIndex_lst['iRFe']])
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    ax1.grid(b=True, which='major', axis = 'both')
    ax1.xaxis.grid(b=True, which='minor', color='0.65')
    plt.ylabel('mm', fontsize=10)
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

    ax2=fig.add_subplot(8,1,2, sharex=ax1)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    try:
        Roobs_m = np.ma.masked_values(flx[flxIndex_lst['iRoobs']], cMF.hnoflo, atol = 0.09)
    #        dgwtobsmin = np.ma.min(dgwtobs_m)
        plt.plot_date(cMF.inputDate,Roobs_m, ls = 'None', color = 'lightblue', marker='o', markeredgecolor = 'lightblue', markerfacecolor = 'None', markersize = 2, label=flxLbl[flxIndex_lst['iRoobs']]) # ls='--', color = 'blue'
    except:
        pass
    if np.sum(np.abs(flx[flxIndex_lst['iRo']])) > 1E-7:
        plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iRo']],'r-', c='blue', linewidth=1.0, label=flxLbl[flxIndex_lst['iRo']])
    colors_nsl = itertools.cycle(clr_lst)
    if np.sum(np.abs(flx[flxIndex_lst['iInf']])) > 1E-7:
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst['iInf']],'--', c=colors_nsl.next(), linewidth=1.5, label=flxLbl[flxIndex_lst['iInf']])
    if np.sum(np.abs(flx[flxIndex_lst['iEsurf']])) > 1E-7:
        plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iEsurf']],'r-', c='deepskyblue', linewidth=0.75, label=flxLbl[flxIndex_lst['iEsurf']])
    if np.sum(np.abs(flx[flxIndex_lst['iSsurf']])) > 1E-7:
        plt.bar(cMF.inputDate, flx[flxIndex_lst['iSsurf']], color='lightblue', linewidth=0, align = 'center', label=flxLbl[flxIndex_lst['iSsurf']])
    if np.sum(np.abs(flx[flxIndex_lst['idSsurf']])) > 1E-7:
        plt.bar(cMF.inputDate, flx[flxIndex_lst['idSsurf']], color='darkblue', linewidth=0, align = 'center', label=flxLbl[flxIndex_lst['idSsurf']])
    plt.ylabel('mm', fontsize=10)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    ax2.grid(b=True, which='major', axis = 'both')
    ax2.xaxis.grid(b=True, which='minor', color='0.65')
    ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax2.get_xticklabels(minor=True), visible=False)

    colors_nsl = itertools.cycle(clr_lst)
    ax3=fig.add_subplot(8,1,3, sharex=ax1)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), fontsize=8)
    # PE
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iPE']],'-', color='lightblue', linewidth=3, label = flxLbl[flxIndex_lst['iPE']])
    if cMF.wel_yn == 1:
        E = flx[flxIndex_lst['iEsoil']] + flx[flxIndex_lst['iEg']]
        # Etot
        plt.plot_date(cMF.inputDate,E,'-', color='darkblue', linewidth=1.5, label = r'$E$')
    # Esoil
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iEsoil']],'--', color='brown', label = flxLbl[flxIndex_lst['iEsoil']])
    if cMF.wel_yn == 1:
        # Eg
        if np.absolute(sum(flx[flxIndex_lst['iEg']]))>1E-6:
            plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iEg']],'-', color='blue', label = flxLbl[flxIndex_lst['iEg']])
    for l in range(nsl):
        if np.absolute(sum(flx[flxIndex_lst['iEsoil_l%d'%(l+1)]])):
            ax3.plot_date(cMF.inputDate, flx[flxIndex_lst['iEsoil_l%d'%(l+1)]], '-', color=colors_nsl.next(), label = flxLbl[flxIndex_lst['iEsoil_l%d'%(l+1)]])
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 3, columnspacing = colspc)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    ax3.grid(b=True, which='major', axis = 'both')
    ax3.xaxis.grid(b=True, which='minor', color='0.65')
    plt.ylabel('mm', fontsize=10)
    ax3.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax3.get_xticklabels(minor=True), visible=False)

    colors_nsl = itertools.cycle(clr_lst)
    ax4=fig.add_subplot(8,1,4, sharex=ax1)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), fontsize=8)    
    # PT
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iPT']],'-', color='lightblue', linewidth=3, label = flxLbl[flxIndex_lst['iPT']])
    if cMF.wel_yn == 1:
        T = flx[flxIndex_lst['iTsoil']] + flx[flxIndex_lst['iTg']]
        # Ttot
        plt.plot_date(cMF.inputDate,T,'-', color='darkblue', linewidth=1.5, label = r'$T$')
    # Tsoil
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iTsoil']],'--', color='brown', label = flxLbl[flxIndex_lst['iTsoil']])
    if cMF.wel_yn == 1:
        # Tg
        if np.absolute(sum(flx[flxIndex_lst['iTg']]))>1E-6:
            plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iTg']],'-', color='blue', label = flxLbl[flxIndex_lst['iTg']])
    for l in range(nsl):
        if np.absolute(sum(flx[flxIndex_lst['iTsoil_l%d'%(l+1)]])):            
            ax4.plot_date(cMF.inputDate, flx[flxIndex_lst['iTsoil_l%d'%(l+1)]], '-', color=colors_nsl.next(), label = flxLbl[flxIndex_lst['iTsoil_l%d'%(l+1)]])
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 3, columnspacing = colspc)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax4.grid(b=True, which='major', axis = 'both')
    ax4.xaxis.grid(b=True, which='minor', color='0.65')
    plt.ylabel('mm', fontsize=10)
    ax4.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax4.get_xticklabels(minor=True), visible=False)

    colors_nsl = itertools.cycle(clr_lst)
    ax6=fig.add_subplot(8,1,5, sharex=ax1)
    plt.setp(ax6.get_xticklabels(), visible=False)
    plt.setp(ax6.get_yticklabels(), fontsize=8)
    if np.sum(np.abs(flx[flxIndex_lst['iEXFg']])) > 1E-7:
        plt.bar(cMF.inputDate,flx[flxIndex_lst['iEXFg']], color='lightblue', linewidth=0, align = 'center', label=flxLbl[flxIndex_lst['iEXFg']])
    if cMF.wel_yn == 1:
        if np.sum(np.abs(flx[flxIndex_lst['iETg']])) > 1E-7:
            plt.plot_date(cMF.inputDate, flx[flxIndex_lst['iETg']],'-', c='blue', linewidth=1.5, label=flxLbl[flxIndex_lst['iETg']])        
    for l in range(nsl):
        if np.sum(np.abs(flx[flxIndex_lst['iRsoil_l%d'%(l+1)]])) > 1E-7:
            ax6.plot_date(cMF.inputDate, -1.0*flx[flxIndex_lst['iRsoil_l%d'%(l+1)]], '-', color=colors_nsl.next(), label=flxLbl[flxIndex_lst['iRsoil_l%d'%(l+1)]])
    if np.sum(np.abs(flx[flxIndex_lst['iRg']])) > 1E-7:
        plt.plot_date(cMF.inputDate, -1.0*flx[flxIndex_lst['iRg']],'-', c='darkblue', linewidth=1.5, label=flxLbl[flxIndex_lst['iRg']])
    colors_nsl = itertools.cycle(clr_lst)
    for l in range(nsl):
        if np.sum(np.abs(flx[flxIndex_lst['iExf_l%d'%(l+1)]])) > 1E-7:
            ax6.plot_date(cMF.inputDate, flx[flxIndex_lst['iExf_l%d'%(l+1)]], '--', color=colors_nsl.next(), label=flxLbl[flxIndex_lst['iExf_l%d'%(l+1)]])        
    plt.legend(labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 3, columnspacing = colspc, loc = 0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax6.grid(b=True, which='major', axis = 'both')
    ax6.xaxis.grid(b=True, which='minor', color='0.65')
    plt.ylabel('mm', fontsize=10)
    ax6.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax6.get_xticklabels(minor=True), visible=False)

    colors_nsl = itertools.cycle(clr_lst)
    ax5=fig.add_subplot(8,1,6, sharex=ax1)
    plt.setp(ax5.get_xticklabels(), visible=False)
    plt.setp(ax5.get_yticklabels(), fontsize=8)
    for l in range(nsl):
        try:
            if flx[flxIndex_lst['iSobs_l%d'%(l+1)]] != []:
                ax5.plot_date(cMF.inputDate, flx[flxIndex_lst['iSobs_l%d'%(l+1)]], ls = 'None', color = 'gray', marker='o', markersize=2, markeredgecolor = colors_nsl.next(), markerfacecolor = 'None', label=flxLbl[flxIndex_lst['iSobs_l%d'%(l+1)]]) #'--', color = color,
        except:
            pass
    sim_tmp = []
    colors_nsl = itertools.cycle(clr_lst)
    for l in range(nsl):
        y = flx[flxIndex_lst['iSsoil_pc_l%d'%(l+1)]]
        sim_tmp.append(y)
        y = np.ma.masked_where(y < 0.0, y)
        ax5.plot_date(cMF.inputDate, y, '-', color = colors_nsl.next(), label = flxLbl[flxIndex_lst['iSsoil_pc_l%d'%(l+1)]])
        del y
    # y axis
    ybuffer=0.1*(max(Sm)-min(Sr))
    plt.ylim(min(Sr) - ybuffer,max(Sm) + ybuffer)
    plt.ylabel('%', fontsize=10)
    # legend
    #lbl_Spcobs = lbl_Sobs + lbl_S
    #plt.legend(lbl_Spcobs, loc=0, labelspacing=lblspc, markerscale=mkscale)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax5.grid(b=True, which='major', axis = 'both')
    ax5.xaxis.grid(b=True, which='minor', color='0.65')
    ax5.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax5.get_xticklabels(minor=True), visible=False)

    ax7=fig.add_subplot(8,1,7, sharex=ax1)
    plt.setp(ax7.get_xticklabels(), fontsize=8)
    plt.setp(ax7.get_yticklabels(), fontsize=8)
    try:
        plt.plot_date(cMF.inputDate,flx[flxIndex_lst['idobs']], ls = 'None', color = 'LightBlue', marker='o', markeredgecolor = 'LightBlue', markerfacecolor = 'None', markersize = 2, label = flxLbl[flxIndex_lst['idobs']]) # ls='--', color = 'blue'
    except:
        pass
    lines = itertools.cycle(['-','--','-.',':','.',',','o','v','^','<','>','s','p','*','h','H','+','x','D','d','|','_'])
    dgwtMFmax = []
    for L in range(cMF.nlay):
        dgwtMF = flx[flxIndex_lst['id_L%d'%(L+1)]]
        dgwtMFmax.append(np.max(dgwtMF))
        plt.plot_date(cMF.inputDate, dgwtMF, lines.next(), color = 'b', label = flxLbl[flxIndex_lst['id_L%d'%(L+1)]])
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['idcorr']],'--', c='g', label = flxLbl[flxIndex_lst['idcorr']])
    # y axis
    plt.ylabel('m', fontsize=10)
    ax7.grid(b=True, which='major', axis = 'both')
    ax7.xaxis.grid(b=True, which='minor', color='0.65')
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
    leg = plt.gca().get_legend()
    dgwtMFmax = np.max(dgwtMFmax)
    if dgwtMFmax > 0:
        ymax = dgwtMFmax * 1.05
    else:
        ymax = 0.25
#    if obs_leg <> None:
#        if dgwtobsmin < 0.0:
#            plt.ylim(ymax = ymax, ymin = dgwtobsmin * 1.05)
#        else:
#            plt.ylim(ymax = ymax)
#    else:
#        plt.ylim(ymax = ymax)
    plt.ylim(ymax = ymax)
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax7.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax7.get_xticklabels(minor=True), visible=True)
    plt.xlabel('Date', fontsize=10)
    labels=ax7.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    ax7.xaxis.set_minor_formatter(dateminorFmt)
    labels=ax7.get_xminorticklabels()
    plt.setp(labels, fontsize=8)
    plt.setp(labels, 'rotation', 90)
    del labels

    ax1.set_xlim(date_ini,date_end)

    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig(plt_export_fn,dpi=150)


#--------------FIGURE PART2----------

    fig = plt.figure(num=None, figsize=(8.27, 11.7), dpi = 60)    #(8.5,15), dpi=30)
    fig.suptitle(plt_suptitle + ' - part 2')
    fig.text(x = 0.5, y = 0.05, s = plt_title, horizontalalignment = 'center', verticalalignment = 'bottom', fontsize = 9)

    ax1=fig.add_subplot(8,1,4)
    plt.setp(ax1.get_xticklabels(), visible = False)
    plt.setp(ax1.get_yticklabels(), fontsize=8)
    try:
        plt.plot_date(cMF.inputDate,flx[flxIndex_lst['ihobs']], ls = 'None', color = 'LightBlue', marker='o', markeredgecolor = 'LightBlue', markerfacecolor = 'None', markersize = 2, label = flxLbl[flxIndex_lst['ihobs']]) # ls='--', color = 'blue'
    except:
        pass
    lines = itertools.cycle(['-','--','-.',':','.',',','o','v','^','<','>','s','p','*','h','H','+','x','D','d','|','_'])
    for L in range(cMF.nlay):
        plt.plot_date(cMF.inputDate,flx[flxIndex_lst['ih_L%d'%(L+1)]],lines.next(), color = 'b',label = flxLbl[flxIndex_lst['ih_L%d'%(L+1)]])
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['ihcorr']],'--', color = 'g', label = flxLbl[flxIndex_lst['ihcorr']])
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['ih_SF']],'-', color = 'r', label = flxLbl[flxIndex_lst['ih_SF']])
    ybuffer=0.1*(hmax-hmin)
    if ybuffer == 0.0:
        ybuffer = 1.0
    plt.ylim((hmin - ybuffer, hmax + ybuffer))
    plt.ylabel('m', fontsize=10)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax1.grid(b=True, which='major', axis = 'both')
    ax1.xaxis.grid(b=True, which='minor', color='0.65')
    ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))

    colors_nsl = itertools.cycle(clr_lst)
    ax8b=fig.add_subplot(8,1,1, sharex=ax1)
    plt.setp(ax8b.get_xticklabels(), visible=False)
    plt.setp(ax8b.get_yticklabels(), fontsize=8)
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iMBsurf']],'-', c='lightblue', label = flxLbl[flxIndex_lst['iMBsurf']])
    MBmin = [min(flx[flxIndex_lst['iMBsurf']])]
    MBmax = [max(flx[flxIndex_lst['iMBsurf']])]
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iMB']],'-', c='r', label = flxLbl[flxIndex_lst['iMB']])
    MBmin = [min(flx[flxIndex_lst['iMB']])]
    MBmax = [max(flx[flxIndex_lst['iMB']])]
    for l in range(nsl):
        ax8b.plot_date(cMF.inputDate, flx[flxIndex_lst['iMB_l%d'%(l+1)]], '-', color=colors_nsl.next(), label=flxLbl[flxIndex_lst['iMB_l%d'%(l+1)]])
        MBmin.append(min(flx[flxIndex_lst['iMB_l%d'%(l+1)]]))
        MBmax.append(max(flx[flxIndex_lst['iMB_l%d'%(l+1)]]))
    # y axis
    plt.ylabel('mm', fontsize=10)
    ax8b.grid(b=True, which='major', axis = 'both')
    ax8b.xaxis.grid(b=True, which='minor', color='0.65')
    MBmax = max(MBmax)
    MBmin = min(MBmin)
    if np.abs(MBmax - MBmin) < 1E-7:
        plt.ylim(-0.1,0.1)
    else:
        minfact = 0.95
        maxfact = 1.05
        if MBmin > 0:
            minfact = 1.05
        if MBmax < 0:
            maxfact = 0.95
        plt.ylim(MBmin*minfact,MBmax*maxfact)
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax8b.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax8b.get_xticklabels(minor=True), visible=False)

    colors_nsl = itertools.cycle(clr_lst)
    ax9a=fig.add_subplot(16,1,5, sharex=ax1)
    plt.setp(ax9a.get_xticklabels(), visible=False)
    plt.setp(ax9a.get_yticklabels(), fontsize=8)
    for l in range(nsl):
        ax9a.plot_date(cMF.inputDate, flx[flxIndex_lst['iSAT_l%d'%(l+1)]], '-', color = colors_nsl.next(), label = flxLbl[flxIndex_lst['iSAT_l%d'%(l+1)]])
    # y axis
    plt.ylim(-0.1,1.1)
    plt.ylabel('SAT', fontsize=10)
    ax9a.yaxis.set_ticks(np.arange(0,1.25,1))
    ax9a.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%1d'))
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax9a.grid(b=True, which='major', axis = 'both')
    ax9a.xaxis.grid(b=True, which='minor', color='0.65')
    plt.setp(ax9a.get_xticklabels(minor=True), visible=False)

    ax10a=fig.add_subplot(16,1,6, sharex=ax1)
    plt.setp(ax10a.get_xticklabels(), visible=False)
    plt.setp(ax10a.get_yticklabels(), fontsize=8)
    uzthick = flx[flxIndex_lst['iuzthick']]
    plt.plot_date(cMF.inputDate,-uzthick,'-', c='brown', label = flxLbl[flxIndex_lst['iuzthick']])
    # y axis
    plt.ylabel('m', fontsize=10)
    ax10a.grid(b=True, which='major', axis = 'both')
    ax10a.xaxis.grid(b=True, which='minor', color='0.65')
    minfact = 0.95
    maxfact = 1.05
    if np.min(uzthick) < 0:
        minfact = 1.05
    if np.max(uzthick) < 0:
        maxfact = 0.95
    if np.min(uzthick) == np.max(uzthick) == 0.0:
        plt.ylim(-0.5, 0.5)
    else:
        plt.ylim(np.min(-uzthick)*minfact, np.max(-uzthick)*maxfact)
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax10a.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax10a.get_xticklabels(minor=True), visible=False)

    colors_nsl = itertools.cycle(clr_lst)
    ax10b=fig.add_subplot(8,1,6, sharex=ax1)
    plt.setp(ax10b.get_xticklabels(), visible = False)
    plt.setp(ax10b.get_yticklabels(), fontsize=8)
    for l in range(nsl):
        y = flx[flxIndex_lst['iSsoil_l%d'%(l+1)]]
        y = np.ma.masked_where( y < 0.0, y)
        ax10b.plot_date(cMF.inputDate, y, '-', color=colors_nsl.next(), label=flxLbl[flxIndex_lst['iSsoil_l%d'%(l+1)]])
    # y axis
    plt.ylim(0,np.max(flx[flxIndex_lst['iSsoil_l%d'%(l+1)]])*1.05)
    plt.ylabel('mm', fontsize=10)
    ax10b.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax10b.grid(b=True, which='major', axis = 'both')
    ax10b.xaxis.grid(b=True, which='minor', color='0.65')

    ax20=fig.add_subplot(8,1,2, sharex=ax1)
    plt.setp(ax20.get_xticklabels(), visible = False)
    plt.setp(ax20.get_yticklabels(), fontsize=8)    
    ymax = None
    try:
    #        dgwtobsmin = np.ma.min(dgwtobs_m)
        plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iRoobs']], ls = 'None', color = 'lightblue', marker='o', markeredgecolor = 'lightblue', markerfacecolor = 'None', markersize = 2, label = flxLbl[flxIndex_lst['iRoobs']]) # ls='--', color = 'blue'
        ymax = np.ma.max(flx[flxIndex_lst['iRoobs']])
    except:
        pass
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iRo']],'r-', c='blue', linewidth=1.0, label = flxLbl[flxIndex_lst['iRo']])
    plt.ylabel('mm', fontsize=10)
    if ymax != None:
        plt.ylim(ymax = ymax)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    ax20.grid(b=True, which='major', axis = 'both')
    
    colors_nsl = itertools.cycle(clr_lst)
    ax5=fig.add_subplot(8,1,7, sharex=ax1)
    plt.setp(ax5.get_xticklabels(), fontsize=8)
    plt.setp(ax5.get_yticklabels(), fontsize=8)
    for l in range(nsl):
        try:
            y = flx[flxIndex_lst['iSobs_l%d'%(l+1)]]
            ax5.plot_date(cMF.inputDate, y, ls = 'None', color = 'gray', marker='o', markersize=2, markeredgecolor = colors_nsl.next(), markerfacecolor = 'None', label=flxLbl[flxIndex_lst['iSobs_l%d'%(l+1)]]) #'--', color = color,
        except:
            pass
    colors_nsl = itertools.cycle(clr_lst)
    for l in range(nsl):
        y = flx[flxIndex_lst['iSsoil_pc_l%d'%(l+1)]]
        y = np.ma.masked_where(y < 0.0, y)
        ax5.plot_date(cMF.inputDate, y, '-', color = colors_nsl.next(), label = flxLbl[flxIndex_lst['iSsoil_pc_l%d'%(l+1)]])
    # y axis
    ybuffer=0.1*(max(Sm)-min(Sr))
    plt.ylim(min(Sr) - ybuffer,max(Sm) + ybuffer)
    plt.ylabel('%', fontsize=10)
    # legend
    #lbl_Spcobs = lbl_Sobs + lbl_S
    #plt.legend(lbl_Spcobs, loc=0, labelspacing=lblspc, markerscale=mkscale)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax5.grid(b=True, which='major', axis = 'both')
    ax5.xaxis.grid(b=True, which='minor', color='0.65')
    ax5.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax7.get_xticklabels(minor=True), visible=True)
    plt.xlabel('Date', fontsize=10)
    labels=ax7.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    ax7.xaxis.set_minor_formatter(dateminorFmt)
    labels=ax7.get_xminorticklabels()
    plt.setp(labels, fontsize=8)
    plt.setp(labels, 'rotation', 90)
    del labels

    ax7=fig.add_subplot(8,1,5, sharex=ax1)
    plt.setp(ax7.get_xticklabels(), visible = False)
    plt.setp(ax7.get_yticklabels(), fontsize=8)
    try:
        if flx[flxIndex_lst['idobs']]!= []:
            plt.plot_date(cMF.inputDate,flx[flxIndex_lst['idobs']], ls = 'None', color = 'LightBlue', marker='o', markeredgecolor = 'LightBlue', markerfacecolor = 'None', markersize = 2, label = flxLbl[flxIndex_lst['idobs']]) # ls='--', color = 'blue'
    except:
        pass
    lines = itertools.cycle(['-','--','-.',':','.',',','o','v','^','<','>','s','p','*','h','H','+','x','D','d','|','_'])
    dgwtMFmax = []
    for L in range(cMF.nlay):
        dgwtMFmax.append(np.max(flx[flxIndex_lst['id_L%d'%(L+1)]]))   
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst['id_L%d'%(L+1)]], lines.next(), color = 'b', label = flxLbl[flxIndex_lst['id_L%d'%(L+1)]])
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['idcorr']],'--', c='g', label = flxLbl[flxIndex_lst['idcorr']])
    # y axis
    plt.ylabel('m', fontsize=10)
    ax7.grid(b=True, which='major', axis = 'both')
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
    leg = plt.gca().get_legend()
    dgwtMFmax = np.max(dgwtMFmax)
    if dgwtMFmax > 0:
        ymax = dgwtMFmax * 1.05
    else:
        ymax = 0.25
#    if obs_leg <> None:
#        if dgwtobsmin < 0.0:
#            plt.ylim(ymax = ymax, ymin = dgwtobsmin * 1.05)
#        else:
#            plt.ylim(ymax = ymax)
#    else:
#        plt.ylim(ymax = ymax)
    plt.ylim(ymax = ymax)
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax7.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax7.get_xticklabels(minor=True), visible=False)
    ax7.xaxis.grid(b=True, which='minor', color='0.65')
    ax7.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))

    ax1.set_xlim(np.min(cMF.inputDate)-15.0, np.max(cMF.inputDate)+15)
    
    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    txt = plt_export_fn.split('.')
    plt_export_fn = txt[0] + '_part2.' + txt[1]
    plt.savefig(plt_export_fn,dpi=150)
#    plt.show()
    plt.clf()
    plt.close('all')

##################

def plotTIMESERIES_flxGW(cMF, flx, flxLbl, flxIndex_lst, plt_export_fn, plt_title, iniMonthHydroYear, date_ini, date_end):
    """
    Plot the time series of the fluxes observed from the whole catchment
    Use Matplotlib
    """
    dateFmt=mpl.dates.DateFormatter('%Y-%b-%d')
    dateminorFmt=mpl.dates.DateFormatter('%b')
    lblspc = 0.05
    mkscale = 1.0
    bdpd = 0.1
    hdltxtpd = 0.05
    colspc = 0.1

    date_ini -= 15
    date_end += 15

    fig = plt.figure(num=None, figsize=(8.27, 11.7), dpi = 60)    #(8.5,15), dpi=30)
    fig.suptitle(plt_title + ' - part 3')
    
    ax0=fig.add_subplot(8,1,1)
    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.setp(ax0.get_yticklabels(), fontsize=8)
    # RF
    ax0.bar(cMF.inputDate,flx[flxIndex_lst['iRF']],color='darkblue', linewidth=0, align = 'center', label = flxLbl[flxIndex_lst['iRF']])
    # RFe
    ax0.bar(cMF.inputDate,flx[flxIndex_lst['iRFe']],color='deepskyblue', linewidth=0, align = 'center', label = flxLbl[flxIndex_lst['iRFe']])
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    ax0.grid(b=True, which='major', axis = 'both')
    ax0.xaxis.grid(b=True, which='minor', color='0.65')
    plt.ylabel('mm', fontsize=10)
    ax0.xaxis.set_major_formatter(dateFmt)
    ax0.xaxis.set_major_locator(mpl.dates.YearLocator(1, month = iniMonthHydroYear, day = 1))
    bymonth = []
    month_tmp = 3
    while len(bymonth)<3:
        if (iniMonthHydroYear+month_tmp) <13:
            bymonth.append(iniMonthHydroYear+month_tmp)
        else:
            bymonth.append(iniMonthHydroYear+month_tmp - 12)
        month_tmp += 3
    del month_tmp
    ax0.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth = bymonth))
    plt.setp(ax0.get_xticklabels(minor=True), visible=False)
    ax0.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))   
    
    # plot Ro
    ax2=fig.add_subplot(8,1,2, sharex=ax0)
    plt.setp(ax2.get_xticklabels(), fontsize=8)
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    # Ro
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iRo']],'r-', c='blue', linewidth=1.0, label = flxLbl[flxIndex_lst['iRo']])
    plt.ylabel('mm', fontsize=10)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    ax2.xaxis.grid(b=True, which='minor', color='0.65')
    ax2.grid(b=True, which='major', axis = 'both')
    ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax2.get_xticklabels(minor=True), visible=True)
    plt.xlabel('Date', fontsize=10)
    labels=ax2.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    ax2.xaxis.set_minor_formatter(dateminorFmt)
    labels=ax2.get_xminorticklabels()
    plt.setp(labels, fontsize=8)
    plt.setp(labels, 'rotation', 90)
    del labels    
    
    # plot gwd
    lines = itertools.cycle(['-','--','-.',':','.',',','o','v','^','<','>','s','p','*','h','H','+','x','D','d','|','_'])
    ax1=fig.add_subplot(8,1,4, sharex = ax0)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), fontsize=8)
    for L in range(cMF.nlay):
        plt.plot_date(cMF.inputDate,flx[flxIndex_lst['id_L%d'%(L+1)]],lines.next(), color = 'b', label = flxLbl[flxIndex_lst['id_L%d'%(L+1)]])
    plt.ylabel('m', fontsize=10)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax1.grid(b=True, which='major', axis = 'both')
    ax1.xaxis.grid(b=True, which='minor', color='0.65')
    ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax1.get_xticklabels(minor=True), visible=False)

    # plot GW fluxes
    ax8=fig.add_subplot(2,1,2, sharex=ax0)
    plt.setp(ax8.get_xticklabels(), fontsize=8)
    plt.setp(ax8.get_yticklabels(), fontsize=8)
    # uzf recharge
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iRg']],'-', c='darkblue', linewidth=1.5, label = flxLbl[flxIndex_lst['iRg']])
    lines = itertools.cycle(['-','--','-.',':','.',',','o','v','^','<','>','s','p','*','h','H','+','x','D','d','|','_'])
    for i in range(flxIndex_lst['idSg_L1'], len(flxIndex_lst)):
        if np.absolute(sum(flx[i])) > 1E-6:
            plt.plot_date(cMF.inputDate,flx[i], lines.next(), color = mpl.colors.rgb2hex(np.random.rand(1,3)[0]), markersize=2, label = flxLbl[i], markeredgecolor = 'None')
    plt.ylabel('mm', fontsize=10)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 4, columnspacing = colspc, numpoints = 3)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax8.grid(b=True, which='major', axis = 'both')
    ax8.xaxis.grid(b=True, which='minor', color='0.65')
    ax8.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax8.get_xticklabels(minor=True), visible=True)
    plt.xlabel('Date', fontsize=10)
    labels=ax8.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    ax8.xaxis.set_minor_formatter(dateminorFmt)
    labels=ax8.get_xminorticklabels()
    plt.setp(labels, fontsize=8)
    plt.setp(labels, 'rotation', 90)
    del labels

    ax0.set_xlim(np.min(cMF.inputDate)-15.0, np.max(cMF.inputDate)+15)

    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    txt = plt_export_fn.split('.')
    plt_export_fn = txt[0] + '_part3MF.' + txt[1]
    plt.savefig(plt_export_fn,dpi=150)
#    plt.show()
    plt.clf()
    plt.close('all')

##################

def plotTIMESERIES_CATCH(cMF, flx, flxLbl, plt_export_fn, plt_title, hmax, hmin, iniMonthHydroYear, date_ini, date_end, flxIndex_lst, obs_catch = None, obs_catch_list = [0, 0, 0], TopSoilAverage = None, MF = None):
    """
    Plot the time series of the fluxes observed from the whole catchment
    Use Matplotlib
    """

    dateFmt=mpl.dates.DateFormatter('%Y-%b-%d')
    dateminorFmt=mpl.dates.DateFormatter('%b')
    lblspc = 0.05
    mkscale = 1.0
    bdpd = 0.1
    hdltxtpd = 0.05
    colspc = 0.1

    #fmt_DH = mpl.dates.DateFormatter('%Y-%m-%d')
    #cUTIL = MMutils.clsUTILITIES(fmt = fmt_DH)

    #date_ini, year_ini = cUTIL.compDATE_INI(cMF.inputDate[0], iniMonthHydroYear)
    #date_end, year_end = cUTIL.compDATE_END(cMF.inputDate[-1], iniMonthHydroYear)
    date_ini -= 15
    date_end += 15

    fig = plt.figure(num=None, figsize=(8.27, 11.7), dpi = 60)    #(8.5,15), dpi=30)
    fig.suptitle(plt_title)

    ax1=fig.add_subplot(8,1,1)
    ax1.set_autoscalex_on(False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), fontsize=8)
    # RF
    ax1.bar(cMF.inputDate,flx[flxIndex_lst['iRF']],color='darkblue', linewidth=0, align = 'center', label = flxLbl[flxIndex_lst['iRF']])
    # RFe
    ax1.bar(cMF.inputDate,flx[flxIndex_lst['iRFe']],color='deepskyblue', linewidth=0, align = 'center', label = flxLbl[flxIndex_lst['iRFe']])
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    ax1.grid(b=True, which='major', axis = 'both')
    ax1.xaxis.grid(b=True, which='minor', color='0.65')
    ax1.set_xlim(date_ini,date_end)
    plt.ylabel('mm', fontsize=10)
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
    plt.setp(ax1.get_xticklabels(minor=True), visible=False)
    ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))

    ax2=fig.add_subplot(8,1,2, sharex=ax1)
    ax2.set_autoscalex_on(False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    # Ro obs
    if obs_catch_list[2] == 1:
        obs_Ro = obs_catch.get('catch')['obs_Ro']
        Roobs_m = np.ma.masked_values(obs_Ro[0], cMF.hnoflo, atol = 0.09)
        plt.plot_date(cMF.inputDate, Roobs_m, markerfacecolor = 'None', marker='o', markeredgecolor = 'lightblue', markersize=2, label = r'$Ro \ obs$')
        print 'RMSE/RSR/NSE/r of obs. at the catch. scale'
        rmse, rsr, nse, r = cMF.cPROCESS.compCalibCrit(flx[flxIndex_lst['iRo']],obs_Ro[0], cMF.hnoflo)
        rmseRo = [rmse]
        rsrRo = [rsr]
        nseRo = [nse]
        rRo = [r]
        del rmse, rsr, nse, r
        if rmseRo[0] != None:
            print 'Ro: %.1f mm / %.2f / %.2f / %.2f' % (rmseRo[0], rsrRo[0], nseRo[0], rRo[0])           
    # Ro
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iRo']],'r-', c='blue', linewidth=1.0, label = flxLbl[flxIndex_lst['iRo']])
    # Inf
    if np.sum(np.abs(flx[flxIndex_lst['iInf']])) > 1E-7:
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst['iInf']], '--', color = 'brown', label= flxLbl[flxIndex_lst['iInf']])
    # Esurf
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iEsurf']],'r-', c='deepskyblue', linewidth=1.0, label = flxLbl[flxIndex_lst['iEsurf']])
    # Ssurf
    plt.bar(cMF.inputDate, flx[flxIndex_lst['iSsurf']], color='lightblue', linewidth=0, align = 'center', label = flxLbl[flxIndex_lst['iSsurf']])
    # DeltaSsurf
    plt.bar(cMF.inputDate, flx[flxIndex_lst['idSsurf']], color='darkblue', linewidth=0, align = 'center', label = flxLbl[flxIndex_lst['idSsurf']])
    plt.ylabel('mm', fontsize=10)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    ax2.xaxis.grid(b=True, which='minor', color='0.65')
    ax2.grid(b=True, which='major', axis = 'both')
    ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax2.get_xticklabels(minor=True), visible=False)

    ax3=fig.add_subplot(8,1,3, sharex=ax1)
    ax3.set_autoscalex_on(False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), fontsize=8)
    # PE
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iPE']],'-', color='lightblue', linewidth=2.5, label = flxLbl[flxIndex_lst['iPE']])
    if cMF.wel_yn == 1:
        E = flx[flxIndex_lst['iEsoil']] + flx[flxIndex_lst['iEg']]
        # Etot
        plt.plot_date(cMF.inputDate,E,'-', color='darkblue', linewidth=1.5, label = r'$E$')
    # Esoil
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iEsoil']],'--', color='brown', label = flxLbl[flxIndex_lst['iEsoil']])    

    if cMF.wel_yn == 1:
        # Eg
        plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iEg']],'-', color='blue', label = flxLbl[flxIndex_lst['iEg']])
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    ax3.grid(b=True, which='major', axis = 'both')
    ax3.xaxis.grid(b=True, which='minor', color='0.65')
    plt.ylabel('mm', fontsize=10)
    ax3.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax3.get_xticklabels(minor=True), visible=False)

    ax4=fig.add_subplot(8,1,4, sharex=ax1)
    ax4.set_autoscalex_on(False)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), fontsize=8)
    # PT
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iPT']],'-', color='lightblue', linewidth=3, label = flxLbl[flxIndex_lst['iPT']])
    if cMF.wel_yn == 1:
        T = flx[flxIndex_lst['iTsoil']] + flx[flxIndex_lst['iTg']]
        # Ttot
        plt.plot_date(cMF.inputDate,T,'-', color='darkblue',  linewidth=1.5, label = r'$T$')
    # Tsoil
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iTsoil']],'--', color='brown', label = flxLbl[flxIndex_lst['iTsoil']])    
    if cMF.wel_yn == 1:
        # Tg
        plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iTg']],'-', color='blue', label = flxLbl[flxIndex_lst['iTg']])
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax4.grid(b=True, which='major', axis = 'both')
    ax4.xaxis.grid(b=True, which='minor', color='0.65')
    plt.ylabel('mm', fontsize=10)
    ax4.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax4.get_xticklabels(minor=True), visible=False)

    ax5=fig.add_subplot(8,1,5, sharex=ax1)
    ax5.set_autoscalex_on(False)
    plt.setp(ax5.get_xticklabels(), fontsize=8)
    plt.setp(ax5.get_yticklabels(), fontsize=8)
    # Rp
    ax5.plot_date(cMF.inputDate, -1.0*flx[flxIndex_lst['iperc']], '-', color = 'brown', label= flxLbl[flxIndex_lst['iRsoil']])
    if cMF != None:
        # Rg
        plt.plot_date(cMF.inputDate,-1.0*flx[flxIndex_lst['iRg']],'-', c='darkblue', linewidth=1.5, label = flxLbl[flxIndex_lst['iRg']])
    if cMF.wel_yn == 1:
        # ETg
        plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iETg']],'-', c='blue', linewidth=1.5, label = flxLbl[flxIndex_lst['iETg']])
    # EXF
    plt.bar(cMF.inputDate,flx[flxIndex_lst['iEXFg']], color='lightblue', linewidth=0, align = 'center', label = flxLbl[flxIndex_lst['iEXFg']])
    plt.legend(loc = 0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax5.grid(b = True, which = 'major', axis = 'both')
    ax5.xaxis.grid(b=True, which='minor', color='0.65')
    plt.xlabel('Date', fontsize=10)    
    plt.ylabel('mm', fontsize=10)
    ax5.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax5.get_xticklabels(minor=True), visible=True)
    labels=ax5.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    ax5.xaxis.set_minor_formatter(dateminorFmt)
    labels=ax5.get_xminorticklabels()
    plt.setp(labels, fontsize=8)
    plt.setp(labels, 'rotation', 90)
    del labels

    ax6=fig.add_subplot(8,1,7, sharex=ax1)
    ax6.set_autoscalex_on(False)
    plt.setp(ax6.get_yticklabels(), fontsize=8)
    # theta
    rmseSM = None
    rsrSM = None
    nseSM = None
    rSM = None
    if obs_catch_list[1] == 1:
        obs_SM = obs_catch.get('catch')['obs_SM']
        Sobs_m = np.ma.masked_values(obs_SM[0], cMF.hnoflo, atol = 0.09)
        ax6.plot_date(cMF.inputDate, Sobs_m, markerfacecolor = 'None', marker='o', markeredgecolor = 'brown', markersize=2, label = r'$\theta \ obs$') 
        rmse, rsr, nse, r = cMF.cPROCESS.compCalibCrit(flx[flxIndex_lst['iSsoil_pc']],obs_SM[0], cMF.hnoflo)
        rmseSM = [100.0*rmse]
        rsrSM = [rsr]
        nseSM = [nse]
        rSM = [r]
        del rmse, rsr, nse, r
        if rmseSM[0] != None:
            print 'SM: %.1f %% / %.2f / %.2f / %.2f' % (rmseSM[0], rsrSM[0], nseSM[0], rSM[0])     
    ax6.plot_date(cMF.inputDate, flx[flxIndex_lst['iSsoil_pc']], '-', color = 'brown', label = flxLbl[flxIndex_lst['iSsoil_pc']])
    # y axis
    plt.ylabel('%', fontsize=10)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, numpoints = 3)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    if cMF == None:
        plt.xlabel('Date', fontsize=10)
        labels=ax6.get_xticklabels()
        plt.setp(labels, fontsize=8)
        plt.setp(labels, 'rotation', 90)
        ax6.xaxis.set_minor_formatter(dateminorFmt)
        labels=ax6.get_xminorticklabels()
        plt.setp(labels, fontsize=8)
        plt.setp(labels, 'rotation', 90)
        del labels
    else:
        plt.setp(ax6.get_xticklabels(), visible=False)
        plt.setp(ax6.get_xticklabels(minor=True), visible=False)
    ax6.grid(b=True, which='major', axis = 'both')
    ax6.xaxis.grid(b=True, which='minor', color='0.65')
    ax6.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))

    if MF != None:
        # compute heads
        if obs_catch_list[0] == 1:
            obs_h = obs_catch.get('catch')['obs_h']
            hobs_m = np.ma.masked_values(obs_h[0], cMF.hnoflo, atol = 0.09)
        # plot GWT
        lines = itertools.cycle(['-','--','-.',':','.',',','o','v','^','<','>','s','p','*','h','H','+','x','D','d','|','_'])
        ax10=fig.add_subplot(8,1,8, sharex=ax1)
        ax10.set_autoscalex_on(False)
        plt.setp(ax10.get_xticklabels(), fontsize=8)
        plt.setp(ax10.get_yticklabels(), fontsize=8)
        if obs_catch_list[0] == 1:
            dobs_m = hobs_m - TopSoilAverage
            ax10.plot_date(cMF.inputDate, dobs_m, markerfacecolor = 'None', marker='o', markeredgecolor = 'LightBlue', markersize=2, label = r'$d \ obs$')
        for l in range(cMF.nlay):
            i = 'id_L%d'%(l+1)
            plt.plot_date(cMF.inputDate,flx[flxIndex_lst[i]],lines.next(), color = 'b', label = flxLbl[flxIndex_lst[i]])
        plt.ylabel('m', fontsize=10)
        plt.xlabel('Date', fontsize=10)
        plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize=8)
        ax10.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))        
        ax10.grid(b=True, which='major', axis = 'both')
        labels=ax10.get_xticklabels()
        plt.setp(labels, 'rotation', 90)
        ax10.xaxis.set_minor_formatter(dateminorFmt)
        labels=ax10.get_xminorticklabels()
        plt.setp(labels, fontsize=8)
        plt.setp(labels, 'rotation', 90)
        del labels
        ax10.xaxis.grid(b=True, which='minor', color='0.65')     
        plt.setp(ax10.get_xticklabels(minor=True), visible=True)

    ax1.set_xlim(date_ini,date_end)

    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig(plt_export_fn,dpi=150)

#--------------FIGURE CATCH PART2----------

    fig = plt.figure(num=None, figsize=(8.27, 11.7), dpi = 60)    #(8.5,15), dpi=30)
    fig.suptitle(plt_title + ' - part 2')
    
    ax0=fig.add_subplot(8,1,1)
    ax0.set_autoscalex_on(False)
    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.setp(ax0.get_yticklabels(), fontsize=8)
    # RF
    ax0.bar(cMF.inputDate,flx[flxIndex_lst['iRF']],color='darkblue', linewidth=0, align = 'center', label = flxLbl[flxIndex_lst['iRF']])
    # RFe
    ax0.bar(cMF.inputDate,flx[flxIndex_lst['iRFe']],color='deepskyblue', linewidth=0, align = 'center', label = flxLbl[flxIndex_lst['iRFe']])
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    ax0.grid(b=True, which='major', axis = 'both')
    ax0.xaxis.grid(b=True, which='minor', color='0.65')
    plt.ylabel('mm', fontsize=10)
    ax0.xaxis.set_major_formatter(dateFmt)
    ax0.xaxis.set_major_locator(mpl.dates.YearLocator(1, month = iniMonthHydroYear, day = 1))
    ax0.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth = bymonth))
    plt.setp(ax0.get_xticklabels(minor=True), visible=False)
    ax0.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))    

    rmseHEADS = None
    rsrHEADS = None
    nseHEADS = None
    rHEADS = None
    if MF != None:
        # plot heads
        lines = itertools.cycle(['-','--','-.',':','.',',','o','v','^','<','>','s','p','*','h','H','+','x','D','d','|','_'])
        ax1=fig.add_subplot(8,1,4)
        ax1.set_autoscalex_on(True)
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.setp(ax1.get_yticklabels(), fontsize=8)
        # RMSE
        if obs_catch_list[0] == 1:
            ax1.plot_date(cMF.inputDate, hobs_m, markerfacecolor = 'None', marker='o', markeredgecolor = 'LightBlue', markersize=2, label = r'$h \ obs$')
            i = 'ih_L%d' % (l+1) 
            if sum(flx[flxIndex_lst[i]]) != 0.0:
                rmse, rsr, nse, r = cMF.cPROCESS.compCalibCrit(flx[flxIndex_lst[i]],obs_h[0], cMF.hnoflo)
                rmseHEADS = [rmse]
                rsrHEADS = [rsr]
                nseHEADS = [nse]
                rHEADS = [r]
                del rmse, rsr, nse, r
                if rmseHEADS[0] != None:
                    print 'h: %.2f m / %.2f / %.2f / %.2f\n-------' % (rmseHEADS[0], rsrHEADS[0], nseHEADS[0], rHEADS[0])             
            else:
                print 'Warning!\nError in computing h calibration criteria'
                rmseHEADS = rsrHEADS = nseHEADS = rHEADS = None
        for l in range(cMF.nlay):
            i = 'ih_L%d' % (l+1) 
            plt.plot_date(cMF.inputDate,flx[flxIndex_lst[i]],lines.next(), color = 'b', label = flxLbl[flxIndex_lst[i]])
        plt.ylim(hmin,hmax)
        plt.ylabel('m', fontsize=10)
        plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize=8 )
        ax1.grid(b=True, which='major', axis = 'both')
        ax1.xaxis.grid(b=True, which='minor', color='0.65')
        ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
        plt.setp(ax1.get_xticklabels(minor=True), visible=False)

        # plot GW fluxes
        ax8=fig.add_subplot(2,1,2, sharex=ax1)
        ax8.set_autoscalex_on(True)
        plt.setp(ax8.get_xticklabels(), fontsize=8)
        plt.setp(ax8.get_yticklabels(), fontsize=8)
        # uzf recharge
        plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iRg']],'-', c='darkblue', linewidth=1.5, label = flxLbl[flxIndex_lst['iRg']])
        lines = itertools.cycle(['-','--','-.',':','.',',','o','v','^','<','>','s','p','*','h','H','+','x','D','d','|','_'])
        for l, (e, lbl) in enumerate(zip(flx[flxIndex_lst['idSg_L1']:], flxLbl[flxIndex_lst['idSg_L1']:])):
            plt.plot_date(cMF.inputDate,e, lines.next(), color = mpl.colors.rgb2hex(np.random.rand(1,3)[0]), markersize=2, label = lbl, markeredgecolor = 'None')
        plt.ylabel('mm', fontsize=10)
        plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 4, columnspacing = colspc, numpoints = 3)
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize=8 )
        ax8.grid(b=True, which='major', axis = 'both')
        ax8.xaxis.grid(b=True, which='minor', color='0.65')
        ax8.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
        plt.setp(ax8.get_xticklabels(minor=True), visible=True)
        plt.xlabel('Date', fontsize=10)
        labels=ax8.get_xticklabels()
        plt.setp(labels, 'rotation', 90)
        ax8.xaxis.set_minor_formatter(dateminorFmt)
        labels=ax8.get_xminorticklabels()
        plt.setp(labels, fontsize=8)
        plt.setp(labels, 'rotation', 90)
        del labels

    # plot Ro
    ax2=fig.add_subplot(8,1,2, sharex=ax1)
    ax2.set_autoscalex_on(True)
    plt.setp(ax2.get_xticklabels(), fontsize=8)
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    # Ro obs
    if obs_catch_list[2] == 1:
        obs_Ro = obs_catch.get('catch')['obs_Ro']
        Roobs_m = np.ma.masked_values(obs_Ro[0], cMF.hnoflo, atol = 0.09)
        plt.plot_date(cMF.inputDate, Roobs_m, markerfacecolor = 'None', marker='o', markeredgecolor = 'lightBlue', markersize=2, label = r'$Ro \ obs$')
    # Ro
    plt.plot_date(cMF.inputDate,flx[flxIndex_lst['iRo']],'r-', c='blue', linewidth=1.0, label = flxLbl[5])
    plt.ylabel('mm', fontsize=10)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad = bdpd, handletextpad = hdltxtpd, ncol = 2, columnspacing = colspc, numpoints = 3)
    if obs_catch_list[2] == 1:
        plt.ylim(ymax = np.ma.max(Roobs_m))
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    ax2.xaxis.grid(b=True, which='minor', color='0.65')
    ax2.grid(b=True, which='major', axis = 'both')
    ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax2.get_xticklabels(minor=True), visible=True)
    plt.xlabel('Date', fontsize=10)
    labels=ax2.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    ax2.xaxis.set_minor_formatter(dateminorFmt)
    labels=ax2.get_xminorticklabels()
    plt.setp(labels, fontsize=8)
    plt.setp(labels, 'rotation', 90)
    del labels

    ax0.set_xlim(np.min(cMF.inputDate)-15.0, np.max(cMF.inputDate)+15)

    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    txt = plt_export_fn.split('.')
    plt_export_fn = txt[0] + '_part2.' + txt[1]
    plt.savefig(plt_export_fn,dpi=150)
#    plt.show()
    plt.clf()
    plt.close('all')
    return rmseHEADS, rmseSM, rsrHEADS, rsrSM, nseHEADS, nseSM, rHEADS, rSM

##################

def plotLAYER(days, str_per, Date, JD, ncol, nrow, nlay, nplot, V, cmap, CBlabel, msg, plt_title, MM_ws, interval_type = 'arange', interval_diff = 1, interval_num = 1, Vmax = 0, Vmin = 0, fmt = None, contours = False, ntick = 1, axisbg = 'silver', points  = None, ptslbl = 0, mask = None, hnoflo = -999.9, animation = 0, pref_plt_title = '_sp_plt', cMF = None):

    # TODO put option to select axes tick as row/col index from MODFLOW or real coordinates (in this last case create it)

    def MinMax(min_, max_, ctrs_):
        if max_ == min_:
            if max_ < 10E-9:
                max_ = 1.0
                min_ = -1.0
            else:
                max_ *= 1.15
                min_ *= 0.85
            ctrs_ = False
        elif max_ < min_:
            tmp = min_
            min_ = max_
            max_ = tmp
            ctrs_ = ctrs_
        else:
            ctrs_ = ctrs_
        return min_, max_, ctrs_

    # Store some arrays for plotting
    x = np.arange(0.5, ncol+1.5, 1)
    y = np.arange(0.5, nrow+1.5, 1)
    xg,yg = np.meshgrid(x,y)

    x = np.arange(1, ncol+1, 1)
    y = np.arange(1, nrow+1, 1)
    xg1,yg1 = np.meshgrid(x,y)

    ims = []
    files_tmp = []
    for i, day in enumerate(days):
        ax= []
        fig = plt.figure(num=None, figsize=(11.7, 8.27), dpi=30)
        figtitle = fig.suptitle('')
        ims.append([])
        if isinstance(Date[i], float):
            figtitle.set_text(plt_title + '\nDate %s, DOY %s, stress period %s, day %d' % (mpl.dates.num2date(Date[i]).isoformat()[:10], JD[i], str_per[i], day+1))
        else:
            figtitle.set_text(plt_title)
        plt.draw()
        for L in range(nplot):
            if mask == None:
                Vtmp = V[i,L,:,:]
            else:
                Vtmp = np.ma.masked_array(V[i,L,:,:], mask[L,:,:])
            Vtmp = np.ma.masked_values(Vtmp, hnoflo, atol = 0.09)
            Vmin_tmp, Vmax_tmp, ctrs_tmp = MinMax(Vmin[i], Vmax[i], contours)
            if fmt == None:
                if Vmax_tmp > 0.0999 or abs(Vmin_tmp)> 0.0999:
                    fmt = '%5.2f'
                else:
                    fmt = '%5.e'
            norm = None
            if interval_type == 'arange':
                ticks = np.arange(Vmin_tmp,Vmax_tmp,interval_diff)
            elif interval_type == 'linspace':
                ticks = np.linspace(Vmin_tmp,Vmax_tmp,interval_num)
            elif interval_type == 'percentile':
                ticks = np.percentile(Vtmp.compressed().flatten(),[0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0])
                norm = mpl.colors.BoundaryNorm(ticks, cmap.N)
            ax.append(fig.add_subplot(1,nlay,L+1, axisbg = axisbg))
            ax[L].xaxis.set_ticks(np.arange(0,ncol+1,ntick))
            ax[L].yaxis.set_ticks(np.arange(0,nrow+1,ntick))
            plt.setp(ax[L].get_xticklabels(), fontsize=8)
            plt.setp(ax[L].get_yticklabels(), fontsize=8)
            plt.ylabel('row i', fontsize=10)
            plt.xlabel('col j', fontsize=10)
            if points <> None:
                for k, (xj,yi,lay, label) in enumerate(zip(points[2],points[1],points[3],points[0])):
                    if lay == L:
                        color = 'dimgrey'
                    else:
                        color = 'lightgrey'
                    ax[L].plot(xj, yi, 'o', linewidth=1, markersize = 4, color = color)
                    if ptslbl>0:
                        ax[L].annotate(label, xy = (xj, yi-0.15), fontsize=8, ha ='center', va = 'bottom')
            if Vmax_tmp >0 and Vmin_tmp<0 and cMF != None:
                cmap = plt.cm.coolwarm_r
                #shifted_cmap = cMF.cUTIL.remappedColorMap(cmap, midpoint=0.75, name='shifted')
                start = 0.0 #(Vmax_tmp-abs(Vmin_tmp))/(2*Vmax_tmp)
                midpoint = abs(Vmin_tmp)/(Vmax_tmp+abs(Vmin_tmp))
                stop = 1.0 #(abs(Vmin_tmp)-Vmax_tmp)/(2*abs(Vmin_tmp))
                shrunk_cmap = cMF.cUTIL.remappedColorMap(cmap, start=start, midpoint=midpoint, stop=stop, name='shrunk')
                cmap = shrunk_cmap                            
            ims[i].append(ax[L].pcolormesh(xg, yg, Vtmp, cmap = cmap, vmin = Vmin_tmp, vmax = Vmax_tmp, norm = norm))
            if ctrs_tmp == True:
                CS = ax[L].contour(xg1, yg1[::-1], Vtmp[::-1], ticks, colors = 'gray')
                plt.draw()
                ax[L].clabel(CS, inline=1, fontsize = 6, fmt=fmt, colors = 'gray')
            if np.ma.max(Vtmp)>np.ma.min(Vtmp):
                ax[L].set_title('layer %d' % (L+1), fontsize = 10)
            else:
                ax[L].set_title('layer %d %s' % (L+1, msg), fontsize = 10)
            ax[L].set_ylim(bottom = np.max(yg1), top = np.min(yg1))
            ax[L].axis('scaled')
            if interval_type == 'percentile':
                if max(x) > max(y):
                    CBorient = 'horizontal'
                else:
                    CBorient = 'vertical'
                CB = fig.colorbar(ims[i][0+L], extend='both', ticks = ticks, format = fmt, orientation = CBorient)
                CB.set_label(CBlabel, fontsize = 10)
                plt.setp(CB.ax.get_xticklabels(), fontsize = 7)
                plt.setp(CB.ax.get_yticklabels(), fontsize = 7)

        if interval_type != 'percentile':
            if max(x) > max(y):
                cax = fig.add_axes([.125, 0.035, 0.75, 0.025])
                CBorient = 'horizontal'
            else:
                cax = fig.add_axes([0.035, 0.125, 0.025, 0.75])
                CBorient = 'vertical'
            CB = fig.colorbar(ims[i][0], extend='both', ticks = ticks, format = fmt, cax = cax,  orientation = CBorient)
            CB.set_label(CBlabel, fontsize = 12)
            if max(x) > max(y):
                cax.xaxis.set_label_position('top')
                plt.setp(CB.ax.get_xticklabels(), fontsize = 7)
            else:
                cax.yaxis.set_label_position('left')
                plt.setp(CB.ax.get_yticklabels(), fontsize = 7)                
            
        if isinstance(Date[i], float):
            plt_export_fn = os.path.join(MM_ws, '%s_%s_day%05d.png' % (pref_plt_title, plt_title, day+1))
        else:
            plt_export_fn = os.path.join(MM_ws, '%s_%s.png' % (pref_plt_title, plt_title))
        plt.savefig(plt_export_fn)
        if len(days)>1 and animation == 1:
            try:
                if i == 0:
                    files_tmp.append(os.path.join(MM_ws,'%05d.png'%(0)))
                    print plt_export_fn, files_tmp[0]
                    shutil.copyfile(plt_export_fn, files_tmp[0])
                else:
                    files_tmp.append(os.path.join(MM_ws,'%05d.png'%(i)))
                    shutil.copyfile(plt_export_fn, files_tmp[i])
            except:
                pass
        for L in range(nplot):
            ax[L].cla()

    if len(days)>1 and animation == 1:
        batch_fn = os.path.join(MM_ws, 'run.bat')
        f = open(batch_fn, 'w')
        f.write('ffmpeg -r 1 -i %s -s:v 1280x720 -c:v libx264 -profile:v high -crf 23 -pix_fmt yuv420p -r 30 -y %s_mov.mp4' % ('%s\%%%%05d.png' % (MM_ws),'%s\%s_%s' % (MM_ws, pref_plt_title, plt_title)))
        f.close()
        run_report_fn = os.path.join(MM_ws, '__FFmpegRunReport.txt')
        run_report = open(run_report_fn, 'w')
        sp.Popen(batch_fn, shell=False, stdout = run_report, stderr = run_report).wait()
        run_report.close()
        os.remove(batch_fn)
        try:
            for f in files_tmp:
                os.remove(f)
        except:
            pass
    
    plt.close()

##################

def plotWBsankey(path, DATE, flx, flxIndex, fn, indexTime, year_lst, cMF, ncell_MM, obspt, fntitle, ibound4Sankey, stdout = None, report = None):

    """ Computes the water balance for a certain time span
    input: ASCII file with water fluxes wrtitten by MM
    """

    # compute fluxes for whole modelled period and hydrological years
    RF=[]
    I=[]
    RFe=[]
    dSsurf=[]
    Ro=[]
    Esurf=[]
    dSsoil=[]
    EXFtotMM=[]
    Exf_l0 = []
    Esoil=[]
    Tsoil=[]
    ETsoil=[]
    Inf = []
    Eg = []
    Tg = []
    Egtot=[]
    Tgtot=[]
    ETg=[]
    Ssurf=[]
    Rp=[]
    dSu=[]
    Rg=[]
    dSg=[]
    FRF=[]
    FFF=[]
    FLF=[]
    EXF=[]
    EXFtotMF=[]
    WEL=[]
    DRN=[]
    GHB=[]
    for k, i in enumerate(indexTime[:-2]):
        if k == 0:
            i = indexTime[1]
            indexend = indexTime[-2]
            mult = 365.0/(indexend-i+1)
        else:
            indexend = indexTime[k+1]-1
            mult = 1.0
        RF.append(mult*np.sum(np.float16(flx[flxIndex['iRF']][i:indexend])))
        I.append(mult*np.sum(np.float16(flx[flxIndex['iI']][i:indexend])))
        RFe.append(mult*np.sum(np.float16(flx[flxIndex['iRFe']][i:indexend])))
        dSsurf.append(mult*np.sum(np.float16(flx[flxIndex['idSsurf']][i:indexend])))
        Ro.append(mult*np.sum(np.float16(flx[flxIndex['iRo']][i:indexend])))
        Esurf.append(mult*np.sum(np.float16(flx[flxIndex['iEsurf']][i:indexend])))
        dSsoil.append(mult*np.sum(np.float16(flx[flxIndex['idSsoil']][i:indexend])))
        EXFtotMM.append(mult*np.sum(np.float16(flx[flxIndex['iEXFg']][i:indexend])))
        Exf_l0.append(mult*np.sum(np.float16(flx[flxIndex['iExf_l1']][i:indexend])))
        Esoil.append(mult*np.sum(np.float16(flx[flxIndex['iEsoil']][i:indexend])))
        Tsoil.append(mult*np.sum(np.float16(flx[flxIndex['iTsoil']][i:indexend])))
        ETsoil.append(mult*np.sum(np.float16(flx[flxIndex['iETsoil']][i:indexend])))
        Inf.append(mult*np.sum(np.float16(flx[flxIndex['iInf']][i:indexend])))
        if cMF.wel_yn == 1:
            Eg.append([])
            Tg.append([])
            for L in range(cMF.nlay):
                Eg[k].append(mult*np.sum(np.float16(flx[flxIndex['iEg_L%d'%(L+1)]][i:indexend])))
                Tg[k].append(mult*np.sum(np.float16(flx[flxIndex['iTg_L%d'%(L+1)]][i:indexend])))
            Egtot.append(mult*np.sum(np.float16(flx[flxIndex['iEg']][i:indexend])))
            Tgtot.append(mult*np.sum(np.float16(flx[flxIndex['iTg']][i:indexend])))
            ETg.append(mult*np.sum(np.float16(flx[flxIndex['iETg']][i:indexend])))
        Ssurf.append(mult*np.sum(np.float16(flx[flxIndex['iSsurf']][i:indexend])))
        Rp.append(mult*np.sum(np.float16(flx[flxIndex['iperc']][i:indexend])))
        dSu.append(mult*np.sum(np.float16(flx[flxIndex['idSu']][i:indexend])))
        Rg.append([])
        dSg.append([])
        FRF.append([])
        FFF.append([])
        FLF.append([])
        EXF.append([])
        WEL.append([])
        DRN.append([])
        GHB.append([])
        for L in range(cMF.nlay):
            Rg[k].append(mult*np.sum(np.float16(flx[flxIndex['iRg_L%d'%(L+1)]][i:indexend])))
            dSg[k].append(mult*np.sum(np.float16(flx[flxIndex['idSg_L%d'%(L+1)]][i:indexend])))
            FRF[k].append(mult*np.sum(np.float16(flx[flxIndex['iFRF_L%d'%(L+1)]][i:indexend])))
            FFF[k].append(mult*np.sum(np.float16(flx[flxIndex['iFFF_L%d'%(L+1)]][i:indexend])))
            if cMF.nlay>1:
                FLF[k].append(mult*np.sum(np.float16(flx[flxIndex['iFLF_L%d'%(L+1)]][i:indexend])))
            EXF[k].append(mult*np.sum(np.float16(flx[flxIndex['iEXFg_L%d'%(L+1)]][i:indexend])))
            if cMF.wel_yn == 1:
                if ncell_MM[L]>0:
                    WEL[k].append(mult*np.sum(np.float16(flx[flxIndex['iWEL_L%d'%(L+1)]][i:indexend])))
                else:
                    WEL[k].append(0)
            if cMF.drn_yn == 1:
                if cMF.drncells[L]>0:
                    DRN[k].append(mult*np.sum(np.float16(flx[flxIndex['iDRN_L%d'%(L+1)]][i:indexend])))
                else:
                    DRN[k].append(0)
            if cMF.ghb_yn == 1:
                if cMF.ghbcells[L] > 0:
                    GHB[k].append(mult*np.sum(np.float16(flx[flxIndex['iGHB_L%d'%(L+1)]][i:indexend])))
                else:
                    GHB[k].append(0)
        EXFtotMF.append(sum(EXF[k]))
#    print "\nWater fluxes imported from file:\n%s" % inputFile_fn

    # Sankey plots
#    prt_test = 0
    for f in [0,1]:
        for k in range(len(RF)):
            #print '-------'
            if k == 0:
                title = "Average of the %d hydrological year(s)" % len(indexTime[1:-2])
            else:
                title = "Hydrological year %d/%d" % (year_lst[k-1], year_lst[k-1] + 1)
            if f == 0:
                ff = 1.0
                fff = (1.5*RF[k])
#                if prt_test == 0:
#                    print '\nWater balance (mm.y-1)'
#                    prt_test += 1
            else:
                ff = RF[k]/100.0
                fff = (1.5*RF[k])/ff
#                if prt_test == 1:
#                    print '\nWater balance (%)'
#                    prt_test += 1
#            lbl_tmp = ''
#            if cMF.wel_yn == 1:
#                lbl_tmp += 'WEL %s' % (np.asarray(WEL[k])/ff)
#            if cMF.drn_yn == 1:
#                lbl_tmp += 'DRN %s' % (np.asarray(DRN[k])/ff)
#            if cMF.ghb_yn == 1:
#                lbl_tmp += 'GHB %s' % (np.asarray(GHB[k])/ff)
            # print '\nMMsurf: RF %s, I %s, RFe %s, Esurf %s, DSsurf %s\nMMsoil: Esoil %s, Tsoil %s, ETsoil %s, Ro %s, Rp %s, DSsoil  %s\nUZF: Rg %s, DSu %s\nMF: Eg %s, Tg %s, ETg %s, %s, EXF  %s, FLF %s, dSg %s' % (RF[k]/ff, I[k]/ff, RFe[k]/ff, Esurf[k]/ff, DSsurf[k]/ff, Esoil[k]/ff, Tsoil[k]/ff, ETsoil[k]/ff, Ro[k]/ff, Rp[k]/ff, DSsoil[k]/ff, (np.asarray(Rg[k])/ff), DSu[k]/ff, (np.asarray(Eg[k])/ff), (np.asarray(Tg[k])/ff), ETg[k]/ff, lbl_tmp, (np.asarray(EXF[k])/ff), (np.asarray(FLF[k])/ff), (np.asarray(dSg[k])/ff))
            treshold = 5E-2
            if k == 0 or k%2 != 0:
                fig = plt.figure(figsize=(8.27, 11.7), dpi = 72)
                if f == 0:
                    fig.suptitle('Water balance ($mm.y^{-1}$) - %s\n'%obspt, fontsize = 10, y = 0.99)
                else:
                    fig.suptitle('Water balance ($\%%$ of yearly rainfall) - %s\n'%obspt, fontsize = 10, y = 0.99)
                ax = []
                ax.append(fig.add_subplot(2, 1, 1, xticks=[], yticks=[]))
                y = 0.530
                p = 0
            else:
                ax.append(fig.add_subplot(2, 1, 2, xticks=[], yticks=[]))
                y = 0.085
                p = 1
            figtitle = fig.text(x = 0.5, y = y, s = title, horizontalalignment = 'center', verticalalignment = 'bottom', fontsize = 8)
            fmt = '%.1f'
            pltsankey = Sankey(ax=ax[p], format=fmt, scale=1.0/fff, offset = 0.25, gap = 0.5, shoulder = 0.0, margin = 0.5)
            pl = 0.5
            tl = 2.0
            # MMveg
            pltsankey.add(patchlabel='veg', facecolor='lightgreen', trunklength =tl/1.5,
                       flows=[RF[k]/ff, -I[k]/ff, -RFe[k]/ff],
                       labels=['$RF$', '$I$', '$RFe$'],
                       orientations=[1,1,-1],
                       pathlengths = [pl, pl, pl])
            # MMsurf
            if np.abs(Exf_l0[k]) > treshold:
               flows=[RFe[k]/ff, -Inf[k]/ff, Exf_l0[k]/ff, -Esurf[k]/ff, -Ro[k]/ff]
               labels=[None, '$Inf$', '$Exf^{l1}$', '$E_{surf}$','$Ro$']
               orientations=[1,-1,-1,1,0]
               pathlengths = [pl, pl, pl, pl, pl]
            else:
               flows=[RFe[k]/ff, -Inf[k]/ff, -Esurf[k]/ff, -Ro[k]/ff]
               labels=[None, '$Inf$', '$E_{surf}$','$Ro$']
               orientations=[1,-1,1,0]
               pathlengths = [pl, pl, pl, pl]
            pltsankey.add(patchlabel = '$\Delta S_{surf}$\n%.1f' % (-dSsurf[k]/ff), label='MMsurf', facecolor='lightblue', trunklength = tl,
               flows=flows,
               labels=labels,
               orientations=orientations,
               pathlengths = pathlengths,
               prior=0, connect=(2,0))
            In = RFe[k] + Exf_l0[k]
            Out = Inf[k]  + Esurf[k] + Ro[k]
            if  dSsurf[k] > 0.0:
                Out += dSsurf[k]
            else:
                In += -dSsurf[k]
            MB_MMsurf = 100*(In - Out)/((In + Out)/2)
            # MMsoil
            if EXFtotMM[k] > treshold:
                if np.abs(-Exf_l0[k]) > treshold:
                    flows = [Inf[k]/ff, -Rp[k]/ff, -Exf_l0[k]/ff, -Esoil[k]/ff, -Tsoil[k]/ff, EXFtotMM[k]/ff]
                    labels = [None,'$R_p$','$Exf^{l1}$','$E_{soil}$','$T_{soil}$','$Exf_g$']
                    orientations = [1,-1,1,1,1,-1]
                    pathlengths = [pl, pl, pl, pl, pl, pl]
                else:
                    flows = [Inf[k]/ff, -Rp[k]/ff, -Esoil[k]/ff, -Tsoil[k]/ff, EXFtotMM[k]/ff]
                    labels = [None,'$R_p$','$E_{soil}$','$T_{soil}$','$Exf_g$']
                    orientations = [1,-1,1,1,-1]
                    pathlengths = [pl, pl, pl, pl, pl]
            else:
                if np.abs(EXF[k][0]) > treshold:
                    flows = [Inf[k]/ff, -Rp[k]/ff, -Exf_l0[k]/ff, -Esoil[k]/ff, -Tsoil[k]/ff]
                    labels = [None,'$R_p$','$Exf^{l1}$','$E_{soil}$','$T_{soil}$']
                    orientations = [1,-1,1,1,1]
                    pathlengths = [pl, pl, pl, pl, pl]
                else:
                    flows = [Inf[k]/ff,-Rp[k]/ff, -Esoil[k]/ff, -Tsoil[k]/ff]
                    labels = [None,'$R_p$','$E_{soil}$','$T_{soil}$']
                    orientations = [1,-1,1,1]
                    pathlengths = [pl, pl, pl, pl]
            pltsankey.add(patchlabel = '$\Delta S_{soil}$\n%.1f' % (dSsoil[k]/ff),
                       label='MMsoil', facecolor='khaki', trunklength = tl,
                       flows=flows,
                       labels=labels,
                       orientations=orientations,
                       pathlengths = pathlengths,
                       prior=1, connect=(1,0))
            In = Inf[k] + EXFtotMM[k]
            Out = Rp[k] + Esoil[k] + Tsoil[k] + Exf_l0[k]
            if  dSsoil[k] > 0.0:
                Out += dSsoil[k]
            else:
                In += -dSsoil[k]
            MB_MMsoil = 100*(In - Out)/((In + Out)/2)
            # MFuzf
            flows=[Rp[k]/ff]
            labels=[None]
            orientations=[1]
            pathlengths = [2*pl]
            for L in range(cMF.nlay):
                if ibound4Sankey[L] > 0:
                    if Rg[k][L]/ff>0.0:
                        flows.append(-Rg[k][L]/ff)
                    else:
                        flows.append(0.0)                        
                    labels.append('$Rg^{L%d}$'%(L+1))
                    orientations.append(-1)
                    pathlengths.append(pl)
            pltsankey.add(patchlabel = '$\Delta S_p$\n%.1f' % (dSu[k]/ff), label='MF_UZF', facecolor='lavender', trunklength = tl,
                       flows = flows,
                       labels=labels,
                       orientations=orientations,
                       pathlengths = pathlengths,
                       prior=2, connect=(1, 0))
            In = Rp[k]
            Out = 0
            for L in range(cMF.nlay):
                Out += Rg[k][L]
            if  dSu[k] > 0.0:
                Out += dSu[k]
            else:
                In += -dSu[k]
            MB_MFuzf = 100*(In - Out)/((In + Out)/2)
            # MF
            # about signs, read MF-2005 manual pag 3-10
            MB_MF = []
            #print '%s'%obspt
            L_act = 0
            for L in range(cMF.nlay):
                if ibound4Sankey[L] > 0:
                    L_act += L
                    In = []
                    Out = []
                    if L_act == 0:
                        if L == (cMF.nlay-1):
                            flows=[Rg[k][L]/ff, -FRF[k][L]/ff, -FFF[k][L]/ff]
                            labels=[None, '$FRF$', '$FFF$']
                            orientations=[1, 0, 0]
                            pathlengths = [pl*4, pl, pl*4]
                        else:
                            flows=[Rg[k][L]/ff, -FLF[k][L]/ff, -FRF[k][L]/ff, -FFF[k][L]/ff]
                            labels=[None, '$FLF$', '$FRF$', '$FFF$']
                            orientations=[1, -1, 0, 0]
                            pathlengths = [pl*4, pl, pl, pl*4]
                        if cMF.nlay>1:
                            if FLF[k][L] > 0.0:
                                Out.append(FLF[k][L])
                            else:
                                In.append(-FLF[k][L])
                    elif L_act == (cMF.nlay-1):
                        flows=[FLF[k][L-1]/ff, Rg[k][L]/ff, -FRF[k][L]/ff, -FFF[k][L]/ff]
                        labels=[None, '$Rg^{L%d}$'%(L+1), '$FRF$', '$FFF$']
                        orientations=[1, 1, 0, 0]
                        pathlengths = [pl*4, pl, pl, pl*4]
                        if FLF[k][L-1] > 0.0:
                            In.append(FLF[k][L-1])
                        else:
                            Out.append(-FLF[k][L-1])
                    else:
                        flows=[FLF[k][L-1]/ff, Rg[k][L]/ff, -FLF[k][L]/ff, -FRF[k][L]/ff, -FFF[k][L]/ff]
                        labels=[None, None, '$\Delta S_g$', '$FLF$', '$FRF$', '$FFF$']
                        orientations=[1, 1, -1, 0, 0]
                        pathlengths = [pl, pl, pl, pl, pl*4]
                        if FLF[k][L] > 0.0:
                            Out.append(FLF[k][L])
                        else:
                            In.append(-FLF[k][L])
                        if FLF[k][L-1] > 0.0:
                            In.append(FLF[k][L-1])
                        else:
                            Out.append(-FLF[k][L-1])
                    # Rg
                    In.append(Rg[k][L])
                    # FRF                
                    if FRF[k][L] > 0.0:
                        Out.append(FRF[k][L])
                    else:
                        In.append(-FRF[k][L])
                    # FFF
                    if FFF[k][L] > 0.0:
                        Out.append(FFF[k][L])
                    else:
                        In.append(-FFF[k][L])
                    # GW STO
                    if  dSg[k][L] > 0.0:
                        In.append(dSg[k][L])
                    else:
                        Out.append(-dSg[k][L])
                    # EXF
                    if np.abs(EXF[k][L])>treshold:
                        flows.append(EXF[k][L]/ff)
                        labels.append('$Exf_g$')
                        orientations.append(-1)
                        pathlengths.append(pl)
                    Out.append(-EXF[k][L])
                    # Eg
                    if cMF.wel_yn == 1:
                        if np.abs(Eg[k][L])>treshold:
                            flows.append(-Eg[k][L]/ff)
                            labels.append('$E_g$')
                            orientations.append(1)
                            pathlengths.append(pl)
                        Out.append(Eg[k][L])
                    # Tg   
                    if cMF.wel_yn == 1:
                        if np.abs(Tg[k][L])>treshold:
                            flows.append(-Tg[k][L]/ff)
                            labels.append('$T_g$')
                            orientations.append(1)
                            pathlengths.append(pl)
                        Out.append(Tg[k][L])
                    # DRN
                    if cMF.drn_yn == 1:
                        if np.abs(DRN[k][L])>treshold:
                            flows.append(DRN[k][L]/ff)
                            labels.append('$DRN$')
                            orientations.append(1)
                            pathlengths.append(pl)
                        Out.append(-DRN[k][L])
                    # GHB
                    if cMF.ghb_yn == 1:
                        if np.abs(GHB[k][L])>treshold:
                            flows.append(GHB[k][L]/ff)
                            labels.append('$GHB$')
                            orientations.append(1)
                            pathlengths.append(pl)
                        if  GHB[k][L] < 0.0:
                            Out.append(-GHB[k][L])
                        else:
                            In.append(GHB[k][L])
                    pltsankey.add(patchlabel = '$\Delta S_g$\n%.1f' % (-dSg[k][L]/ff), label='MFL%d'%(L+1), facecolor='LightSteelBlue', trunklength = tl, flows = flows, labels = labels, orientations = orientations, pathlengths = pathlengths, prior=3+L_act, connect=(1, 0))
                    MB_MF.append(100*(sum(In) - sum(Out))/((sum(In) + sum(Out))/2))
                    L_act == 0
                else:
                    L_act -= 1
                    MB_MF.append(np.nan)
            # plot all patches
            diagrams = pltsankey.finish()
            diagrams[-1].patch.set_hatch('/')
            for d in diagrams:
                d.patch.set_edgecolor('gray')
                for t in d.texts:
                    t.set_fontsize(6)
                    t.set_color('navy')
                    t.set_fontweight('bold')
                d.text.set_fontsize(6)
                d.text.set_color('navy')
                d.text.set_fontweight('bold')
            # legend
            if p ==0:
                lblspc = 0.05
                mkscl = 0.1
                handletextpad = 0.25
                borderaxespad = 0.5
                colspc = 0.05
                handlelength = 1.0
                plt.legend(loc='best', labelspacing=lblspc, markerscale=mkscl, handlelength = handlelength, handletextpad = handletextpad, borderaxespad = borderaxespad, ncol = 2, columnspacing = colspc) #loc='upper right'
                leg = plt.gca().get_legend()
                ltext  = leg.get_texts()
                plt.setp(ltext, fontsize=6)
            # water balance box
            #msg = 'Water balance\nclosure (%%)\n-----------------\nMMsurf =%3.1f\nMMsoil =%3.1f\nMFUZF =%3.1f' % (MB_MMsurf, MB_MMsoil, MB_MFuzf) #
            msg = 'Water balance closure (%%)\nMMsurf=%3.1f, MMsoil=%3.1f, MFUZF=%3.1f' % (MB_MMsurf, MB_MMsoil, MB_MFuzf) #
            for L in range(cMF.nlay):
                msg += ' // MFL%d=%3.1f' % (L+1, MB_MF[L]) #%2.1f
            xy = (figtitle.get_position()[0], figtitle.get_position()[1]-0.005)   
            plt.annotate(msg, xy, horizontalalignment='center', verticalalignment='top', fontsize = 6, xycoords='figure fraction')
            if k == 0:
                #print 'Plot of the whole modelled period done!'
                msg = '_0whole'
            else:
                #print 'Plot of hydrological year %d/%d done!' % (year_lst[k-1], year_lst[k-1] + 1)
                msg = '_%d_%d' % (year_lst[k-1], year_lst[k-1] + 1)
            # export png
            if k%2 == 0 or k == len(RF)-1:
                plt.subplots_adjust(left=0.05, bottom=0.10, right=0.95, top=0.95, wspace=0.01, hspace=0.1)
                if f == 0:
                    plt_export_fn = os.path.join(path, '_%s_WBsanley%s.png' % (fntitle, msg))
                else:
                    plt_export_fn = os.path.join(path, '_%s_WBsanley%s_pc.png' % (fntitle,msg))
                plt.savefig(plt_export_fn,dpi=150)
    #print '-------'

##################

def plotCALIBCRIT(calibcritSM, calibcritSMobslst, calibcritHEADS, calibcritHEADSobslst, calibcritHEADSc, calibcritHEADScobslst, plt_export_fn, plt_title, calibcrit, calibcritSMmax = None, calibcritHEADSmax = None, ymin = None, units = '', hnoflo = -9999.9):
    # plot RMSE
    fig = plt.figure()
    fig.suptitle(plt_title, fontsize=10)

    if calibcritSM <> []:
        ax1=fig.add_subplot(2,1,1)
        plt.setp(ax1.get_xticklabels(), fontsize=8)
        plt.setp(ax1.get_yticklabels(), fontsize=8)
        ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        plt.ylabel('%s soil moisture %s' % (calibcrit, units[1]), fontsize=10, horizontalalignment = 'center')
        plt.grid(True)
        xserie = []
        yserie = []
        labels = []
        if calibcritSMmax == None:
            tmp = np.max(list(itertools.chain.from_iterable(calibcritSM)))
            if tmp > 0:
                max_tmp = 1.2*tmp
            else:
                max_tmp = 1.0
            del tmp
        else:
            max_tmp = calibcritSMmax
        n = 0
        for e in calibcritSM:
            if len(e) > 1:
                numtick = len(e)
                if n == 0:
                    newx = 1.0
                else:
                    newx = 1.0 + xserie[-1]
                ee = 0
                lst = range(1,numtick/2+1)
                lst.reverse()
                for i in lst:
                    xserie.append(newx-i/10.0)
                    yserie.append(e[ee])
                    ee += 1
                    labels.append('')
                xserie.append(newx)
                labels.append('%s' % calibcritSMobslst[n])
                if numtick %2 <> 0:
                    yserie.append(e[ee])
                    ee += 1
                else:
                    yserie.append(max_tmp*2.0)
                lst = range(1,numtick/2+1)
                for i in lst:
                    xserie.append(newx+i/10.0)
                    yserie.append(e[ee])
                    ee += 1
                    labels.append('')
            else:
                if n == 0:
                    xserie.append(1.0)
                    yserie.append(e[0])
                    labels.append('%s' % calibcritSMobslst[n])
                else:
                    xserie.append(1+xserie[-1])
                    yserie.append(e[0])
                    labels.append('%s' % calibcritSMobslst[n])
            n += 1
        offset = (max_tmp)*0.05
        for i in range(len(xserie)):
            plt.scatter(xserie[i], yserie[i], marker='o', c = 'orange', s = 15)
            if yserie[i] < max_tmp:
                plt.text(xserie[i], yserie[i]+offset, '%.1f' % yserie[i], fontsize=6, ha = 'center', va = 'center')
        plt.xticks(xserie, labels)
        if ymin == None:
            tmp = np.min(np.ma.masked_where(np.asarray(yserie).flatten() == hnoflo,np.asarray(yserie).flatten()))
            if tmp > 0:
                ymin_tmp = 0
            else:
                ymin_tmp = 1.2*tmp
            del tmp
        else:
            ymin_tmp = ymin
        plt.ylim(ymin_tmp, max_tmp)
        ax1.set_xlim(0, int(max(xserie))+1.0)

    if calibcritHEADS <> []:
        ax2=fig.add_subplot(2,1,2)
        plt.setp(ax2.get_xticklabels(), fontsize=8)
        plt.setp(ax2.get_yticklabels(), fontsize=8)
        ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        plt.ylabel('%s hydraulic heads %s' % (calibcrit, units[0]), fontsize=10, horizontalalignment = 'center')
        plt.grid(True)
        xserie = range(1,len(calibcritHEADSobslst)+1)
        yserie_txt = list(itertools.chain.from_iterable(calibcritHEADS))
        yseriec_txt = list(itertools.chain.from_iterable(calibcritHEADSc))
        if calibcritHEADSmax == None:
            tmp = np.max(list(itertools.chain.from_iterable(calibcritHEADS)))
            if tmp > 0:
                max_tmp = 1.2*tmp
            else:
                max_tmp = 1.0
            del tmp
        else:
            max_tmp = calibcritHEADSmax
#        offset = (max_tmp)*0.05
        for i in range(len(xserie)):
            plt.scatter(xserie[i], calibcritHEADS[i], marker='o', c = 'orange', s = 15)
            for ii, cc in enumerate(calibcritHEADScobslst):
                if calibcritHEADSobslst[i] == cc:
                    plt.scatter(xserie[i], calibcritHEADSc[ii], marker='o', c = 'green', s = 10)
                    if yseriec_txt[ii] < max_tmp:
                        plt.text(xserie[i]-.05, yseriec_txt[ii], '%.1f' % yseriec_txt[ii], fontsize=6, ha = 'right', va = 'center')
            del cc, ii
            if yserie_txt[i] < max_tmp:
                plt.text(xserie[i]+.05, yserie_txt[i], '%.1f' % yserie_txt[i], fontsize=6, ha = 'left', va = 'center')
        plt.xticks(xserie, calibcritHEADSobslst)
        if ymin == None:
            tmp = np.min(np.ma.masked_where(np.asarray(calibcritHEADS).flatten() == hnoflo,np.asarray(calibcritHEADS).flatten()))
            if tmp< 0:
                ymin_tmp = 1.2*tmp
            else:
                ymin_tmp = 0.8*tmp
            del tmp
        else:
            ymin_tmp = ymin
        plt.ylim(ymin_tmp, max_tmp)
        ax2.set_xlim(0, len(calibcritHEADS)+1)

    plt.savefig(plt_export_fn, dpi = 150)

##################
if __name__ == "__main__":
    print '\nWARNING!\nStart MARMITES-MODFLOW models using the script startMARMITES_v3.py\n'

# EOF