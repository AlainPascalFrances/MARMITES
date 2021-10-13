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

def plotTIMESERIES(cMF, i, j, flx, flxLbl, flxIndex_lst, Sm, Sr, plt_export_fn, plt_suptitle, plt_title, clr_lst, hmax,
                   hmin, obs_name, l_obs, nsl, iniMonthHydroYear, date_ini, date_end):
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

    # fmt_DH = mpl.dates.DateFormatter('%Y-%m-%d %H:%M')
    #    cUTIL = MMutils.clsUTILITIES(fmt = fmt_DH)
    #    date_ini, year_ini = cUTIL.compDATE_INI(cMF.inputDate[0], iniMonthHydroYear)
    #    date_end, year_end = cUTIL.compDATE_END(cMF.inputDate[-1], iniMonthHydroYear)
    global l
    date_ini -= 15
    date_end += 15

    dateFmt = mpl.dates.DateFormatter('%Y-%b-%d')
    dateminorFmt = mpl.dates.DateFormatter('%b')
    bymonth = []
    month_tmp = 3
    while len(bymonth) < 3:
        if (iniMonthHydroYear + month_tmp) < 13:
            bymonth.append(iniMonthHydroYear + month_tmp)
        else:
            bymonth.append(iniMonthHydroYear + month_tmp - 12)
        month_tmp += 3
    del month_tmp

    lblspc = 0.05
    mkscale = 1.0
    bdpd = 0.1
    hdltxtpd = 0.05
    colspc = 0.1

    fig = plt.figure(num=None, figsize=(8.27, 11.7), dpi=150)  # (8.5,15), dpi=30)
    fig.suptitle(plt_suptitle)
    fig.text(x=0.5, y=0.05, s=plt_title, horizontalalignment='center', verticalalignment='bottom', fontsize=9)

    ax1 = fig.add_subplot(8, 1, 1)
    ax1.bar(cMF.inputDate, flx[flxIndex_lst[b'iRF']], color='darkblue', linewidth=0, align='center',
            label=flxLbl[flxIndex_lst[b'iRF']])
    ax1.bar(cMF.inputDate, flx[flxIndex_lst[b'iRFe']], color='deepskyblue', linewidth=0, align='center',
            label=flxLbl[flxIndex_lst[b'iRFe']])
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # y
    plt.ylabel('mm', fontsize=10)
    ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax1.get_yticklabels(), fontsize=8)
    # x
    ax1.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax1.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax1.get_xticklabels(which='both'), visible=False)
    # grd
    ax1.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax1.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax2 = fig.add_subplot(8, 1, 2, sharex=ax1)
    if np.sum(np.abs(flx[flxIndex_lst[b'iRo']])) > 1E-7:
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iRo']], '-', c='blue', linewidth=1.0,
                      label=flxLbl[flxIndex_lst[b'iRo']])
    colors_nsl = itertools.cycle(clr_lst)
    if np.sum(np.abs(flx[flxIndex_lst[b'iInf']])) > 1E-7:
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iInf']], '--', c=next(colors_nsl), linewidth=1.5,
                      label=flxLbl[flxIndex_lst[b'iInf']])
    if np.sum(np.abs(flx[flxIndex_lst[b'iEsurf']])) > 1E-7:
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iEsurf']], '-', c='deepskyblue', linewidth=0.75,
                      label=flxLbl[flxIndex_lst[b'iEsurf']])
    if np.sum(np.abs(flx[flxIndex_lst[b'iSsurf']])) > 1E-7:
        plt.bar(cMF.inputDate, flx[flxIndex_lst[b'iSsurf']], color='lightblue', linewidth=0, align='center',
                label=flxLbl[flxIndex_lst[b'iSsurf']])
    if np.sum(np.abs(flx[flxIndex_lst[b'idSsurf']])) > 1E-7:
        plt.bar(cMF.inputDate, flx[flxIndex_lst[b'idSsurf']], color='darkblue', linewidth=0, align='center',
                label=flxLbl[flxIndex_lst[b'idSsurf']])
    try:
        if b'iRoobs' in flxIndex_lst:
            Roobs_m = np.ma.masked_values(flx[flxIndex_lst[b'iRoobs']], cMF.hnoflo, atol=0.09)
            #        dgwtobsmin = np.ma.min(dgwtobs_m)
            plt.plot_date(cMF.inputDate, Roobs_m,  'o', ls='None', color='lightblue',markeredgecolor='lightblue',
                          markerfacecolor='None', markersize=2,
                          label=flxLbl[flxIndex_lst[b'iRoobs']])  # ls='--', color = 'blue'
    except:
        print("ERROR plotting Ro obs")
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc, numpoints=3)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # grd
    ax2.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax2.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')
    # x
    ax2.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax2.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax2.get_xticklabels(which='both'), visible=False)
    # y
    plt.ylabel('mm', fontsize=10)
    ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax2.get_yticklabels(), fontsize=8)

    colors_nsl = itertools.cycle(clr_lst)
    ax3 = fig.add_subplot(8, 1, 3, sharex=ax1)
    # PE
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iPE']], '-', color='lightblue', linewidth=3,
                  label=flxLbl[flxIndex_lst[b'iPE']])
    if cMF.wel_yn == 1:
        E = flx[flxIndex_lst[b'iEsoil']] + flx[flxIndex_lst[b'iEg']]
        # Etot
        plt.plot_date(cMF.inputDate, E, '-', color='darkblue', linewidth=1.5, label=r'$E$')
    # Esoil
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iEsoil']], '--', color='brown',
                  label=flxLbl[flxIndex_lst[b'iEsoil']])
    if cMF.wel_yn == 1:
        # Eg
        if np.absolute(sum(flx[flxIndex_lst[b'iEg']])) > 1E-6:
            plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iEg']], '-', color='blue',
                          label=flxLbl[flxIndex_lst[b'iEg']])
    for l in range(nsl):
        if np.absolute(sum(flx[flxIndex_lst[b'iEsoil_%d' % (l + 1)]])):
            ax3.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iEsoil_%d' % (l + 1)]], '-', color=next(colors_nsl),
                          label=flxLbl[flxIndex_lst[b'iEsoil_%d' % (l + 1)]])
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=3,
               columnspacing=colspc)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # grd
    ax3.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax3.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')
    # x
    ax3.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax3.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax3.get_xticklabels(which='both'), visible=False)
    # y
    plt.ylabel('mm', fontsize=10)
    ax3.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax3.get_yticklabels(), fontsize=8)

    colors_nsl = itertools.cycle(clr_lst)
    ax4 = fig.add_subplot(8, 1, 4, sharex=ax1)
    # PT
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iPT']], '-', color='lightblue', linewidth=3,
                  label=flxLbl[flxIndex_lst[b'iPT']])
    if cMF.wel_yn == 1:
        T = flx[flxIndex_lst[b'iTsoil']] + flx[flxIndex_lst[b'iTg']]
        # Ttot
        plt.plot_date(cMF.inputDate, T, '-', color='darkblue', linewidth=1.5, label=r'$T$')
    # Tsoil
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iTsoil']], '--', color='brown',
                  label=flxLbl[flxIndex_lst[b'iTsoil']])
    if cMF.wel_yn == 1:
        # Tg
        if np.absolute(sum(flx[flxIndex_lst[b'iTg']])) > 1E-6:
            plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iTg']], '-', color='blue',
                          label=flxLbl[flxIndex_lst[b'iTg']])
    for l in range(nsl):
        if np.absolute(sum(flx[flxIndex_lst[b'iTsoil_%d' % (l + 1)]])):
            ax4.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iTsoil_%d' % (l + 1)]], '-', color=next(colors_nsl),
                          label=flxLbl[flxIndex_lst[b'iTsoil_%d' % (l + 1)]])
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=3,
               columnspacing=colspc)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # grd
    ax4.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax4.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')
    # x
    ax4.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax4.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax4.get_xticklabels(which='both'), visible=False)
    # y
    plt.ylabel('mm', fontsize=10)
    ax4.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax4.get_yticklabels(), fontsize=8)

    colors_nsl = itertools.cycle(clr_lst)
    ax6 = fig.add_subplot(8, 1, 5, sharex=ax1)
    if np.sum(np.abs(flx[flxIndex_lst[b'iEXFg']])) > 1E-7:
        plt.bar(cMF.inputDate, flx[flxIndex_lst[b'iEXFg']], color='lightblue', linewidth=0, align='center',
                label=flxLbl[flxIndex_lst[b'iEXFg']])
    if cMF.wel_yn == 1:
        if np.sum(np.abs(flx[flxIndex_lst[b'iETg']])) > 1E-7:
            plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iETg']], '-', c='blue', linewidth=1.5,
                          label=flxLbl[flxIndex_lst[b'iETg']])
    for l in range(nsl):
        if np.sum(np.abs(flx[flxIndex_lst[b'iRsoil_%d' % (l + 1)]])) > 1E-7:
            ax6.plot_date(cMF.inputDate, -1.0 * flx[flxIndex_lst[b'iRsoil_%d' % (l + 1)]], '-', color=next(colors_nsl),
                          label=flxLbl[flxIndex_lst[b'iRsoil_%d' % (l + 1)]])
    if np.sum(np.abs(flx[flxIndex_lst[b'iRg']])) > 1E-7:
        plt.plot_date(cMF.inputDate, -1.0 * flx[flxIndex_lst[b'iRg']], '-', c='darkblue', linewidth=1.5,
                      label=flxLbl[flxIndex_lst[b'iRg']])
    colors_nsl = itertools.cycle(clr_lst)
    for l in range(nsl):
        if np.sum(np.abs(flx[flxIndex_lst[b'iExf_%d' % (l + 1)]])) > 1E-7:
            ax6.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iExf_%d' % (l + 1)]], '--', color=next(colors_nsl),
                          label=flxLbl[flxIndex_lst[b'iExf_%d' % (l + 1)]])
            # y
    ax6.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.ylabel('mm', fontsize=10)
    plt.setp(ax6.get_yticklabels(), fontsize=8)
    # leg
    plt.legend(labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=3,
               columnspacing=colspc, loc=0)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # grd
    ax6.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax6.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')
    # x
    ax6.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax6.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax6.get_xticklabels(which='both'), visible=False)

    colors_nsl = itertools.cycle(clr_lst)
    ax5 = fig.add_subplot(8, 1, 6, sharex=ax1)
    for l in range(nsl):
        try:
            if b'iSobs_%d' % (l + 1) in flxIndex_lst:
                if list(flx[flxIndex_lst[b'iSobs_%d' % (l + 1)]]):
                    ax5.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iSobs_%d' % (l + 1)]], 'o', ls='None', color='gray',
                                  markersize=2, markeredgecolor=next(colors_nsl), markerfacecolor='None',
                                  label=flxLbl[flxIndex_lst[b'iSobs_%d' % (l + 1)]])  # '--', color = color,  markevery = 2
        except:
            print("ERROR plotting SM obs")
    sim_tmp = []
    colors_nsl = itertools.cycle(clr_lst)
    for l in range(nsl):
        y = flx[flxIndex_lst[b'iSsoil_pc_s_%d' % (l + 1)]]
        sim_tmp.append(y)
        y = np.ma.masked_where(y < 0.0, y)
        ax5.plot_date(cMF.inputDate, y, '-', color=next(colors_nsl),
                      label=flxLbl[flxIndex_lst[b'iSsoil_pc_s_%d' % (l + 1)]])
        del y
    # y axis
    ybuffer = 0.1 * (max(Sm) - min(Sr))
    plt.ylim(min(Sr) - ybuffer, max(Sm) + ybuffer)
    ax5.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.ylabel('%', fontsize=10)
    plt.setp(ax5.get_yticklabels(), fontsize=8)
    # legend
    # lbl_Spcobs = lbl_Sobs + lbl_S
    # plt.legend(lbl_Spcobs, loc=0, labelspacing=lblspc, markerscale=mkscale)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc, numpoints=3)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # grd
    ax5.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax5.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')
    # x
    ax5.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax5.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax5.get_xticklabels(which='both'), visible=False)

    ax7 = fig.add_subplot(8, 1, 7, sharex=ax1)
    lines = itertools.cycle(
        ['-', '--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|',
         '_'])
    dgwtMFmax = []
    for L in range(cMF.nlay):
        dgwtMF = flx[flxIndex_lst[b'id_%d' % (L + 1)]]
        dgwtMFmax.append(np.max(dgwtMF))
        plt.plot_date(cMF.inputDate, dgwtMF, next(lines), color='b', markersize=2, markevery=7,
                      label=flxLbl[flxIndex_lst[b'id_%d' % (L + 1)]])
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'idcorr']], '--', c='g', markersize=2, markevery=7,
                  label=flxLbl[flxIndex_lst[b'idcorr']])
    try:
        if b'idobs' in flxIndex_lst:
            plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'idobs']], 'o', ls='None', color='LightBlue',
                          markeredgecolor='LightBlue', markerfacecolor='None', markersize=2,
                          label=flxLbl[flxIndex_lst[b'idobs']])  # ls='--', color = 'blue'  markevery = 7,
    except:
        print("ERROR plotting dgw obs")
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc, numpoints=3)
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
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # y axis
    ax7.set_ylabel('m', fontsize=10)
    ax7.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    ax7.set_ylim(ymax=ymax)
    plt.setp(ax7.get_yticklabels(), fontsize=8)
    # x axis
    ax7.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax7.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    ax7.xaxis.set_major_formatter(dateFmt)
    ax7.xaxis.set_minor_formatter(dateminorFmt)
    labels = ax7.get_xticklabels(which='both')
    plt.setp(labels, rotation=90, fontsize=8)
    del labels
    ax7.set_xlabel('Date', fontsize=10)
    # grd
    ax7.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax7.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax1.set_xlim(date_ini, date_end)

    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig(plt_export_fn, dpi=150)
    plt.close()

    # --------------FIGURE PART2----------

    fig = plt.figure(num=None, figsize=(8.27, 11.7), dpi=150)  # (8.5,15), dpi=30)
    fig.suptitle(plt_suptitle + ' - part 2')
    fig.text(x=0.5, y=0.05, s=plt_title, horizontalalignment='center', verticalalignment='bottom', fontsize=9)

    ax1 = fig.add_subplot(8, 1, 4)
    lines = itertools.cycle(
        ['-', '--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|',
         '_'])
    for L in range(cMF.nlay):
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'ih_%d' % (L + 1)]], next(lines), color='b',
                      markersize=2, markevery=7, label=flxLbl[flxIndex_lst[b'ih_%d' % (L + 1)]])
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'ihcorr']], '--', color='g', markersize=2, markevery=7,
                  label=flxLbl[flxIndex_lst[b'ihcorr']])
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'ih_SF']], '-', color='r', markersize=2, markevery=7,
                  label=flxLbl[flxIndex_lst[b'ih_SF']])
    try:
        if b'ihobs' in flxIndex_lst:
            plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'ihobs']], 'o', ls='None', color='LightBlue',
                          markeredgecolor='LightBlue', markerfacecolor='None', markersize=2,
                          label=flxLbl[flxIndex_lst[b'ihobs']])  # ls='--', color = 'blue' markevery = 7,
    except:
        print("ERROR plotting h obs")
    # y
    ybuffer = 0.1 * (hmax - hmin)
    if ybuffer == 0.0:
        ybuffer = 1.0
    plt.ylim((hmin - ybuffer, hmax + ybuffer))
    plt.ylabel('m', fontsize=10)
    ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax1.get_yticklabels(), fontsize=8)
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc, numpoints=3)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # x
    ax1.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax1.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax1.get_xticklabels(which='both'), visible=False)
    # grd
    ax1.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax1.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    colors_nsl = itertools.cycle(clr_lst)
    ax8b = fig.add_subplot(8, 1, 1, sharex=ax1)
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iMBsurf']], '-', c='lightblue',
                  label=flxLbl[flxIndex_lst[b'iMBsurf']])
    MBmin = [min(flx[flxIndex_lst[b'iMBsurf']])]
    MBmax = [max(flx[flxIndex_lst[b'iMBsurf']])]
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iMB']], '-', c='r', label=flxLbl[flxIndex_lst[b'iMB']])
    MBmin.append(min(flx[flxIndex_lst[b'iMB']]))
    MBmax.append(max(flx[flxIndex_lst[b'iMB']]))
    for l in range(nsl):
        ax8b.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iMB_s_%d' % (l + 1)]], '-', color=next(colors_nsl),
                       label=flxLbl[flxIndex_lst[b'iMB_s_%d' % (l + 1)]])
        MBmin.append(min(flx[flxIndex_lst[b'iMB_s_%d' % (l + 1)]]))
        MBmax.append(max(flx[flxIndex_lst[b'iMB_s_%d' % (l + 1)]]))
    # y axis
    plt.ylabel('mm', fontsize=10)
    MBmax = max(MBmax)
    MBmin = min(MBmin)
    if np.abs(MBmax - MBmin) < 1E-7:
        plt.ylim(-0.1, 0.1)
    else:
        minfact = 0.95
        maxfact = 1.05
        if MBmin > 0:
            minfact = 1.05
        if MBmax < 0:
            maxfact = 0.95
        plt.ylim(MBmin * minfact, MBmax * maxfact)
    plt.setp(ax8b.get_yticklabels(), fontsize=8)
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    ax8b.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    # x
    ax8b.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax8b.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax8b.get_xticklabels(which='both'), visible=False)
    # grd
    ax8b.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax8b.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax20 = fig.add_subplot(8, 1, 2, sharex=ax1)
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iRo']], '-', c='blue', linewidth=1.0,
                  label=flxLbl[flxIndex_lst[b'iRo']])
    ymax = None
    try:
        #        dgwtobsmin = np.ma.min(dgwtobs_m)
        if b'iRoobs' in flxIndex_lst:
            plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iRoobs']], 'o', ls='None', color='lightblue',
                          markeredgecolor='lightblue', markerfacecolor='None', markersize=2,
                          label=flxLbl[flxIndex_lst[b'iRoobs']])  # ls='--', color = 'blue'
            ymax = np.ma.max(flx[flxIndex_lst[b'iRoobs']])
    except:
        print("ERROR plotting Ro obs")
    # y
    plt.ylabel('mm', fontsize=10)
    plt.setp(ax20.get_yticklabels(), fontsize=8)
    if ymax is not None:
        plt.ylim(ymax=ymax)
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc, numpoints=3)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # x
    ax20.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax20.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax20.get_xticklabels(which='both'), visible=False)
    # grd
    ax20.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax20.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    colors_nsl = itertools.cycle(clr_lst)
    ax9a = fig.add_subplot(16, 1, 5, sharex=ax1)
    for l in range(nsl):
        ax9a.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iSAT_%d' % (l + 1)]], '-', color=next(colors_nsl),
                       label=flxLbl[flxIndex_lst[b'iSAT_%d' % (l + 1)]])
    # y axis
    plt.ylim(-0.1, 1.1)
    plt.ylabel('SAT', fontsize=10)
    ax9a.yaxis.set_ticks(np.arange(0, 1.25, 1))
    ax9a.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%1d'))
    plt.setp(ax9a.get_yticklabels(), fontsize=8)
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # x
    ax9a.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax9a.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax9a.get_xticklabels(which='both'), visible=False)
    # grd
    ax9a.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax9a.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax10a = fig.add_subplot(16, 1, 6, sharex=ax1)
    uzthick = flx[flxIndex_lst[b'iuzthick']]
    plt.plot_date(cMF.inputDate, -uzthick, '-', c='brown', label=flxLbl[flxIndex_lst[b'iuzthick']])
    # y axis
    plt.ylabel('m', fontsize=10)
    minfact = 0.95
    maxfact = 1.05
    if np.min(uzthick) < 0:
        minfact = 1.05
    if np.max(uzthick) < 0:
        maxfact = 0.95
    if np.min(uzthick) == np.max(uzthick) == 0.0:
        plt.ylim(-0.5, 0.5)
    else:
        plt.ylim(np.min(-uzthick) * minfact, np.max(-uzthick) * maxfact)
    plt.setp(ax10a.get_yticklabels(), fontsize=8)
    ax10a.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # x
    ax10a.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax10a.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax10a.get_xticklabels(which='both'), visible=False)
    # grd
    ax10a.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax10a.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax7 = fig.add_subplot(8, 1, 5, sharex=ax1)
    lines = itertools.cycle(
        ['-', '--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|',
         '_'])
    dgwtMFmax = []
    for L in range(cMF.nlay):
        dgwtMFmax.append(np.max(flx[flxIndex_lst[b'id_%d' % (L + 1)]]))
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'id_%d' % (L + 1)]], next(lines), color='b', markersize=2,
                      markevery=7,
                      label=flxLbl[flxIndex_lst[b'id_%d' % (L + 1)]])
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'idcorr']], '--', c='g', markersize=2, markevery=7,
                  label=flxLbl[flxIndex_lst[b'idcorr']])
    try:
        if b'idobs' in flxIndex_lst:
            if flx[flxIndex_lst[b'idobs']].size > 0:   # array.size>0
                plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'idobs']], 'o', ls='None', color='LightBlue',
                              markeredgecolor='LightBlue', markerfacecolor='None', markersize=2,
                              label=flxLbl[flxIndex_lst[b'idobs']])  # ls='--', color = 'blue' markevery = 7,
    except:
        print("ERROR plotting dgw obs")
    # y axis
    plt.ylabel('m', fontsize=10)
    plt.setp(ax7.get_yticklabels(), fontsize=8)
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc, numpoints=3)
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
    plt.ylim(ymax=ymax)
    ax7.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # x
    ax7.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax7.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax7.get_xticklabels(which='both'), visible=False)
    # grd
    ax7.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax7.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    colors_nsl = itertools.cycle(clr_lst)
    ax10b = fig.add_subplot(8, 1, 6, sharex=ax1)
    for l in range(nsl):
        y = flx[flxIndex_lst[b'iSsoil_%d' % (l + 1)]]
        y = np.ma.masked_where(y < 0.0, y)
        ax10b.plot_date(cMF.inputDate, y, '-', color=next(colors_nsl),
                        label=flxLbl[flxIndex_lst[b'iSsoil_%d' % (l + 1)]])
    # y axis
    plt.ylim(0, np.max(flx[flxIndex_lst[b'iSsoil_%d' % (l + 1)]]) * 1.05)
    plt.ylabel('mm', fontsize=10)
    ax10b.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax10b.get_yticklabels(), fontsize=8)
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # x
    ax10b.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax10b.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax10b.get_xticklabels(which='both'), visible=False)
    # grd
    ax10b.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax10b.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    colors_nsl = itertools.cycle(clr_lst)
    ax5 = fig.add_subplot(8, 1, 7, sharex=ax1)
    colors_nsl = itertools.cycle(clr_lst)
    for l in range(nsl):
        y = flx[flxIndex_lst[b'iSsoil_pc_s_%d' % (l + 1)]]
        y = np.ma.masked_where(y < 0.0, y)
        ax5.plot_date(cMF.inputDate, y, '-', color=next(colors_nsl),
                      label=flxLbl[flxIndex_lst[b'iSsoil_pc_s_%d' % (l + 1)]])
    for l in range(nsl):
        try:
            if b'iSobs_%d' % (l + 1) in flxIndex_lst:
                if list(flx[flxIndex_lst[b'iSobs_%d' % (l + 1)]]):
                    y = flx[flxIndex_lst[b'iSobs_%d' % (l + 1)]]
                    ax5.plot_date(cMF.inputDate, y, 'o', ls='None', color='gray', markersize=2,
                                  markeredgecolor=next(colors_nsl), markerfacecolor='None',
                                  label=flxLbl[flxIndex_lst[b'iSobs_%d' % (l + 1)]])  # '--', color = color, markevery = 2,
        except:
            print("ERROR plotting SM obs")
    # y axis
    ybuffer = 0.1 * (max(Sm) - min(Sr))
    plt.ylim(min(Sr) - ybuffer, max(Sm) + ybuffer)
    plt.ylabel('%', fontsize=10)
    ax5.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax5.get_yticklabels(), fontsize=8)
    # legend
    # lbl_Spcobs = lbl_Sobs + lbl_S
    # plt.legend(lbl_Spcobs, loc=0, labelspacing=lblspc, markerscale=mkscale)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc, numpoints=3)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # x axis
    ax5.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax5.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    ax5.xaxis.set_major_formatter(dateFmt)
    ax5.xaxis.set_minor_formatter(dateminorFmt)
    labels = ax5.get_xticklabels(which='both')
    plt.setp(labels, rotation=90, fontsize=8)
    del labels
    ax5.set_xlabel('Date', fontsize=10)
    # grd
    ax5.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax5.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax1.set_xlim(np.min(cMF.inputDate) - 15.0, np.max(cMF.inputDate) + 15)

    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    txt = plt_export_fn.split('.')
    plt_export_fn = txt[0] + '_part2.' + txt[1]
    plt.savefig(plt_export_fn, dpi=150)
    #    plt.show()
    plt.clf()
    plt.close()


##################

def plotTIMESERIES_flxGW(cMF, flx, flxLbl, flxIndex_lst, plt_export_fn, plt_title, iniMonthHydroYear, date_ini,
                         date_end):
    """
    Plot the time series of the fluxes observed from the whole catchment
    Use Matplotlib
    """
    dateFmt = mpl.dates.DateFormatter('%Y-%b-%d')
    dateminorFmt = mpl.dates.DateFormatter('%b')
    bymonth = []
    month_tmp = 3
    while len(bymonth) < 3:
        if (iniMonthHydroYear + month_tmp) < 13:
            bymonth.append(iniMonthHydroYear + month_tmp)
        else:
            bymonth.append(iniMonthHydroYear + month_tmp - 12)
        month_tmp += 3
    del month_tmp

    lblspc = 0.05
    mkscale = 1.0
    bdpd = 0.1
    hdltxtpd = 0.05
    colspc = 0.1

    date_ini -= 15
    date_end += 15

    fig = plt.figure(num=None, figsize=(8.27, 11.7), dpi=150)  # (8.5,15), dpi=30)
    fig.suptitle(plt_title + ' - part 3')

    ax0 = fig.add_subplot(8, 1, 1)
    # RF
    ax0.bar(cMF.inputDate, flx[flxIndex_lst[b'iRF']], color='darkblue', linewidth=0, align='center',
            label=flxLbl[flxIndex_lst[b'iRF']])
    # RFe
    ax0.bar(cMF.inputDate, flx[flxIndex_lst[b'iRFe']], color='deepskyblue', linewidth=0, align='center',
            label=flxLbl[flxIndex_lst[b'iRFe']])
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # y
    plt.ylabel('mm', fontsize=10)
    ax0.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax0.get_yticklabels(), fontsize=8)
    # x
    ax0.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax0.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax0.get_xticklabels(which='both'), visible=False)
    # grd
    ax0.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax0.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    # plot Ro
    ax2 = fig.add_subplot(8, 1, 2, sharex=ax0)
    # Ro
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iRo']], '-', c='blue', linewidth=1.0,
                  label=flxLbl[flxIndex_lst[b'iRo']])
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc, numpoints=3)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # y
    plt.ylabel('mm', fontsize=10)
    ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    # x axis
    ax2.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax2.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    ax2.xaxis.set_major_formatter(dateFmt)
    ax2.xaxis.set_minor_formatter(dateminorFmt)
    labels = ax2.get_xticklabels(which='both')
    plt.setp(labels, rotation=90, fontsize=8)
    del labels
    ax2.set_xlabel('Date', fontsize=10)
    # grd
    ax2.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax2.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    # plot gwd
    lines = itertools.cycle(
        ['-', '--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|',
         '_'])
    ax1 = fig.add_subplot(8, 1, 4, sharex=ax0)
    for L in range(cMF.nlay):
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'id_%d' % (L + 1)]], next(lines), color='b', markersize=2,
                      markevery=7,
                      label=flxLbl[flxIndex_lst[b'id_%d' % (L + 1)]])
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc, numpoints=3)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # y
    ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.ylabel('m', fontsize=10)
    plt.setp(ax1.get_yticklabels(), fontsize=8)
    # x
    ax1.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax1.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax1.get_xticklabels(which='both'), visible=False)
    # grd
    ax1.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax1.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    # plot GW fluxes
    ax8 = fig.add_subplot(2, 1, 2, sharex=ax0)
    # uzf recharge
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iRg']], '-', c='darkblue', linewidth=1.5,
                  label=flxLbl[flxIndex_lst[b'iRg']])
    lines = itertools.cycle(
        ['-', '--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|',
         '_'])
    for i in range(flxIndex_lst[b'idSg_1'], len(flxIndex_lst)):
        if np.absolute(sum(flx[i])) > 1E-6:
            plt.plot_date(cMF.inputDate, flx[i], next(lines), color=mpl.colors.rgb2hex(np.random.rand(1, 3)[0]),
                          markersize=2, markevery=7, label=flxLbl[i], markeredgecolor='None')
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=4,
               columnspacing=colspc, numpoints=3)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # y
    plt.ylabel('mm', fontsize=10)
    ax8.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax8.get_yticklabels(), fontsize=8)
    # x axis
    ax8.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax8.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    ax8.xaxis.set_major_formatter(dateFmt)
    ax8.xaxis.set_minor_formatter(dateminorFmt)
    labels = ax8.get_xticklabels(which='both')
    plt.setp(labels, rotation=90, fontsize=8)
    del labels
    ax8.set_xlabel('Date', fontsize=10)
    # grd
    ax8.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax8.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax0.set_xlim(np.min(cMF.inputDate) - 15.0, np.max(cMF.inputDate) + 15)

    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    txt = plt_export_fn.split('.')
    plt_export_fn = txt[0] + '_part3MF.' + txt[1]
    plt.savefig(plt_export_fn, dpi=150)
    #    plt.show()
    plt.clf()
    plt.close('all')


##################

def plotTIMESERIES_CATCH(cMF, flx, flxLbl, plt_export_fn, plt_title, hmax, hmin, iniMonthHydroYear, date_ini, date_end,
                         flxIndex_lst, obs_catch=None, obs_catch_list=[0, 0, 0], TopSoilAverage=None, MF=None):
    """
    Plot the time series of the fluxes observed from the whole catchment
    Use Matplotlib
    """

    # global Roobs_m, obs_h, l, hobs_m
    dateFmt = mpl.dates.DateFormatter('%Y-%b-%d')
    dateminorFmt = mpl.dates.DateFormatter('%b')
    bymonth = []
    month_tmp = 3
    while len(bymonth) < 3:
        if (iniMonthHydroYear + month_tmp) < 13:
            bymonth.append(iniMonthHydroYear + month_tmp)
        else:
            bymonth.append(iniMonthHydroYear + month_tmp - 12)
        month_tmp += 3
    del month_tmp

    lblspc = 0.05
    mkscale = 1.0
    bdpd = 0.1
    hdltxtpd = 0.05
    colspc = 0.1

    # fmt_DH = mpl.dates.DateFormatter('%Y-%m-%d')
    # cUTIL = MMutils.clsUTILITIES(fmt = fmt_DH)

    # date_ini, year_ini = cUTIL.compDATE_INI(cMF.inputDate[0], iniMonthHydroYear)
    # date_end, year_end = cUTIL.compDATE_END(cMF.inputDate[-1], iniMonthHydroYear)

    fig = plt.figure(num=None, figsize=(8.27, 11.7), dpi=150)  # (8.5,15), dpi=30)
    fig.suptitle(plt_title)

    ax1 = fig.add_subplot(8, 1, 1)
    # ax1.set_autoscalex_on(False)
    # RF
    ax1.bar(cMF.inputDate, flx[flxIndex_lst[b'iRF']], color='darkblue', linewidth=0, align='center',
            label=flxLbl[flxIndex_lst[b'iRF']])
    # RFe
    ax1.bar(cMF.inputDate, flx[flxIndex_lst[b'iRFe']], color='deepskyblue', linewidth=0, align='center',
            label=flxLbl[flxIndex_lst[b'iRFe']])
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # y
    plt.ylabel('mm', fontsize=10)
    ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax1.get_yticklabels(), fontsize=8)
    # x
    ax1.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax1.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax1.get_xticklabels(which='both'), visible=False)
    # grd
    ax1.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax1.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax2 = fig.add_subplot(8, 1, 2, sharex=ax1)
    # ax2.set_autoscalex_on(False)
    # Ro
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iRo']], '-', c='blue', linewidth=1.0,
                  label=flxLbl[flxIndex_lst[b'iRo']])  # 'r-', c='blue',
    # Inf
    if np.sum(np.abs(flx[flxIndex_lst[b'iInf']])) > 1E-7:
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iInf']], '--', color='brown',
                      label=flxLbl[flxIndex_lst[b'iInf']])
    # Esurf
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iEsurf']], '-', c='deepskyblue', linewidth=1.0,
                  label=flxLbl[flxIndex_lst[b'iEsurf']])
    # Ssurf
    plt.bar(cMF.inputDate, flx[flxIndex_lst[b'iSsurf']], color='lightblue', linewidth=0, align='center',
            label=flxLbl[flxIndex_lst[b'iSsurf']])
    # DeltaSsurf
    plt.bar(cMF.inputDate, flx[flxIndex_lst[b'idSsurf']], color='darkblue', linewidth=0, align='center',
            label=flxLbl[flxIndex_lst[b'idSsurf']])
    # Ro obs
    try:
        if obs_catch_list[2] == 1:
            obs_Ro = obs_catch.get('catch')['obs_Ro']
            Roobs_m = np.ma.masked_values(obs_Ro[0], cMF.hnoflo, atol=0.09)
            plt.plot_date(cMF.inputDate, Roobs_m, 'o', markerfacecolor='None', markeredgecolor='lightblue',
                          markersize=2, label=r'$Ro \ obs$')  # markevery = 7,
            print('RMSE/RSR/NSE/r of obs. at the catch. scale')
            rmse, rsr, nse, r = cMF.cPROCESS.compCalibCrit(flx[flxIndex_lst[b'iRo']], obs_Ro[0], cMF.hnoflo)
            rmseRo = [rmse]
            rsrRo = [rsr]
            nseRo = [nse]
            rRo = [r]
            del rmse, rsr, nse, r
            if rmseRo[0] is not None:
                print('Ro: %.1f mm / %.2f / %.2f / %.2f' % (rmseRo[0], rsrRo[0], nseRo[0], rRo[0]))
    except:
       print("ERROR plotting Ro obs. catchment")
    # y
    plt.ylabel('mm', fontsize=10)
    ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc, numpoints=3)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # x
    ax2.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax2.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax2.get_xticklabels(which='both'), visible=False)
    # grd
    ax2.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax2.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax3 = fig.add_subplot(8, 1, 3, sharex=ax1)
    # ax3.set_autoscalex_on(False)
    # PE
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iPE']], '-', color='lightblue', linewidth=2.5,
                  label=flxLbl[flxIndex_lst[b'iPE']])
    if cMF.wel_yn == 1:
        E = flx[flxIndex_lst[b'iEsoil']] + flx[flxIndex_lst[b'iEg']]
        # Etot
        plt.plot_date(cMF.inputDate, E, '-', color='darkblue', linewidth=1.5, label=r'$E$')
    # Esoil
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iEsoil']], '--', color='brown',
                  label=flxLbl[flxIndex_lst[b'iEsoil']])
    if cMF.wel_yn == 1:
        # Eg
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iEg']], '-', color='blue', label=flxLbl[flxIndex_lst[b'iEg']])
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # y
    ax3.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.ylabel('mm', fontsize=10)
    plt.setp(ax3.get_yticklabels(), fontsize=8)
    # x
    ax3.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax3.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax3.get_xticklabels(which='both'), visible=False)
    # grd
    ax3.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax3.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax4 = fig.add_subplot(8, 1, 4, sharex=ax1)
    # ax4.set_autoscalex_on(False)
    # PT
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iPT']], '-', color='lightblue', linewidth=3,
                  label=flxLbl[flxIndex_lst[b'iPT']])
    if cMF.wel_yn == 1:
        T = flx[flxIndex_lst[b'iTsoil']] + flx[flxIndex_lst[b'iTg']]
        # Ttot
        plt.plot_date(cMF.inputDate, T, '-', color='darkblue', linewidth=1.5, label=r'$T$')
    # Tsoil
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iTsoil']], '--', color='brown',
                  label=flxLbl[flxIndex_lst[b'iTsoil']])
    if cMF.wel_yn == 1:
        # Tg
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iTg']], '-', color='blue', label=flxLbl[flxIndex_lst[b'iTg']])
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # y
    plt.ylabel('mm', fontsize=10)
    ax4.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax4.get_yticklabels(), fontsize=8)
    # x
    ax4.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax4.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax4.get_xticklabels(which='both'), visible=False)
    # grd
    ax4.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax4.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax5 = fig.add_subplot(8, 1, 5, sharex=ax1)
    # ax5.set_autoscalex_on(False)
    # Rp
    ax5.plot_date(cMF.inputDate, -1.0 * flx[flxIndex_lst[b'iperc']], '-', color='brown',
                  label=flxLbl[flxIndex_lst[b'iRsoil']])
    if cMF is not None:
        # Rg
        plt.plot_date(cMF.inputDate, -1.0 * flx[flxIndex_lst[b'iRg']], '-', c='darkblue', linewidth=1.5,
                      label=flxLbl[flxIndex_lst[b'iRg']])
    if cMF.wel_yn == 1:
        # ETg
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iETg']], '-', c='blue', linewidth=1.5,
                      label=flxLbl[flxIndex_lst[b'iETg']])
    # EXF
    plt.bar(cMF.inputDate, flx[flxIndex_lst[b'iEXFg']], color='lightblue', linewidth=0, align='center',
            label=flxLbl[flxIndex_lst[b'iEXFg']])
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # y
    plt.ylabel('mm', fontsize=10)
    ax5.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax5.get_yticklabels(), fontsize=8)
    # x axis
    ax5.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax5.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    ax5.xaxis.set_major_formatter(dateFmt)
    ax5.xaxis.set_minor_formatter(dateminorFmt)
    labels = ax5.get_xticklabels(which='both')
    plt.setp(labels, rotation=90, fontsize=8)
    del labels
    ax5.set_xlabel('Date', fontsize=10)
    # grd
    ax5.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax5.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax6 = fig.add_subplot(8, 1, 7, sharex=ax1)
    # ax6.set_autoscalex_on(False)
    plt.setp(ax6.get_yticklabels(), fontsize=8)
    # theta
    rmseSM = None
    rsrSM = None
    nseSM = None
    rSM = None
    ax6.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iSsoil_pc']], '-', color='brown',
                  label=flxLbl[flxIndex_lst[b'iSsoil_pc']])
    if obs_catch_list[1] == 1:
        obs_SM = obs_catch.get('catch')['obs_SM']
        Sobs_m = np.ma.masked_values(obs_SM[0], cMF.hnoflo, atol=0.09)
        ax6.plot_date(cMF.inputDate, Sobs_m, 'o', markerfacecolor='None', markeredgecolor='brown', markersize=2,
                      label=r'$\theta \ obs$')  # markevery = 7,
        rmse, rsr, nse, r = cMF.cPROCESS.compCalibCrit(flx[flxIndex_lst[b'iSsoil_pc']], obs_SM[0], cMF.hnoflo)
        rmseSM = [100.0 * rmse]
        rsrSM = [rsr]
        nseSM = [nse]
        rSM = [r]
        del rmse, rsr, nse, r
        if rmseSM[0] is not None:
            print('SM: %.1f %% / %.2f / %.2f / %.2f' % (rmseSM[0], rsrSM[0], nseSM[0], rSM[0]))
    # y axis
    plt.ylabel('%', fontsize=10)
    ax6.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    plt.setp(ax6.get_yticklabels(), fontsize=8)
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, numpoints=3)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    if cMF == None:
        # x axis
        ax6.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
        ax6.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
        ax6.xaxis.set_major_formatter(dateFmt)
        ax6.xaxis.set_minor_formatter(dateminorFmt)
        labels = ax6.get_xticklabels(which='both')
        plt.setp(labels, rotation=90, fontsize=8)
        del labels
        ax6.set_xlabel('Date', fontsize=10)
    else:
        plt.setp(ax6.get_xticklabels(which='both'), visible=False)
    # grd
    ax6.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax6.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    if MF is not None:
        # compute heads
        # plot GWT
        lines = itertools.cycle(
            ['-', '--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|',
             '_'])
        ax10 = fig.add_subplot(8, 1, 8, sharex=ax1)
        # ax10.set_autoscalex_on(False)
        for l in range(cMF.nlay):
            i = b'id_%d' % (l + 1)
            plt.plot_date(cMF.inputDate, flx[flxIndex_lst[i]], next(lines), color='b', markersize=2, markevery=7,
                          label=flxLbl[flxIndex_lst[i]])
        if obs_catch_list[0] == 1:
            obs_h = obs_catch.get('catch')['obs_h']
            hobs_m = np.ma.masked_values(obs_h[0], cMF.hnoflo, atol=0.09)
            dobs_m = hobs_m - TopSoilAverage
            ax10.plot_date(cMF.inputDate, dobs_m, 'o', markerfacecolor='None', markeredgecolor='LightBlue',
                           markersize=2, markevery=7, label=r'$d \ obs$')
        # y
        plt.ylabel('m', fontsize=10)
        plt.setp(ax10.get_yticklabels(), fontsize=8)
        # leg
        plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
                   columnspacing=colspc, numpoints=3)
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=8)
        ax10.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
        # x axis
        ax10.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
        ax10.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
        ax10.xaxis.set_major_formatter(dateFmt)
        ax10.xaxis.set_minor_formatter(dateminorFmt)
        labels = ax10.get_xticklabels(which='both')
        plt.setp(labels, rotation=90, fontsize=8)
        del labels
        ax10.set_xlabel('Date', fontsize=10)
        # grd
        ax10.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
        ax10.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax1.set_xlim(np.min(cMF.inputDate) - 15.0, np.max(cMF.inputDate) + 15)

    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig(plt_export_fn, dpi=150)
    plt.close()

    # --------------FIGURE CATCH PART2----------

    fig = plt.figure(num=None, figsize=(8.27, 11.7), dpi=150)  # (8.5,15), dpi=30)
    fig.suptitle(plt_title + ' - part 2')

    ax0 = fig.add_subplot(8, 1, 1)
    # ax0.set_autoscalex_on(False)
    # RF
    ax0.bar(cMF.inputDate, flx[flxIndex_lst[b'iRF']], color='darkblue', linewidth=0, align='center',
            label=flxLbl[flxIndex_lst[b'iRF']])
    # RFe
    ax0.bar(cMF.inputDate, flx[flxIndex_lst[b'iRFe']], color='deepskyblue', linewidth=0, align='center',
            label=flxLbl[flxIndex_lst[b'iRFe']])
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # y
    plt.ylabel('mm', fontsize=10)
    plt.setp(ax0.get_yticklabels(), fontsize=8)
    ax0.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    # x
    ax0.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax0.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    plt.setp(ax0.get_xticklabels(which='both'), visible=False)
    # grd
    ax0.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax0.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    # plot Ro
    ax2 = fig.add_subplot(8, 1, 2, sharex=ax0)
    # ax2.set_autoscalex_on(True)
    # Ro
    plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iRo']], '-', c='blue', linewidth=1.0, label=flxLbl[5])
    # Ro obs
    if obs_catch_list[2] == 1:
        obs_Ro = obs_catch.get('catch')['obs_Ro']
        Roobs_m = np.ma.masked_values(obs_Ro[0], cMF.hnoflo, atol=0.09)
        plt.plot_date(cMF.inputDate, Roobs_m, 'o', markerfacecolor='None', markeredgecolor='lightBlue',
                      markersize=2, label=r'$Ro \ obs$')  # markevery = 7,
    # y
    plt.ylabel('mm', fontsize=10)
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
    # leg
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
               columnspacing=colspc, numpoints=3)
    if obs_catch_list[2] == 1:
        plt.ylim(ymax=10.0*np.ma.max(Roobs_m))
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    # x axis
    ax2.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
    ax2.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
    ax2.xaxis.set_major_formatter(dateFmt)
    ax2.xaxis.set_minor_formatter(dateminorFmt)
    labels = ax2.get_xticklabels(which='both')
    plt.setp(labels, rotation=90, fontsize=8)
    del labels
    ax2.set_xlabel('Date', fontsize=10)
    # grd
    ax2.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
    ax2.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    rmseHEADS = None
    rsrHEADS = None
    nseHEADS = None
    rHEADS = None
    if MF is not None:
        # plot heads
        lines = itertools.cycle(
            ['-', '--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|',
             '_'])
        ax1 = fig.add_subplot(8, 1, 4)
        # ax1.set_autoscalex_on(True)
        # RMSE
        if obs_catch_list[0] == 1:
            i = b'ih_%d' % (l + 1)
            if sum(flx[flxIndex_lst[i]]) != 0.0:
                rmse, rsr, nse, r = cMF.cPROCESS.compCalibCrit(flx[flxIndex_lst[i]], obs_h[0], cMF.hnoflo)
                rmseHEADS = [rmse]
                rsrHEADS = [rsr]
                nseHEADS = [nse]
                rHEADS = [r]
                del rmse, rsr, nse, r
                if rmseHEADS[0] is not None:
                    print('h: %.2f m / %.2f / %.2f / %.2f\n-------' % (
                        rmseHEADS[0], rsrHEADS[0], nseHEADS[0], rHEADS[0]))
            else:
                print('Warning!\nERROR computing h calibration criteria')
                rmseHEADS = rsrHEADS = nseHEADS = rHEADS = None
        for l in range(cMF.nlay):
            i = b'ih_%d' % (l + 1)
            plt.plot_date(cMF.inputDate, flx[flxIndex_lst[i]], next(lines), color='b', markersize=2, markevery=7,
                          label=flxLbl[flxIndex_lst[i]])
        if obs_catch_list[0] == 1:
            ax1.plot_date(cMF.inputDate, hobs_m, 'o', markerfacecolor='None', markeredgecolor='LightBlue',
                          markersize=2, markevery=7, label=r'$h \ obs$')
        # y
        plt.ylim(hmin, hmax)
        plt.ylabel('m', fontsize=10)
        plt.setp(ax1.get_yticklabels(), fontsize=8)
        # leg
        plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=2,
                   columnspacing=colspc, numpoints=3)
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=8)
        # x
        ax1.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
        ax1.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
        plt.setp(ax1.get_xticklabels(which='both'), visible=False)
        # grd
        ax1.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
        ax1.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

        # plot GW fluxes
        ax8 = fig.add_subplot(2, 1, 2, sharex=ax1)
        # ax8.set_autoscalex_on(True)
        # uzf recharge
        plt.plot_date(cMF.inputDate, flx[flxIndex_lst[b'iRg']], '-', c='darkblue', linewidth=1.5,
                      label=flxLbl[flxIndex_lst[b'iRg']])
        lines = itertools.cycle(
            ['-', '--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|',
             '_'])
        for l, (e, lbl) in enumerate(zip(flx[flxIndex_lst[b'idSg_1']:], flxLbl[flxIndex_lst[b'idSg_1']:])):
            plt.plot_date(cMF.inputDate, e, next(lines), color=mpl.colors.rgb2hex(np.random.rand(1, 3)[0]),
                          markersize=2, markevery=7, label=lbl, markeredgecolor='None')
        # y
        plt.ylabel('mm', fontsize=10)
        plt.setp(ax8.get_yticklabels(), fontsize=8)
        ax8.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
        # leg
        plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale, borderpad=bdpd, handletextpad=hdltxtpd, ncol=4,
                   columnspacing=colspc, numpoints=3)
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=8)
        # x axis
        ax8.xaxis.set_major_locator(mpl.dates.YearLocator(1, month=iniMonthHydroYear, day=1))
        ax8.xaxis.set_minor_locator(mpl.dates.MonthLocator(bymonth=bymonth))
        ax8.xaxis.set_major_formatter(dateFmt)
        ax8.xaxis.set_minor_formatter(dateminorFmt)
        labels = ax8.get_xticklabels(which='both')
        plt.setp(labels, rotation=90, fontsize=8)
        del labels
        ax8.set_xlabel('Date', fontsize=10)
        # grd
        ax8.grid(b=True, which='major', axis='both', linestyle=':', color='darkgray')
        ax8.grid(b=True, which='minor', axis='x', linestyle=':', color='gainsboro')

    ax0.set_xlim(np.min(cMF.inputDate) - 15.0, np.max(cMF.inputDate) + 15)
    ax1.set_xlim(np.min(cMF.inputDate) - 15.0, np.max(cMF.inputDate) + 15)

    plt.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    txt = plt_export_fn.split('.')
    plt_export_fn = txt[0] + '_part2.' + txt[1]
    plt.savefig(plt_export_fn, dpi=150)
    #    plt.show()
    plt.clf()
    plt.close('all')
    return rmseHEADS, rmseSM, rsrHEADS, rsrSM, nseHEADS, nseSM, rHEADS, rSM


##################

def plotLAYER(days, str_per, Date, JD, ncol, nrow, nlay, nplot, V, cmap, CBlabel, msg, plt_title, MM_ws,
              interval_type=b'arange', interval_diff=1, interval_num=1, Vmax=0, Vmin=0, fmt=None, contours=False,
              ntick=1, facecolor='silver', points=None, ptslbl=0, mask=None, hnoflo=-999.9, animation=0,
              pref_plt_title='_sp_plt', cMF=None):
    # TODO put axes tick as row/col index from MODFLOW AND real coordinates

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
    x = np.arange(0.5, ncol + 1.5, 1)
    y = np.arange(0.5, nrow + 1.5, 1)
    xg, yg = np.meshgrid(x, y)

    x = np.arange(1, ncol + 1, 1)
    y = np.arange(1, nrow + 1, 1)
    xg1, yg1 = np.meshgrid(x, y)

    ims = []
    files_tmp = []
    # print(plt_title)
    # print("nplot: %s"% (nplot))
    for i, day in enumerate(days):
        # DEFINE HERE SUBPLOT FORMAT
        if nrow > ncol:
            figsize = (8.27, 11.7)
            NrowPage = 1
            if nlay > 1:
                NcolPage = 2
            else:
                NcolPage = 1
        else:
            figsize = (11.7, 8.27)
            NcolPage = 1
            if nlay > 1:
                NrowPage = 2
            else:
                NrowPage = 1
        NPage = int(np.ceil(nplot / 2.0))
        L = 0
        Vmin_tmp, Vmax_tmp, ctrs_tmp = MinMax(np.min(Vmin), np.max(Vmax), contours)
        if fmt == None:
            if Vmax_tmp > 0.0999 or abs(Vmin_tmp) > 0.0999:
                fmt = '%5.2f'
            else:
                fmt = '%5.e'
        if Vmax_tmp > 0 and Vmin_tmp < 0 and cMF is not None:
            cmap = plt.cm.coolwarm_r
            # shifted_cmap = cMF.cUTIL.remappedColorMap(cmap, midpoint=0.75, name='shifted')
            start = 0.0  # (Vmax_tmp-abs(Vmin_tmp))/(2*Vmax_tmp)
            midpoint = abs(Vmin_tmp) / (Vmax_tmp + abs(Vmin_tmp))
            stop = 1.0  # (abs(Vmin_tmp)-Vmax_tmp)/(2*abs(Vmin_tmp))
            shrunk_cmap = cMF.cUTIL.remappedColorMap(cmap, start=start, midpoint=midpoint, stop=stop, name='shrunk')
            cmap = shrunk_cmap
        norm = None
        if interval_type == b'arange':
            ticks = np.arange(Vmin_tmp, Vmax_tmp, interval_diff)
            levels = mpl.ticker.MaxNLocator(nbins=cmap.N).tick_values(Vmin_tmp, Vmax_tmp)
            norm = mpl.colors.BoundaryNorm(levels, cmap.N)
        elif interval_type == b'linspace':
            ticks = np.linspace(Vmin_tmp, Vmax_tmp, interval_num)
            levels = mpl.ticker.MaxNLocator(nbins=cmap.N).tick_values(Vmin_tmp, Vmax_tmp)
            norm = mpl.colors.BoundaryNorm(levels, cmap.N)  # , vmin=Vmin_tmp, vmax=Vmax_tmp)
        # print(Vmin_tmp, Vmax_tmp, interval_num, interval_diff, ticks, cmap.N)
        for F in range(NPage):
            fig = plt.figure(num=None, figsize=figsize, dpi=90)
            figtitle = fig.suptitle('')
            ims.append([])
            try:
                plt_title = plt_title.encode("utf-8")
            except (UnicodeDecodeError, AttributeError):
                pass
            if isinstance(Date[i], float):
                figtitle.set_text((plt_title + bytes('\nDate %s, DOY %s, stress period %s, day %d, Page %d/%d' % (
                mpl.dates.num2date(Date[i]).isoformat()[:10], JD[i], str_per[i], day + 1, F + 1, NPage),
                                                     "utf-8")).decode("utf-8"))
            else:
                figtitle.set_text((plt_title + bytes(', Page %d/%d' % (F + 1, NPage), "utf-8")).decode("utf-8"))
            # plt.draw()  # TODO confirmar impacto desta linha
            ax = []
            for l in range(2):
                if L < nplot:
                    # print("layer %d, Layer %d, Page %d" % (l, L, F))
                    if mask.all() == None:
                        Vtmp = V[i, L, :, :]
                    else:
                        Vtmp = np.ma.masked_array(V[i, L, :, :], mask[L, :, :])
                    Vtmp = np.ma.masked_values(Vtmp, hnoflo, atol=0.09)
                    if interval_type == b'percentile':
                        ticks = np.percentile(Vtmp.compressed().flatten(),
                                              [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0])
                        levels = mpl.ticker.MaxNLocator(nbins=cmap.N).tick_values(0.0, 100.0)
                        norm = mpl.colors.BoundaryNorm(ticks, cmap.N)  # , vmin=Vmin_tmp, vmax=Vmax_tmp)
                    ax.append(fig.add_subplot(NrowPage, NcolPage, l + 1, facecolor=facecolor))
                    ax[l].xaxis.set_ticks(np.arange(0, ncol + 1, ntick))
                    plt.setp(ax[l].get_xticklabels(), fontsize=8)
                    if l < 1:
                        ax[l].yaxis.set_ticks(np.arange(0, nrow + 1, ntick))
                        plt.setp(ax[l].get_yticklabels(), fontsize=8)
                        plt.ylabel('row i', fontsize=10)
                        ax[l].yaxis.set_label_position("right")
                    else:
                        ax[l].set_yticklabels([])
                    ax[l].yaxis.tick_right()
                    ax[l].yaxis.set_ticks_position('both')
                    plt.xlabel('col j', fontsize=10)
                    ax[l].xaxis.set_label_position("top")
                    ax[l].xaxis.tick_top()
                    ax[l].xaxis.set_ticks_position('both')
                    if points is not None:
                        for k, (xj, yi, lay, label) in enumerate(zip(points[2], points[1], points[3], points[0])):
                            if lay == L:
                                color = 'dimgrey'
                            else:
                                color = 'lightgrey'
                            ax[l].plot(xj, yi, 'o', linewidth=1, markersize=4, color=color)
                            if ptslbl > 0:
                                ax[l].annotate(label, xy=(xj, yi - 0.15), fontsize=8, ha='center', va='bottom')
                    ims[F].append(ax[l].pcolormesh(xg, yg, Vtmp, cmap=cmap, norm=norm))
                    if ctrs_tmp == True:
                        try:
                            CS = ax[l].contour(xg1, yg1[::-1], Vtmp[::-1], ticks, colors='gray')
                            plt.draw()
                            ax[l].clabel(CS, inline=1, fontsize=6, fmt=fmt, colors='gray')
                        except:
                            print('Error in drawing contours for map %s' % plt_title)
                    if np.ma.max(Vtmp) > np.ma.min(Vtmp):
                        ax[l].set_title('layer %d' % (L + 1), fontsize=10, y=-0.1, fontweight='bold')
                    else:
                        ax[l].set_title('layer %d %s' % (L + 1, msg), fontsize=10, y=-0.1, fontweight='bold')
                    ax[l].set_ylim(bottom=np.max(yg1), top=np.min(yg1))
                    ax[l].axis('scaled')
                    axl, axb, axw, axh = ax[l].get_position().bounds
                    if interval_type == b'percentile':
                        if max(x) > max(y):
                            if l == 0:
                                # rect : sequence of float
                                # The dimensions [left, bottom, width, height] of the new axes. All quantities are in fractions of figure width and height.
                                # cax = fig.add_axes([0.125, 0.035, 0.75, 0.025])
                                cax = fig.add_axes([axl, 0.035, axw, 0.025])
                            else:
                                # cax = fig.add_axes([0.625, 0.035, 0.75, 0.025])
                                cax = fig.add_axes([axl, 0.035, axw, 0.025])
                            CBorient = 'horizontal'
                            # cax.xaxis.set_label_position('top')
                        else:
                            if l == 0:
                                # cax = fig.add_axes([0.005, 0.125, 0.025, 0.75])
                                cax = fig.add_axes([0.035, axb, 0.025, axh])
                                # cax.yaxis.set_label_position('left')
                            else:
                                # cax = fig.add_axes([0.925, 0.125, 0.025, 0.75])
                                cax = fig.add_axes([axl + axw + 0.015, axb, 0.025, axh])
                                # cax.yaxis.set_label_position('right')
                            CBorient = 'vertical'
                        CB = fig.colorbar(ims[F][0 + l], ticks=ticks, extend='both', format=fmt, cax=cax,
                                          orientation=CBorient)  # , shrink=0.8)
                        if l == 0:
                            CB.set_label(CBlabel, fontsize=10)  # , loc='center')
                        plt.setp(CB.ax.get_xticklabels(), fontsize=7)
                        plt.setp(CB.ax.get_yticklabels(), fontsize=7)
                        del cax
                    L += 1

            if interval_type != b'percentile':
                if max(x) > max(y):
                    # cax = fig.add_axes([0.125, 0.035, 0.75, 0.025])
                    cax = fig.add_axes([axl, 0.035, axw, 0.025])
                    CBorient = 'horizontal'
                else:
                    # cax = fig.add_axes([0.035, 0.125, 0.025, 0.75])
                    cax = fig.add_axes([0.035, axb, 0.025, axh])
                    CBorient = 'vertical'
                CB = fig.colorbar(ims[F][0], ticks=ticks, extend='both', format=fmt, cax=cax,
                                  orientation=CBorient)  # , shrink=0.8)
                # print(ticks)
                CB.set_label(CBlabel, fontsize=10)
                if max(x) > max(y):
                    cax.xaxis.set_label_position('top')
                    plt.setp(CB.ax.get_xticklabels(), fontsize=7)
                else:
                    cax.yaxis.set_label_position('left')
                    plt.setp(CB.ax.get_yticklabels(), fontsize=7)
                del cax

            if isinstance(Date[i], float):
                plt_export_fn = os.path.join(MM_ws, '%s_%s_day%05d_%s_%s.png' % (
                    pref_plt_title, plt_title, day + 1, F + 1, NPage))
            else:
                plt_export_fn = os.path.join(MM_ws, '%s_%s_%s_%s.png' % (pref_plt_title, plt_title, F + 1, NPage))
            plt.savefig(plt_export_fn)
            # print("Printed %s" % plt_export_fn)
            if len(days) > 1 and animation == 1:
                try:
                    if i == 0:
                        files_tmp.append(os.path.join(MM_ws, '%05d.png' % (0)))
                        print(plt_export_fn, files_tmp[0])
                        shutil.copyfile(plt_export_fn, files_tmp[0])
                    else:
                        files_tmp.append(os.path.join(MM_ws, '%05d.png' % (i)))
                        shutil.copyfile(plt_export_fn, files_tmp[i])
                except:
                    pass
            for l in range(len(ax)):
                ax[l].cla()
    # TODO correct to produce movies for each pages
    if len(days) > 1 and animation == 1:
        batch_fn = os.path.join(MM_ws, 'run.bat')
        f = open(batch_fn, 'w')
        f.write(
            'ffmpeg -r 1 -i %s -s:v 1280x720 -c:v libx264 -profile:v high -crf 23 -pix_fmt yuv420p -r 30 -y %s_mov.mp4' % (
                '%s\%%%%05d.png' % (MM_ws), '%s\%s_%s' % (MM_ws, pref_plt_title, plt_title)))
        f.close()
        run_report_fn = os.path.join(MM_ws, '__FFmpegRunReport.txt')
        run_report = open(run_report_fn, 'w')
        sp.Popen(batch_fn, shell=False, stdout=run_report, stderr=run_report).wait()
        run_report.close()
        os.remove(batch_fn)
        try:
            for f in files_tmp:
                os.remove(f)
        except:
            pass

    plt.close('all')


##################

def plotWBsankey(path, DATE, flx, flxIndex, fn, indexTime, year_lst, cMF, ncell_MM, obspt, fntitle, ibound4Sankey,
                 stdout=None, report=None):
    """ Computes the water balance for a certain time span
    input: ASCII file with water fluxes wrtitten by MM
    """

    # compute fluxes for whole modelled period and hydrological years
    RF = []
    I = []
    RFe = []
    dSsurf = []
    Ro = []
    Esurf = []
    dSsoil = []
    EXFtotMM = []
    Exf_l0 = []
    Esoil = []
    Tsoil = []
    ETsoil = []
    Inf = []
    Eg = []
    Tg = []
    Egtot = []
    Tgtot = []
    ETg = []
    Ssurf = []
    Rp = []
    dSu = []
    Rg = []
    dSg = []
    FRF = []
    FFF = []
    FLF = []
    EXF = []
    #EXFtotMF = []
    WEL = []
    DRN = []
    GHB = []
    CH = []
    for k, i in enumerate(indexTime[:-2]):
        if k == 0:
            i = indexTime[1]
            indexend = indexTime[-2]
            mult = 365.0 / (indexend - i + 1)
        else:
            indexend = indexTime[k + 1] - 1
            mult = 1.0
        RF.append(mult * np.sum(np.float16(flx[flxIndex[b'iRF']][i:indexend])))
        I.append(mult * np.sum(np.float16(flx[flxIndex[b'iI']][i:indexend])))
        RFe.append(mult * np.sum(np.float16(flx[flxIndex[b'iRFe']][i:indexend])))
        dSsurf.append(mult * np.sum(np.float16(flx[flxIndex[b'idSsurf']][i:indexend])))
        Ro.append(mult * np.sum(np.float16(flx[flxIndex[b'iRo']][i:indexend])))
        Esurf.append(mult * np.sum(np.float16(flx[flxIndex[b'iEsurf']][i:indexend])))
        dSsoil.append(mult * np.sum(np.float16(flx[flxIndex[b'idSsoil']][i:indexend])))
        EXFtotMM.append(mult * np.sum(np.float16(flx[flxIndex[b'iEXFg']][i:indexend])))
        Exf_l0.append(mult * np.sum(np.float16(flx[flxIndex[b'iExf_1']][i:indexend])))
        Esoil.append(mult * np.sum(np.float16(flx[flxIndex[b'iEsoil']][i:indexend])))
        Tsoil.append(mult * np.sum(np.float16(flx[flxIndex[b'iTsoil']][i:indexend])))
        ETsoil.append(mult * np.sum(np.float16(flx[flxIndex[b'iETsoil']][i:indexend])))
        Inf.append(mult * np.sum(np.float16(flx[flxIndex[b'iInf']][i:indexend])))
        if cMF.wel_yn == 1:
            Eg.append(np.zeros(cMF.Mnlay))
            Tg.append(np.zeros(cMF.Mnlay))
            for ii, (ML, L) in enumerate(zip(cMF.Mlay, range(cMF.nlay))):
                Eg[k][ML-1] += mult * np.sum(np.float16(flx[flxIndex[b'iEg_%d' % (L + 1)]][i:indexend]))
                Tg[k][ML-1] += mult * np.sum(np.float16(flx[flxIndex[b'iTg_%d' % (L + 1)]][i:indexend]))
            Egtot.append(mult * np.sum(np.float16(flx[flxIndex[b'iEg']][i:indexend])))
            Tgtot.append(mult * np.sum(np.float16(flx[flxIndex[b'iTg']][i:indexend])))
            ETg.append(mult * np.sum(np.float16(flx[flxIndex[b'iETg']][i:indexend])))
        Ssurf.append(mult * np.sum(np.float16(flx[flxIndex[b'iSsurf']][i:indexend])))
        Rp.append(mult * np.sum(np.float16(flx[flxIndex[b'iperc']][i:indexend])))
        dSu.append(mult * np.sum(np.float16(flx[flxIndex[b'idSu']][i:indexend])))
        Rg.append(np.zeros(cMF.Mnlay))
        dSg.append(np.zeros(cMF.Mnlay))
        FRF.append(np.zeros(cMF.Mnlay))
        FFF.append(np.zeros(cMF.Mnlay))
        FLF.append(np.zeros(cMF.Mnlay))
        EXF.append(np.zeros(cMF.Mnlay))
        WEL.append(np.zeros(cMF.Mnlay))
        DRN.append(np.zeros(cMF.Mnlay))
        GHB.append(np.zeros(cMF.Mnlay))
        CH.append(np.zeros(cMF.Mnlay))
        for ii, (ML, L) in enumerate(zip(cMF.Mlay, range(cMF.nlay))):
            Rg[k][ML - 1] += mult * np.sum(np.float16(flx[flxIndex[b'iRg_%d' % (L + 1)]][i:indexend]))
            dSg[k][ML - 1] += mult * np.sum(np.float16(flx[flxIndex[b'idSg_%d' % (L + 1)]][i:indexend]))
            FRF[k][ML - 1] += mult * np.sum(np.float16(flx[flxIndex[b'iFRF_%d' % (L + 1)]][i:indexend]))
            FFF[k][ML - 1] += mult * np.sum(np.float16(flx[flxIndex[b'iFFF_%d' % (L + 1)]][i:indexend]))
            if cMF.nlay > 1:
                FLF[k][ML - 1] += mult * np.sum(np.float16(flx[flxIndex[b'iFLF_%d' % (L + 1)]][i:indexend]))
            EXF[k][ML - 1] += mult * np.sum(np.float16(flx[flxIndex[b'iEXFg_%d' % (L + 1)]][i:indexend]))
            if cMF.wel_yn == 1:
                if ncell_MM[L] > 0:
                    WEL[k][ML - 1] += mult * np.sum(np.float16(flx[flxIndex[b'iWEL_%d' % (L + 1)]][i:indexend]))
            if cMF.drn_yn == 1:
                if cMF.drncells[L] > 0:
                    DRN[k][ML - 1] += mult * np.sum(np.float16(flx[flxIndex[b'iDRN_%d' % (L + 1)]][i:indexend]))
            if cMF.ghb_yn == 1:
                if cMF.ghbcells[L] > 0:
                    GHB[k][ML - 1] += mult * np.sum(np.float16(flx[flxIndex[b'iGHB_%d' % (L + 1)]][i:indexend]))
            if len(cMF.ibound[cMF.ibound < 0]) > 0:
                CH[k][ML - 1] += mult * np.sum(np.float16(flx[flxIndex[b'iCH_%d' % (L + 1)]][i:indexend]))
        #EXFtotMF[k] += sum(EXF[k])
    #    print "\nWater fluxes imported from file:\n%s" % inputFile_fn

    # Sankey plots
    #    prt_test = 0
    for f in [0, 1]:
        for k in range(len(RF)):
            # print '-------'
            if k == 0:
                title = "Average of the %d hydrological year(s)" % len(indexTime[1:-2])
            else:
                title = "Hydrological year %d/%d" % (year_lst[k - 1], year_lst[k - 1] + 1)
            if f == 0:
                ff = 1.0
                fff = (1.5 * RF[k])
            #                if prt_test == 0:
            #                    print '\nWater balance (mm.y-1)'
            #                    prt_test += 1
            else:
                ff = RF[k] / 100.0
                fff = (1.5 * RF[k]) / ff
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
            ##print(title)
            treshold = 5E-2  # 5E-2
            if k == 0 or k % 2 != 0:
                fig = plt.figure(figsize=(8.27, 11.7), dpi=72)
                if f == 0:
                    fig.suptitle('Water balance ($mm.y^{-1}$) - %s\n' % obspt, fontsize=10, y=0.99)
                else:
                    fig.suptitle('Water balance ($\%%$ of yearly rainfall) - %s\n' % obspt, fontsize=10, y=0.99)
                ax = []
                ax.append(fig.add_subplot(2, 1, 1, xticks=[], yticks=[]))
                y = 0.530
                p = 0
            else:
                ax.append(fig.add_subplot(2, 1, 2, xticks=[], yticks=[]))
                y = 0.085
                p = 1
            figtitle = fig.text(x=0.5, y=y, s=title, horizontalalignment='center', verticalalignment='bottom',
                                fontsize=8)
            fmt = '%.1f'
            pltsankey = Sankey(ax=ax[p], format=fmt, scale=1.0 / fff, offset=0.25, gap=0.5, shoulder=0.0, margin=0.5)
            pl = 0.5
            tl = 2.0
            # MMveg
            pltsankey.add(patchlabel='veg', facecolor='lightgreen', trunklength=tl / 1.5,
                          flows=[RF[k] / ff, -I[k] / ff, -RFe[k] / ff],
                          labels=['$RF$', '$I$', '$RFe$'],
                          orientations=[1, 1, -1],
                          pathlengths=[pl, pl, pl])
            # MMsurf
            flows = [RFe[k] / ff, -Inf[k] / ff, Exf_l0[k] / ff, -Esurf[k] / ff, -Ro[k] / ff]
            if np.abs(Exf_l0[k]) > treshold:
                #flows = [RFe[k] / ff, -Inf[k] / ff, Exf_l0[k] / ff, -Esurf[k] / ff, -Ro[k] / ff]
                labels = [None, '$Inf$', '$Exf_1$', '$E_{surf}$', '$Ro$']
                orientations = [1, -1, -1, 1, 0]
                pathlengths = [pl, pl, pl, pl, pl]
            else:
                #flows = [RFe[k] / ff, -Inf[k] / ff, 0.0, -Esurf[k] / ff, -Ro[k] / ff]
                labels = [None, '$Inf$', '', '$E_{surf}$', '$Ro$']
                orientations = [1, -1, -1, 1, 0]
                pathlengths = [pl, pl, 0, pl, pl]
            pltsankey.add(patchlabel='$\Delta S_{surf}$\n%.1f' % (-dSsurf[k] / ff), label='MMsurf',
                          facecolor='lightblue', trunklength=tl,
                          flows=flows,
                          labels=labels,
                          orientations=orientations,
                          pathlengths=pathlengths,
                          prior=0, connect=(2, 0))
            In = RFe[k] + Exf_l0[k]
            Out = Inf[k] + Esurf[k] + Ro[k]
            if dSsurf[k] > 0.0:
                Out += dSsurf[k]
            else:
                In += -dSsurf[k]
            MB_MMsurf = 100 * (In - Out) / ((In + Out) / 2)
            # MMsoil
            if EXFtotMM[k] > treshold:
                flows = [Inf[k] / ff, -Rp[k] / ff, -Exf_l0[k] / ff, -Esoil[k] / ff, -Tsoil[k] / ff,
                         EXFtotMM[k] / ff]
                if np.abs(-Exf_l0[k]) > treshold:
                    #flows = [Inf[k] / ff, -Rp[k] / ff, -Exf_l0[k] / ff, -Esoil[k] / ff, -Tsoil[k] / ff,
                             #EXFtotMM[k] / ff]
                    labels = [None, '$R_p$', '$Exf_1$', '$E_{soil}$', '$T_{soil}$', '$Exf_g$']
                    orientations = [1, -1, 1, 1, 1, -1]
                    pathlengths = [pl, pl, pl, pl, pl, pl]
                else:
                    #flows = [Inf[k] / ff, -Rp[k] / ff, 0.0, -Esoil[k] / ff, -Tsoil[k] / ff, EXFtotMM[k] / ff]
                    labels = [None, '$R_p$', '', '$E_{soil}$', '$T_{soil}$', '$Exf_g$']
                    orientations = [1, -1, 1, 1, 1, -1]
                    pathlengths = [pl, pl, 0, pl, pl, pl]
            else:
                flows = [Inf[k] / ff, -Rp[k] / ff, -Exf_l0[k] / ff, -Esoil[k] / ff, -Tsoil[k] / ff]
                if np.abs(EXF[k][0]) > treshold:
                    #flows = [Inf[k] / ff, -Rp[k] / ff, -Exf_l0[k] / ff, -Esoil[k] / ff, -Tsoil[k] / ff]
                    labels = [None, '$R_p$', '$Exf_1$', '$E_{soil}$', '$T_{soil}$']
                    orientations = [1, -1, 1, 1, 1]
                    pathlengths = [pl, pl, pl, pl, pl]
                else:
                    #flows = [Inf[k] / ff, -Rp[k] / ff, 0.0, -Esoil[k] / ff, -Tsoil[k] / ff]
                    labels = [None, '$R_p$', '', '$E_{soil}$', '$T_{soil}$']
                    orientations = [1, -1, 1, 1, 1]
                    pathlengths = [pl, pl, 0, pl, pl]
            pltsankey.add(patchlabel='$\Delta S_{soil}$\n%.1f' % (dSsoil[k] / ff),
                          label='MMsoil', facecolor='khaki', trunklength=tl,
                          flows=flows,
                          labels=labels,
                          orientations=orientations,
                          pathlengths=pathlengths,
                          prior=1, connect=(1, 0))
            In = Inf[k] + EXFtotMM[k]
            Out = Rp[k] + Esoil[k] + Tsoil[k] + Exf_l0[k]
            if dSsoil[k] > 0.0:
                Out += dSsoil[k]
            else:
                In += -dSsoil[k]
            MB_MMsoil = 100 * (In - Out) / ((In + Out) / 2)
            # MFuzf
            flows = [Rp[k] / ff]
            labels = [None]
            orientations = [1]
            pathlengths = [2 * pl]
            for L in range(cMF.Mnlay):
                if ibound4Sankey[L] > 0:
                    flows.append(-Rg[k][L] / ff)
                    if Rg[k][L] / ff > treshold:
                        #flows.append(-Rg[k][L] / ff)
                        labels.append('$Rg_%d$' % (L + 1))
                    else:
                        #flows.append(0.0)
                        labels.append('')
                    orientations.append(-1)
                    pathlengths.append(pl)
            pltsankey.add(patchlabel='$\Delta S_p$\n%.1f' % (dSu[k] / ff), label='MF_UZF', facecolor='lavender',
                          trunklength=tl,
                          flows=flows,
                          labels=labels,
                          orientations=orientations,
                          pathlengths=pathlengths,
                          prior=2, connect=(1, 0))
            In = Rp[k]
            Out = 0
            for L in range(cMF.Mnlay):
                Out += Rg[k][L]
            if dSu[k] > 0.0:
                Out += dSu[k]
            else:
                In += -dSu[k]
            MB_MFuzf = 100 * (In - Out) / ((In + Out) / 2)
            # MF
            # about signs, read MF-2005 manual pag 3-10
            MB_MF = []
            ##print('\nSANKEY TIME!--------\n%s'%obspt)
            ##print("ibound4Sankey: %s" % ibound4Sankey)
            L_act = 0
            lst_colors = []
            for L in range(cMF.Mnlay):
                lst_colors.append(getattr(mpl.cm, 'Blues')(L / cMF.Mnlay))
            colors = itertools.cycle(lst_colors)
            for L in range(cMF.Mnlay):
                if ibound4Sankey[L] > 0:
                    ##print("----")
                    ##print("L: %s, L_act: %s" % (L, L_act))
                    In = []
                    Out = []
                    # first layer
                    if L_act == 0:
                        ##print("First layer")
                        if L == (cMF.Mnlay - 1):
                            if cMF.Mnlay > 1:
                                flows = [Rg[k][L] / ff, -FLF[k][L] / ff, -FRF[k][L] / ff, -FFF[k][L] / ff]
                                labels = [None, '$FLF$', '$FRF$', '$FFF$']
                                orientations = [1, -1, 0, 0]
                                pathlengths = [pl * 4, pl, pl, pl * 4]
                            else:
                                flows = [Rg[k][L] / ff, -FRF[k][L] / ff, -FFF[k][L] / ff]
                                labels = [None, '$FRF$', '$FFF$']
                                orientations = [1, 0, 0]
                                pathlengths = [pl * 4, pl, pl * 4]
                        else:
                            flows = [Rg[k][L] / ff, -FLF[k][L] / ff, -FRF[k][L] / ff, -FFF[k][L] / ff]
                            labels = [None, '$FLF$', '$FRF$', '$FFF$']
                            orientations = [1, -1, 0, 0]
                            pathlengths = [pl * 4, pl, pl, pl * 4]
                        connect = (1, 0)
                        if cMF.Mnlay > 1:
                            if FLF[k][L] > 0.0:
                                Out.append(FLF[k][L])
                            else:
                                In.append(-FLF[k][L])
                        facecolor = next(colors)
                        ##print("labels: %s" % labels)
                        ##print("flows: %s" % flows)
                    # last layer
                    elif L_act == (cMF.Mnlay - 1):
                        ##print("Last layer")
                        flows = [FLF[k][L - 1] / ff, Rg[k][L] / ff, -FRF[k][L] / ff, -FFF[k][L] / ff]
                        if np.abs(Rg[k][L]) > treshold:
                            #flows = [FLF[k][L - 1] / ff, Rg[k][L] / ff, -FRF[k][L] / ff, -FFF[k][L] / ff]
                            labels = [None, '$Rg_%d$' % (L + 1), '$FRF$', '$FFF$']
                            orientations = [1, 1, 0, 0]
                            pathlengths = [pl * 4, pl, pl, pl * 4]
                            #connect = (2, 0)
                            connect = (1, 0)
                        else:
                            #flows = [FLF[k][L - 1] / ff, 0.0, -FRF[k][L] / ff, -FFF[k][L] / ff]
                            labels = [None, '', '$FRF$', '$FFF$']
                            orientations = [1, 1, 0, 0]
                            pathlengths = [pl * 4, pl, pl, pl * 4]
                            #if np.abs(Rg[k][L - 1]) > treshold:
                            #    connect = (2, 0)
                            #else:
                            #    connect = (1, 0)
                        if FLF[k][L - 1] > 0.0:
                            In.append(FLF[k][L - 1])
                        else:
                            Out.append(-FLF[k][L - 1])
                        ##print("labels: %s" % labels)
                        ##print("flows: %s" % flows)
                        facecolor = next(colors)
                    else:
                        # intermediary layers
                        ##print("intermediary layer")
                        flows = [FLF[k][L - 1] / ff, Rg[k][L] / ff, -FLF[k][L] / ff, -FRF[k][L] / ff,
                                 -FFF[k][L] / ff]
                        if np.abs(Rg[k][L]) > treshold:
                            #flows = [FLF[k][L - 1] / ff, Rg[k][L] / ff, -FLF[k][L] / ff, -FRF[k][L] / ff,
                                     #-FFF[k][L] / ff]
                            labels = [None, None, '$FLF$', '$FRF$',
                                      '$FFF$']  # labels=[None, None, '$\Delta S_g$', '$FLF$', '$FRF$', '$FFF$']
                            orientations = [1, 1, -1, 0, 0]
                            pathlengths = [pl, pl, pl, pl, pl * 4]
                            #rch = True
                        else:
                            #flows = [FLF[k][L - 1] / ff, 0.0, -FLF[k][L] / ff, -FRF[k][L] / ff, -FFF[k][L] / ff]
                            labels = [None, '$FLF$', '$FRF$',
                                      '$FFF$']  # labels=[None, None, '$\Delta S_g$', '$FLF$', '$FRF$', '$FFF$']
                            orientations = [1, 1, -1, 0, 0]
                            pathlengths = [pl, pl, pl, pl, pl * 4]
                            #rch = False
                        if L_act == 1:
                            connect = (1, 0)  # (prior, this)
                        else:
                            if np.abs(Rg[k][L - 1]) > treshold:
                                connect = (2, 0)
                            else:
                                connect = (1, 0)
                        if FLF[k][L] > 0.0:
                            Out.append(FLF[k][L])
                        else:
                            In.append(-FLF[k][L])
                        if FLF[k][L - 1] > 0.0:
                            In.append(FLF[k][L - 1])
                        else:
                            Out.append(-FLF[k][L - 1])
                        facecolor = next(colors)
                        ##print("labels: %s" % labels)
                        ##print("flows: %s" % flows)
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
                    if dSg[k][L] > 0.0:
                        In.append(dSg[k][L])
                    else:
                        Out.append(-dSg[k][L])
                    # EXF
                    flows.append(EXF[k][L] / ff)
                    orientations.append(1)
                    if np.abs(EXF[k][L]) > treshold:
                        #flows.append(EXF[k][L] / ff)
                        labels.append('$Exf_g$')
                        #orientations.append(1)
                        pathlengths.append(pl)
                    else:
                        labels.append('')
                        pathlengths.append(pl)
                    Out.append(-EXF[k][L])
                    # Eg
                    if cMF.wel_yn == 1:
                        flows.append(-Eg[k][L] / ff)
                        orientations.append(1)
                        if np.abs(Eg[k][L]) > treshold:
                            #flows.append(-Eg[k][L] / ff)
                            labels.append('$E_g$')
                            #orientations.append(1)
                            pathlengths.append(pl)
                        else:
                            labels.append('')
                            pathlengths.append(pl)
                        Out.append(Eg[k][L])
                    # Tg   
                    if cMF.wel_yn == 1:
                        flows.append(-Tg[k][L] / ff)
                        orientations.append(1)
                        if np.abs(Tg[k][L]) > treshold:
                            #flows.append(-Tg[k][L] / ff)
                            labels.append('$T_g$')
                            #orientations.append(1)
                            pathlengths.append(pl)
                        else:
                            labels.append('')
                            pathlengths.append(pl)
                        Out.append(Tg[k][L])
                    # DRN
                    if cMF.drn_yn == 1:
                        flows.append(DRN[k][L] / ff)
                        orientations.append(1)
                        if np.abs(DRN[k][L]) > treshold:
                            #flows.append(DRN[k][L] / ff)
                            labels.append('$DRN$')
                            #orientations.append(1)
                            pathlengths.append(pl)
                        else:
                            labels.append('')
                            pathlengths.append(pl)
                        Out.append(-DRN[k][L])
                    # GHB
                    if cMF.ghb_yn == 1:
                        flows.append(GHB[k][L] / ff)
                        orientations.append(1)
                        if np.abs(GHB[k][L]) > treshold:
                            #flows.append(GHB[k][L] / ff)
                            labels.append('$GHB$')
                            #orientations.append(1)
                            pathlengths.append(pl)
                        else:
                            labels.append('')
                            pathlengths.append(pl)
                        if GHB[k][L] < 0.0:
                            Out.append(-GHB[k][L])
                        else:
                            In.append(GHB[k][L])
                    # CH
                    if len(cMF.ibound[cMF.ibound < 0]) > 0:
                        flows.append(CH[k][L] / ff)
                        orientations.append(1)
                        if np.abs(CH[k][L]) > treshold:
                            #flows.append(CH[k][L] / ff)
                            labels.append('$CH$')
                            #orientations.append(1)
                            pathlengths.append(pl)
                        else:
                            labels.append('')
                            pathlengths.append(pl)
                        if CH[k][L] < 0.0:
                            Out.append(-CH[k][L])
                        else:
                            In.append(CH[k][L])
                    ##print("labels: %s" % labels)
                    ##print("flows: %s" % flows)
                    ##print("connect: %s" % (connect,))
                    ##print("facecolor: %s" % (facecolor,))
                    ##print("dSg %s" % (dSg[k][L]))
                    pltsankey.add(patchlabel='$\Delta S_g$\n%.1f' % (-dSg[k][L] / ff), label='MF_ML%d' % (L + 1),
                                  facecolor=facecolor, trunklength=tl, flows=flows, labels=labels,
                                  orientations=orientations, pathlengths=pathlengths, prior=3 + L_act, connect=connect)
                    MB_MF.append(100 * (sum(In) - sum(Out)) / ((sum(In) + sum(Out)) / 2))
                    L_act += 1
                else:
                    MB_MF.append(np.nan)
                ##print("L_act: %s" % L_act)
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
            if p == 0:
                lblspc = 0.05
                mkscl = 0.1
                handletextpad = 0.25
                borderaxespad = 0.5
                colspc = 0.05
                handlelength = 1.0
                plt.legend(loc='best', labelspacing=lblspc, markerscale=mkscl, handlelength=handlelength,
                           handletextpad=handletextpad, borderaxespad=borderaxespad, ncol=2,
                           columnspacing=colspc)  # loc='upper right'
                leg = plt.gca().get_legend()
                ltext = leg.get_texts()
                plt.setp(ltext, fontsize=6)
            # water balance box
            # msg = 'Water balance\nclosure (%%)\n-----------------\nMMsurf =%3.1f\nMMsoil =%3.1f\nMFUZF =%3.1f' % (MB_MMsurf, MB_MMsoil, MB_MFuzf) #
            msg = 'Water balance closure (%%)\nMMsurf=%3.1f, MMsoil=%3.1f, MFUZF=%3.1f' % (
                MB_MMsurf, MB_MMsoil, MB_MFuzf)  #
            for L in range(cMF.Mnlay):
                msg += ' // MFL%d=%3.1f' % (L + 1, MB_MF[L])  # %2.1f
            xy = (figtitle.get_position()[0], figtitle.get_position()[1] - 0.005)
            plt.annotate(msg, xy, horizontalalignment='center', verticalalignment='top', fontsize=6,
                         xycoords='figure fraction')
            if k == 0:
                # print('Plot of the whole modelled period done!')
                msg = '_0whole'
            else:
                # print('Plot of hydrological year %d/%d done!' % (year_lst[k-1], year_lst[k-1] + 1))
                msg = '_%d_%d' % (year_lst[k - 1], year_lst[k - 1] + 1)
            # export png
            if k % 2 == 0 or k == len(RF) - 1:
                plt.subplots_adjust(left=0.05, bottom=0.10, right=0.95, top=0.95, wspace=0.01, hspace=0.1)
                if f == 0:
                    plt_export_fn = os.path.join(path, '_%s_WBsanley%s.png' % (fntitle, msg))
                else:
                    plt_export_fn = os.path.join(path, '_%s_WBsanley%s_pc.png' % (fntitle, msg))
                plt.savefig(plt_export_fn, dpi=150)
        plt.close('all')
    # print '-------'


##################

def plotCALIBCRIT(calibcritSM, calibcritSMobslst, calibcritHEADS, calibcritHEADSobslst, calibcritHEADSc,
                  calibcritHEADScobslst, plt_export_fn, plt_title, calibcrit, calibcritSMmax=None,
                  calibcritHEADSmax=None, ymin=None, units='', hnoflo=-9999.9):
    # plot RMSE

    num_plt_max = 6
    if len(calibcritSM) > num_plt_max or len(calibcritHEADS) > num_plt_max:
        num_plt_tmp = []
        num_plt_tmp.append(int(np.ceil(len(calibcritSM) / float(num_plt_max))))
        num_plt_tmp.append(int(np.ceil(len(calibcritHEADS) / float(num_plt_max))))
        num_plt = np.max(num_plt_tmp)
        del num_plt_tmp
    else:
        num_plt = 1

    p_mult = 0
    for p in range(num_plt):
        fig = plt.figure()
        fig.suptitle(plt_title, fontsize=10)
        test = 0
        calibcritSM_tmp = calibcritSM[p_mult:(p_mult + 1) * num_plt_max]
        calibcritSMobslst_tmp = calibcritSMobslst[p_mult:(p_mult + 1) * num_plt_max]
        if len(calibcritSM_tmp) > 0:
            test += 1
            ax1 = fig.add_subplot(2, 1, 1)
            ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            plt.ylabel('%s soil moisture %s' % (calibcrit, units[1]), fontsize=10, horizontalalignment='center')
            plt.grid(True)
            xserie = []
            yserie = []
            labels = []
            if calibcritSMmax == None:
                tmp = np.max(list(itertools.chain.from_iterable(calibcritSM_tmp)))
                if tmp > 0:
                    max_tmp = 1.2 * tmp
                else:
                    max_tmp = 1.0
                del tmp
            else:
                max_tmp = calibcritSMmax
            n = 0
            for e in calibcritSM_tmp:
                if len(e) > 1:
                    numtick = len(e)
                    if n == 0:
                        newx = 1.0
                    else:
                        newx = 1.0 + xserie[-1]
                    ee = 0
                    lst = list(range(1, int(numtick / 2 + 1)))
                    lst.reverse()
                    for i in lst:
                        xserie.append(newx - i / 10.0)
                        yserie.append(e[ee])
                        ee += 1
                        labels.append('')
                    xserie.append(newx)
                    labels.append('%s' % calibcritSMobslst_tmp[n])
                    if numtick % 2 != 0:
                        yserie.append(e[ee])
                        ee += 1
                    else:
                        yserie.append(max_tmp * 2.0)
                    lst = list(range(1, int(numtick / 2 + 1)))
                    for i in lst:
                        xserie.append(newx + i / 10.0)
                        yserie.append(e[ee])
                        ee += 1
                        labels.append('')
                else:
                    if n == 0:
                        xserie.append(1.0)
                        yserie.append(e[0])
                        labels.append('%s' % calibcritSMobslst_tmp[n])
                    else:
                        xserie.append(1 + xserie[-1])
                        yserie.append(e[0])
                        labels.append('%s' % calibcritSMobslst_tmp[n])
                n += 1
            offset = (max_tmp) * 0.05
            for i in range(len(xserie)):
                plt.scatter(xserie[i], yserie[i], marker = 'o', c='orange', s=15)
                if yserie[i] < max_tmp:
                    plt.text(xserie[i], yserie[i] + offset, '%.1f' % yserie[i], fontsize=6, ha='center', va='center')
            plt.xticks(xserie, labels)
            if ymin == None:
                tmp = np.min(np.ma.masked_where(np.asarray(yserie).flatten() == hnoflo, np.asarray(yserie).flatten()))
                if tmp > 0:
                    ymin_tmp = 0
                else:
                    ymin_tmp = 1.2 * tmp
                del tmp
            else:
                ymin_tmp = ymin
            plt.ylim(ymin_tmp, max_tmp)
            plt.setp(ax1.get_yticklabels(), fontsize=6)
            ax1.set_xlim(0, int(max(xserie)) + 1.0)
            labels = ax1.get_xticklabels(which='both')
            plt.setp(labels, rotation=90, fontsize=6)
            del labels

        calibcritHEADS_tmp = calibcritHEADS[p_mult:(p_mult + 1) * num_plt_max]
        calibcritHEADSc_tmp = calibcritHEADSc[p_mult:(p_mult + 1) * num_plt_max]
        calibcritHEADSobslst_tmp = calibcritHEADSobslst[p_mult:(p_mult + 1) * num_plt_max]
        calibcritHEADScobslst_tmp = calibcritHEADScobslst[p_mult:(p_mult + 1) * num_plt_max]
        if len(calibcritHEADS_tmp) > 0:
            test += 1
            ax2 = fig.add_subplot(2, 1, 2)
            ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            plt.ylabel('%s hydraulic heads %s' % (calibcrit, units[0]), fontsize=10, horizontalalignment='center')
            plt.grid(True)
            xserie = list(range(1, len(calibcritHEADSobslst_tmp) + 1))
            yserie_txt = list(itertools.chain.from_iterable(calibcritHEADS_tmp))
            yseriec_txt = list(itertools.chain.from_iterable(calibcritHEADSc_tmp))
            if calibcritHEADSmax == None:
                tmp = np.max(list(itertools.chain.from_iterable(calibcritHEADS_tmp)))
                if tmp > 0:
                    max_tmp = 1.2 * tmp
                else:
                    max_tmp = 1.0
                del tmp
            else:
                max_tmp = calibcritHEADSmax
            #        offset = (max_tmp)*0.05
            for i in range(len(xserie)):
                plt.scatter(xserie[i], calibcritHEADS_tmp[i], marker = 'o', c='orange', s=15)
                for ii, cc in enumerate(calibcritHEADScobslst_tmp):
                    if calibcritHEADSobslst_tmp[i] == cc:
                        plt.scatter(xserie[i], calibcritHEADSc_tmp[ii], marker = 'o', c='green', s=10)
                        if yseriec_txt[ii] < max_tmp:
                            plt.text(xserie[i] - .05, yseriec_txt[ii], '%.1f' % yseriec_txt[ii], fontsize=6, ha='right',
                                     va='center')
                if yserie_txt[i] < max_tmp:
                    plt.text(xserie[i] + .05, yserie_txt[i], '%.1f' % yserie_txt[i], fontsize=6, ha='left', va='center')
            plt.xticks(xserie, calibcritHEADSobslst_tmp)
            if ymin == None:
                tmp = np.min(np.ma.masked_where(np.asarray(calibcritHEADS_tmp).flatten() == hnoflo,
                                                np.asarray(calibcritHEADS_tmp).flatten()))
                if tmp < 0:
                    ymin_tmp = 1.2 * tmp
                else:
                    ymin_tmp = 0.8 * tmp
                del tmp
            else:
                ymin_tmp = ymin
            plt.ylim(ymin_tmp, max_tmp)
            plt.setp(ax2.get_yticklabels(), fontsize=6)
            ax2.set_xlim(0, len(calibcritHEADS_tmp) + 1)
            labels = ax2.get_xticklabels(which='both')
            plt.setp(labels, rotation=90, fontsize=6)
            del labels
        if test > 0:
            plt_export_fn_tmp = "%s_%d.png" % (plt_export_fn.split(".")[0], p)
            plt.savefig(plt_export_fn_tmp, dpi=150)
            plt.close('all')
        p_mult += num_plt_max

        # TODO export a table with calib. crit. results


##################
if __name__ == "__main__":
    print('\nWARNING!\nStart MARMITES-MODFLOW models using the script startMARMITES_v3.py\n')

# EOF
