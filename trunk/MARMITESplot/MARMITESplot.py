# -*- coding: cp1252 -*-

import os
import datetime
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def allPLOT(DateInput, P, PET, PE, Pe, SUST, Qs, Eu, Tu, Eg, Tg, S, Rp, R, Es, h_MF, h_SF, hmeas, Smeas, Sm, Sr, hnoflo, plot_export_fn, colors_nsl, hmax, hmin):
    """
    allGRAPH: GRAPH the computed data
    Use Matplotlib
    _______________________________________________________________________________

    INPUTS
            STATE VARIABLES
                TS              Time step
                P               Daily rainfall
                PET             Daily evapotranspiration
                Pe              Daily Excess rainfall
                Eu              Daily evaporation (bare soil)
                Tu              Daily transpiration
                S               Daily soil moisture
                Rp              Daily percolation
                SUST            Daily ponding
                Qs              Daily runoff
                R               Daily recharge
                h               Daily water level
                hmeas           Daily measured water level
    ______________________________________________________________________________
    ______________________________________________________________________________
    """

    monthsFmt=matplotlib.dates.DateFormatter('%y-%m-%d')
    lblspc = 0.05
    mkscale = 0.5

    fig = plt.figure(figsize=(2*11.7, 2*8.27))

    nsl = len(Tu)
    lbl_S = []
    lbl_S.append('Smeas')
    lbl_Rp = []
    lbl_Eu = []
    lbl_Tu = []
    lbl_Rp.append('R')
    lbl_Eu.append('PE')
    lbl_Eu.append('E_tot')
    lbl_Eu.append('Eu_tot')
    lbl_Tu.append('PET')
    lbl_Tu.append('T_tot')
    lbl_Tu.append('Tu_tot')
    for l in range(nsl):
        lbl_Rp.append('Rp_l'+str(l+1))
        lbl_S.append('Ssim_l'+str(l+1))
        lbl_Eu.append('Eu_l'+str(l+1))
        lbl_Tu.append('Tu_l'+str(l+1))
    lbl_Tu.append('Tg')
    lbl_Eu.append('Eg')

# First column of plots
    ax5=fig.add_subplot(4,2,7)
    plt.setp( ax5.get_xticklabels(), fontsize=8)
    plt.setp( ax5.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput, Smeas, 'mo', markersize=2, c = 'lime', markeredgecolor='lime')
    for l, (y, color, lbl) in enumerate(zip(S, colors_nsl, lbl_S[1:len(lbl_S)])) :
        ax5.plot_date(DateInput, y, '-', color=color, label=lbl)
    # x axis
    ax5.xaxis.set_major_formatter(monthsFmt)
    plt.xlim(DateInput[0],DateInput[len(P)-1])
    labels=ax5.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    plt.xlabel('Date', fontsize=10)
    # y axis
    ybuffer=0.1*(max(Sm)-min(Sr))
    plt.ylim(min(Sr) - ybuffer,max(Sm) + ybuffer)
    plt.ylabel('%', fontsize=10)
    ax5.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.2f'))
    # legend
    plt.legend(lbl_S, loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.grid(True)

    Eu_tot = np.zeros([len(Eu[0])], dtype = float)
    for l in range(len(Eu)):
        Eu_tot = Eu_tot + Eu[l]
    E_tot = Eu_tot + Eg
    ax1=fig.add_subplot(4,2,3, sharex=ax5)
    plt.setp( ax1.get_xticklabels(), visible=False)
    plt.setp( ax1.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,PE,'b-', color='lightblue', linewidth=2)
    plt.plot_date(DateInput,E_tot,'b-', color='darkblue', linewidth=1.5)
    plt.plot_date(DateInput,Eu_tot,'b-', color=colors_nsl[len(colors_nsl)-1])
    for l, (y, color, lbl) in enumerate(zip(Eu, colors_nsl, lbl_Eu[2:len(lbl_Eu)])):
        ax1.plot_date(DateInput, y, '--', color=color, label=lbl)
    plt.plot_date(DateInput,Eg,'b-', color='blue')
    plt.legend(lbl_Eu, loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    plt.grid(True)
    plt.ylabel('mm', fontsize=10)

    Tu_tot = np.zeros([len(Tu[0])], dtype = float)
    for l in range(len(Tu)):
        Tu_tot = Tu_tot + Tu[l]
    T_tot = Tu_tot + Tg
    ax3=fig.add_subplot(4,2,5, sharex=ax5)
    plt.setp( ax3.get_xticklabels(), visible=False)
    plt.setp( ax3.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,PET,'b-', color='lightblue', linewidth=2)
    plt.plot_date(DateInput,T_tot,'b-', color='darkblue',  linewidth=1.5)
    plt.plot_date(DateInput,Tu_tot,'b-', color=colors_nsl[len(colors_nsl)-1])
    for l, (y, color, lbl) in enumerate(zip(Tu, colors_nsl, lbl_Tu[2:len(lbl_Tu)])):
        ax3.plot_date(DateInput, y, '--', color=color, label=lbl)
    plt.plot_date(DateInput,Tg,'b-', color='blue')
    plt.legend(lbl_Tu, loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.grid(True)
    plt.ylabel('mm', fontsize=10)

    ax4=fig.add_subplot(4,2,1, sharex=ax5)
    plt.setp( ax4.get_xticklabels(), visible=False)
    plt.setp( ax4.get_yticklabels(), fontsize=8)
    ax4.bar(DateInput,P,color='darkblue', linewidth=0, align = 'edge', label='RF')
    ax4.bar(DateInput,Pe,color='deepskyblue', linewidth=0, align = 'edge', label='RFe')
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    plt.grid(True)
    plt.ylabel('mm', fontsize=10)

# Second column of plots
    ax2=fig.add_subplot(4,2,2, sharex=ax5)
    plt.setp( ax2.get_xticklabels(), visible=False)
    plt.setp( ax2.get_yticklabels(), fontsize=8)
    ax2.bar(DateInput,P,color='darkblue', linewidth=0, align = 'edge', label='RF')
    ax2.bar(DateInput,Pe,color='deepskyblue', linewidth=0, align = 'edge', label='RFe')
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    plt.grid(True)
    plt.ylabel('mm', fontsize=10)

    ax7=fig.add_subplot(4,2,4, sharex=ax5)
    plt.setp(ax7.get_xticklabels(), visible=False)
    plt.setp(ax7.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,Qs,'r--', c='lightblue', linewidth=2)
    plt.plot_date(DateInput,Es,'r-', c='darkblue', linewidth=1)
    plt.bar(DateInput, SUST, color='blue', linewidth=0, align = 'edge')
    ax7.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.1f'))
    labels=ax7.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
    plt.ylabel('mm', fontsize=10)
    plt.legend(['Qs','Es','SUST'], loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    plt.grid(True)

    ax6=fig.add_subplot(4,2,6, sharex=ax5)
    plt.setp( ax6.get_xticklabels(), visible=False)
    plt.setp( ax6.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,R,'-', c='blue', linewidth=2)
    for l, (y, color, lbl) in enumerate(zip(Rp, colors_nsl, lbl_Rp[1:len(lbl_Rp)])) :
        ax6.plot_date(DateInput, y, '--', color=color, label=lbl)
    plt.legend(lbl_Rp, loc = 0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.grid(True)
    plt.ylabel('mm', fontsize=10)

    ax8=fig.add_subplot(4,2,8, sharex=ax5)
    plt.setp( ax8.get_xticklabels(), fontsize=8)
    plt.setp( ax8.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,hmeas,'ro', markersize=2, c='lime', markeredgecolor='lime', )
    plt.plot_date(DateInput,h_MF,'-', color = 'b')
    plt.plot_date(DateInput,h_SF,'-', color = 'r')
    ax8.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
    labels=ax8.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
    ax8.xaxis.set_major_formatter(monthsFmt)
    labels=ax8.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    ybuffer=0.1*(hmax-hmin)
    if ybuffer == 0.0:
        ybuffer = 1.0
    plt.ylim((hmin - ybuffer, hmax + ybuffer))
    plt.ylabel('m', fontsize=10)
    plt.legend((r'hmeas',r'hsim_MF',r'hsim_SF'), loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.xlabel('Date', fontsize=10)
    plt.grid(True)

    plt.subplots_adjust(left=0.05, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig(plot_export_fn,dpi=150)
#    plt.show()
    plt.close()
    del fig

def plotMBerror(DateInput, MB, FLOOD, SATpart, POND, Runoff, plot_export_fn):
    """
    Plot mass balance error
    """
    monthsFmt = matplotlib.dates.DateFormatter('%y-%m-%d')
    lblspc = 0.05
    mkscale = 0.5

    fig = plt.figure(figsize=(11.7, 8.27))

    ax5=fig.add_subplot(5,1,5)
    plt.setp(ax5.get_xticklabels(), fontsize=8)
    plt.setp(ax5.get_yticklabels(), fontsize=8)
    ax5.bar(DateInput,Runoff,color='r', linewidth=0, align = 'edge')
    plt.xlabel('Date', fontsize=10)
    ax5.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1d'))
    ax5.yaxis.set_ticks((0,1))
    plt.legend(['Runoff'], loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.ylim(0,1)
    plt.ylabel('Occurence', fontsize=10)
    plt.setp(ltext, fontsize=8 )
    plt.grid(True)

    ax4=fig.add_subplot(5,1,4, sharex = ax5)
    plt.setp( ax4.get_xticklabels(), visible=False)
    plt.setp( ax4.get_yticklabels(), fontsize=8)
    ax4.bar(DateInput,POND,color='r', linewidth=0, align = 'edge')
    ax4.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1d'))
    ax4.yaxis.set_ticks((0,1))
    plt.legend(['POND'], loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.ylim(0,1)
    plt.ylabel('Occurence', fontsize=10)
    plt.grid(True)

    ax3=fig.add_subplot(5,1,3, sharex = ax5)
    plt.setp( ax3.get_xticklabels(), visible=False)
    plt.setp( ax3.get_yticklabels(), fontsize=8)
    ax3.bar(DateInput,FLOOD, color='r', linewidth=0, align = 'edge')
    ax3.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1d'))
    ax3.yaxis.set_ticks((0,1))
    plt.legend(['FLOOD'], loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.ylim(0,1)
    plt.ylabel('Occurence', fontsize=10)
    plt.grid(True)

    ax2=fig.add_subplot(5,1,2, sharex = ax5)
    plt.setp( ax2.get_xticklabels(), visible=False)
    plt.setp( ax2.get_yticklabels(), fontsize=8)
    ax2.bar(DateInput,SATpart,color='r', linewidth=0, align = 'edge')
    ax2.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1d'))
    ax2.yaxis.set_ticks((0,1))
    plt.legend(['SATpart'], loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.ylim(0,1)
    plt.ylabel('Occurence', fontsize=10)
    plt.grid(True)

    ax1=fig.add_subplot(5,1,1, sharex = ax5)
    plt.setp( ax1.get_xticklabels(), visible=False)
    plt.setp( ax1.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,MB,'-', c='r')
    # y axis
    plt.ylabel('mm', fontsize=10)
    # legend
    plt.legend(['MB'], loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.grid(True)

    # x axes
    ax5.xaxis.set_major_formatter(monthsFmt)
    plt.xlim(DateInput[0],DateInput[len(MB)-1])
    labels=ax5.get_xticklabels()
    plt.setp(labels, 'rotation', 90)

    plt.subplots_adjust(left=0.05, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig(plot_export_fn,dpi=150)
#    plt.show()
    plt.close()
    del fig

def plotLAYER(TS, ncol, nrow, nlay, nplot, V, cmap, CBlabel, msg, plttitle, MM_ws, interval_type = 'arange', interval_diff = 1, interval_num = 1, Vmax = 0, Vmin = 0, fmt = '%.2f'):

    # Store some arrays for plotting
    x = np.arange(0.5, ncol+1.5, 1)
    y = np.arange(0.5, nrow+1.5, 1)
    xg,yg = np.meshgrid(x,y)

    x = np.arange(1, ncol+1, 1)
    y = np.arange(1, nrow+1, 1)
    xg1,yg1 = np.meshgrid(x,y)

    ax = []
    fig = plt.figure()
    for L in range(nplot):
        if interval_type == 'arange':
            ticks = np.arange(Vmin,Vmax,interval_diff)
        elif interval_type == 'linspace':
            ticks = np.linspace(Vmin,Vmax,interval_num)
        ax.append(fig.add_subplot(1,nlay,L+1, axisbg='silver'))
        plt.setp(ax[L].get_xticklabels(), fontsize=8)
        plt.setp(ax[L].get_yticklabels(), fontsize=8)
        plt.ylabel('row i', fontsize=10)
        plt.grid(True)
        plt.xlabel('col j', fontsize=10)
        ax[L].xaxis.set_ticks(np.arange(1,ncol+1))
        ax[L].yaxis.set_ticks(np.arange(1,nrow+1))
        if Vmax>Vmin:
            PC = plt.pcolor(xg, yg, V[L], cmap = cmap, vmin = Vmin, vmax = Vmax)
            CS = plt.contour(xg1, yg1[::-1], V[L][::-1],ticks, colors = 'gray')
            plt.clabel(CS, inline=1, fontsize=8, fmt=fmt, colors = 'gray')
            if L==nplot-1:
                CB = plt.colorbar(PC, shrink=0.8, extend='both', ticks = ticks, format = fmt)
                CB.set_label(CBlabel, fontsize = 8)
                plt.setp(CB.ax.get_yticklabels(), fontsize=8)
            plt.title('layer ' + str(L+1)+', time step ' + str(TS+1), fontsize = 10)
        else:
            plt.title('layer ' + str(L+1)+', time step ' + str(TS+1) + ': ' + msg, fontsize = 10)
        plt.ylim(plt.ylim()[::-1])
        plt.axis('scaled')
    plot_export_fn = os.path.join(MM_ws, plttitle + '_TS'+str(TS+1) + '.png')
    plt.savefig(plot_export_fn)
    plt.close()
    del fig
    del ax

##    #Make a cross-sectional figure of layers 1, 2, and 10
##    plt.figure()
##    plt.plot(xg[0,:],valg[:,50,0],label='Top layer')
##    plt.plot(xg[0,:],valg[:,50,1],label='Second layer')
##    plt.legend(loc='best')