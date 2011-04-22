# -*- coding: utf-8 -*-

import os
import datetime
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def allPLOT(DateInput, P, PET, PE, Pe, dPOND, POND, Ro, Eu, Tu, Eg, Tg, S, dS, Spc, Rp, SEEPAGE, R, Rn, Es, MB, h_MF, h_SF, hmeas, Smeas, Sm, Sr, hnoflo, plot_export_fn, colors_nsl, hmax, hmin):
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
                POND            Daily ponding
                Ro              Daily runoff
                R               Daily recharge
                h               Daily water level
                hmeas           Daily measured water level
    ______________________________________________________________________________
    ______________________________________________________________________________
    """

    monthsFmt=matplotlib.dates.DateFormatter('%y-%m-%d')
    lblspc = 0.05
    mkscale = 0.5

    fig = plt.figure(num=None, figsize=(2*8.27, 2*11.7), dpi=30, facecolor='gray', edgecolor='r', linewidth=1.0)
    #(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    #(figsize=(25, 15), edgecolor = 'black')  #(2*8.27, 2*11.7))

    nsl = len(Tu)
    lbl_Spc = []
    lbl_Spc.append('Smeas')
    lbl_S = []
    lbl_dS = []
    lbl_Rp = []
    lbl_Eu = []
    lbl_Tu = []
    lbl_Rp.append('R')
    lbl_Rp.append('Rn')
    lbl_Eu.append('PE')
    lbl_Eu.append('E_tot')
    lbl_Eu.append('Eu_tot')
    lbl_Tu.append('PET')
    lbl_Tu.append('T_tot')
    lbl_Tu.append('Tu_tot')
    for l in range(nsl):
        lbl_Spc.append('S_l'+str(l+1))
        lbl_S.append('S_l'+str(l+1))
        lbl_dS.append('dS_l'+str(l+1))
        lbl_Eu.append('Eu_l'+str(l+1))
        lbl_Tu.append('Tu_l'+str(l+1))
    for l in range(nsl-1):
        lbl_Rp.append('Rp_l'+str(l+1))
    lbl_Rp.append('SEEPAGE')
    lbl_Tu.append('Tg')
    lbl_Eu.append('Eg')

    ax1=fig.add_subplot(10,1,1)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), fontsize=8)
    ax1.bar(DateInput,P,color='darkblue', linewidth=0, align = 'center', label='RF')
    ax1.bar(DateInput,Pe,color='deepskyblue', linewidth=0, align = 'center', label='RFe')
    plt.xlim(DateInput[0]-1,DateInput[len(P)-1]+1)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    plt.grid(True)
    plt.ylabel('mm', fontsize=10)
    ax1.xaxis.set_major_formatter(monthsFmt)
    ax1.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.1f'))

    ax2=fig.add_subplot(10,1,2, sharex=ax1)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,Ro,'r-', c='darkblue', linewidth=2, label = 'Ro')
    plt.plot_date(DateInput,Es,'r-', c='darkblue', linewidth=1, label = 'Es')
    plt.bar(DateInput, POND, color='lightblue', linewidth=0, align = 'center', label = 'POND')
    plt.bar(DateInput, dPOND, color='blue', width=0.60, linewidth=0, align = 'center', label = 'dPOND')
    ax2.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.1f'))
    plt.xlim(DateInput[0]-1,DateInput[len(P)-1]+1)
    plt.ylabel('mm', fontsize=10)
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    plt.grid(True)
    ax2.xaxis.set_major_formatter(monthsFmt)

    Eu_tot = np.zeros([len(Eu[0])], dtype = float)
    for l in range(len(Eu)):
        Eu_tot = Eu_tot + Eu[l]
    E_tot = Eu_tot + Eg
    ax3=fig.add_subplot(10,1,3, sharex=ax1)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,PE,'-', color='lightblue', linewidth=3)
    plt.plot_date(DateInput,E_tot,'-', color='darkblue', linewidth=1)
    plt.plot_date(DateInput,Eu_tot,'-.', color=colors_nsl[len(colors_nsl)-1])
    for l, (y, color, lbl) in enumerate(zip(Eu, colors_nsl, lbl_Eu[2:len(lbl_Eu)])):
        ax3.plot_date(DateInput, y, '-', color=color, label=lbl)
    plt.plot_date(DateInput,Eg,'b-', color='blue')
    plt.xlim(DateInput[0]-1,DateInput[len(P)-1]+1)
    plt.legend(lbl_Eu, loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    plt.grid(True)
    plt.ylabel('mm', fontsize=10)
    ax3.xaxis.set_major_formatter(monthsFmt)
    ax3.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.1f'))

    Tu_tot = np.zeros([len(Tu[0])], dtype = float)
    for l in range(len(Tu)):
        Tu_tot = Tu_tot + Tu[l]
    T_tot = Tu_tot + Tg
    ax4=fig.add_subplot(10,1,4, sharex=ax1)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,PET,'b-', color='lightblue', linewidth=3)
    plt.plot_date(DateInput,T_tot,'b-', color='darkblue',  linewidth=1)
    plt.plot_date(DateInput,Tu_tot,'-.', color=colors_nsl[len(colors_nsl)-1])
    for l, (y, color, lbl) in enumerate(zip(Tu, colors_nsl, lbl_Tu[2:len(lbl_Tu)])):
        ax4.plot_date(DateInput, y, '-', color=color, label=lbl)
    plt.plot_date(DateInput,Tg,'b-', color='blue')
    plt.xlim(DateInput[0]-1,DateInput[len(P)-1]+1)
    plt.legend(lbl_Tu, loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.grid(True)
    plt.ylabel('mm', fontsize=10)
    ax4.xaxis.set_major_formatter(monthsFmt)
    ax4.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.1f'))

    ax5a=fig.add_subplot(20,1,9, sharex=ax1)
    plt.setp(ax5a.get_xticklabels(), visible=False)
    plt.setp(ax5a.get_yticklabels(), fontsize=8)
    for l, (y, color, lbl) in enumerate(zip(dS, colors_nsl, lbl_dS)) :
        ax5a.plot_date(DateInput, y, '-', color=color, label=lbl)
    # x axis
    plt.xlim(DateInput[0]-1,DateInput[len(P)-1]+1)
    # y axis
    plt.ylabel('mm', fontsize=10)
    ax5a.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.1f'))
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.grid(True)
    ax5a.xaxis.set_major_formatter(monthsFmt)

    ax5b=fig.add_subplot(20,1,10, sharex=ax1)
    plt.setp(ax5b.get_xticklabels(), visible=False)
    plt.setp(ax5b.get_yticklabels(), fontsize=8)
    for l, (y, color, lbl) in enumerate(zip(S, colors_nsl, lbl_S)) :
        ax5b.plot_date(DateInput, y, '-', color=color, label=lbl)
    # x axis
    plt.xlim(DateInput[0]-1,DateInput[len(P)-1]+1)
    # y axis
    plt.ylabel('mm', fontsize=10)
    ax5b.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.1f'))
    # legend
    plt.legend(loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.grid(True)
    ax5b.xaxis.set_major_formatter(monthsFmt)

    ax6=fig.add_subplot(10,1,6, sharex=ax1)
    plt.setp(ax6.get_xticklabels(), visible=False)
    plt.setp(ax6.get_yticklabels(), fontsize=8)
    plt.bar(DateInput,SEEPAGE, color='lightblue', linewidth=0, align = 'center', label='SEEPAGE')
    plt.plot_date(DateInput,R,'-', c='darkblue', linewidth=2)
    plt.plot_date(DateInput,Rn,'-', c='blue', linewidth=1)
    for l, (y, color, lbl) in enumerate(zip(Rp, colors_nsl, lbl_Rp[2:len(lbl_Rp)])) :
        ax6.plot_date(DateInput, y, '-', color=color, label=lbl)
    plt.xlim(DateInput[0]-1,DateInput[len(P)-1]+1)
    plt.legend(lbl_Rp, loc = 0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.grid(True)
    plt.ylabel('mm', fontsize=10)
    ax6.xaxis.set_major_formatter(monthsFmt)
    ax6.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.1f'))

    ax7=fig.add_subplot(10,1,7, sharex=ax1)
    plt.setp(ax7.get_xticklabels(), visible=False)
    plt.setp(ax7.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput, Smeas, 'mo', markersize=2, c = 'lime', markeredgecolor='lime')
    for l, (y, color, lbl) in enumerate(zip(Spc, colors_nsl, lbl_Spc[1:len(lbl_Spc)])) :
        ax7.plot_date(DateInput, y, '-', color=color, label=lbl)
    # x axis
    plt.xlim(DateInput[0]-1,DateInput[len(P)-1]+1)
    # y axis
    ybuffer=0.1*(max(Sm)-min(Sr))
    plt.ylim(min(Sr) - ybuffer,max(Sm) + ybuffer)
    plt.ylabel('%', fontsize=10)
    ax7.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.2f'))
    # legend
    plt.legend(lbl_Spc, loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.grid(True)
    ax7.xaxis.set_major_formatter(monthsFmt)

    ax8=fig.add_subplot(10,1,8, sharex=ax1)
    plt.setp(ax8.get_xticklabels(), fontsize=8)
    plt.setp(ax8.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,hmeas,'ro', markersize=2, c='lime', markeredgecolor='lime', )
    plt.plot_date(DateInput,h_MF,'-', color = 'b')
    plt.plot_date(DateInput,h_SF,'-', color = 'r')
    plt.xlim(DateInput[0]-1,DateInput[len(P)-1]+1)
    ax8.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
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

    ax10=fig.add_subplot(10,1,10, sharex=ax1)
    plt.setp(ax10.get_xticklabels(), visible=False)
    plt.setp(ax10.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,MB,'-', c='r')
    # y axis
    plt.ylabel('mm', fontsize=10)
    plt.grid(True)
    plt.xlim(DateInput[0]-1,DateInput[len(MB)-1]+1)
    if MB.max()<0.001 and MB.min()>-0.001:
        plt.ylim(-0.1,0.1)
    # legend
    plt.legend(['MB'], loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    ax10.xaxis.set_major_formatter(monthsFmt)
    ax10.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2g'))

    plt.subplots_adjust(left=0.05, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig(plot_export_fn,dpi=150)
#    plt.show()
    plt.close()
    del fig

##################

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
        if np.nanmax(V[L])>np.nanmin(V[L]):
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
    plot_export_fn = os.path.join(MM_ws, '00_' + plttitle + '_TS'+str(TS+1) + '.png')
    plt.savefig(plot_export_fn)
    plt.close()
    del fig
    del ax

##    #Make a cross-sectional figure of layers 1, 2, and 10
##    plt.figure()
##    plt.plot(xg[0,:],valg[:,50,0],label='Top layer')
##    plt.plot(xg[0,:],valg[:,50,1],label='Second layer')
##    plt.legend(loc='best')

##################

def plotGWbudget(flxlst, flxlbl, colors_flx, plot_export_fn, fluxmax, fluxmin):
    """
    Plot GW budget
    """

    def autolabel(rects):
    # attach some text labels
        for rect in rects:
            y = rect.get_y()
            height = rect.get_height()
            if y>=0:
                ytext = height + 0.1
                va = 'bottom'
            else:
                ytext = -height - 0.1
                height = - height
                va = 'top'
            ax1.text(rect.get_x()+rect.get_width()/2., ytext, '%.2f'%float(height),
                    ha='center', va=va)
    ##################

    lblspc = 0.05
    mkscale = 0.5

    fig = plt.figure(figsize=(2*11.7, 2*8.27))

    ax1=fig.add_subplot(2,1,1)
    plt.setp( ax1.get_xticklabels(), fontsize=10)
    plt.setp( ax1.get_yticklabels(), fontsize=10)
    x = np.arange(len(flxlst))
    width = 0.8
    rects = plt.bar(x , flxlst, color=colors_flx, linewidth=0.5, width=width, align = 'edge', label=flxlbl)
    # y axis
    plt.ylabel('mm/d', fontsize=10)
    plt.grid(True)
    # fmt xaxis
    ax1.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%1.2f'))
    ax1.set_xticks(x+width/2)
    ax1.set_xticklabels(flxlbl)
    labels=ax1.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    autolabel(rects)
    plt.xlim(0,len(flxlst))
    plt.ylim(fluxmin, fluxmax)
    ax1.axhline(0, color='black', lw=1)

    plt.subplots_adjust(left=0.05, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig(plot_export_fn,dpi=150)
#    plt.show()
    plt.close()
    del fig

