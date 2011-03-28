# -*- coding: cp1252 -*-

import numpy as np
from matplotlib.ticker import FormatStrFormatter
from matplotlib.dates import DateFormatter
import matplotlib.pyplot as plt

def allPLOT(DateInput, P, PET, PE, Pe, SUST, Qs, Eu, Tu, S, Rp, R, Es, h_MF, h_SF, hmeas, Smeas, Sm, Sr, hnoflo, plot_export_fn, colors_nsl):
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

    monthsFmt=DateFormatter('%y-%m-%d')
    lblspc = 0.05
    mkscale = 0.5

#__________________Create outputs plots______________________#
#    ioff()


    fig = plt.figure(figsize=(2*11.7, 2*8.27))

    nsl = len(Tu)
    lbl_S = []
    lbl_S.append('Smeas')
    lbl_Rp = []
    lbl_Eu = []
    lbl_Tu = []
    lbl_Rp.append('R')
    lbl_Eu.append('PE')
    lbl_Eu.append('Eu_tot')
    lbl_Tu.append('PET')
    lbl_Tu.append('Tu_tot')
    for l in range(nsl):
        lbl_Rp.append('Rp_l'+str(l+1))
        lbl_S.append('Ssim_l'+str(l+1))
        lbl_Eu.append('Eu_l'+str(l+1))
        lbl_Tu.append('Tu_l'+str(l+1))

# First column of plots
    ax5=fig.add_subplot(4,2,7)
    plt.setp( ax5.get_xticklabels(), fontsize=8)
    plt.setp( ax5.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput, Smeas, 'mo', markersize=2, c = 'lime', markeredgecolor='lime')
    for l, (y, color, lbl) in enumerate(zip(S, colors_nsl, lbl_S[1:len(lbl_S)])) :
        ax5.plot_date(DateInput, y, '-', color=color, label=lbl)
    ybuffer=0.1*(max(Sm)-min(Sr))
    plt.ylim(min(Sr) - ybuffer,max(Sm) + ybuffer)
    plt.ylabel('%')
    ax5.yaxis.set_major_formatter(FormatStrFormatter('%1.2f'))
    ax5.xaxis.set_major_formatter(monthsFmt)
    labels=ax5.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    plt.xlabel('Date', fontsize=10)
    plt.legend(lbl_S, loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.grid(True)

    Eu_tot = np.zeros([len(Eu[0])], dtype = float)
    for l in range(len(Eu)):
        Eu_tot = Eu_tot + Eu[l]
    ax1=fig.add_subplot(4,2,3, sharex=ax5)
    plt.setp( ax1.get_xticklabels(), visible=False)
    plt.setp( ax1.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,PE,'b--', color='lightblue', linewidth=2)
    plt.plot_date(DateInput,Eu_tot,'b--', color='blue')
    for l, (y, color, lbl) in enumerate(zip(Eu, colors_nsl, lbl_Eu[2:len(lbl_Eu)])):
        ax1.plot_date(DateInput, y, '-', color=color, label=lbl)
    plt.legend(lbl_Eu, loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8)
    plt.grid(True)
    plt.ylabel('mm', fontsize=10)

    Tu_tot = np.zeros([len(Tu[0])], dtype = float)
    for l in range(len(Tu)):
        Tu_tot = Tu_tot + Tu[l]
    ax3=fig.add_subplot(4,2,5, sharex=ax5)
    plt.setp( ax3.get_xticklabels(), visible=False)
    plt.setp( ax3.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,PET,'b--', color='lightblue', linewidth=2)
    plt.plot_date(DateInput,Tu_tot,'b--', color='blue')
    for l, (y, color, lbl) in enumerate(zip(Tu, colors_nsl, lbl_Tu[2:len(lbl_Tu)])):
        ax3.plot_date(DateInput, y, '-', color=color, label=lbl)
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
    labels=ax4.get_xticklabels()
    plt.setp(labels, 'horizontalalignment', 'right')


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
    labels=ax2.get_xticklabels()
    plt.setp(labels, 'horizontalalignment', 'right')

    ax7=fig.add_subplot(4,2,4, sharex=ax5)
    plt.setp(ax7.get_xticklabels(), visible=False)
    plt.setp(ax7.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,Qs,'r--', c='lightblue', linewidth=2)
    plt.plot_date(DateInput,Es,'r-', c='darkblue', linewidth=1)
    plt.bar(DateInput, SUST, color='blue', linewidth=0, align = 'edge')
    ax7.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))
    labels=ax7.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
    plt.xlim((DateInput[0],DateInput[len(P)-1]))
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
    ax8.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    labels=ax8.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
    ax8.xaxis.set_major_formatter(monthsFmt)
#    ax2.autoscale_view()
    labels=ax8.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    plt.xlim((DateInput[0],DateInput[len(P)-1]))
    hmax=[]
    hmin=[]
    if h_MF[0]==hnoflo:
        hmax0=-999
        hmin0=9999
    else:
        hmax0=hmin0=h_MF[0]
    if hmeas[0]==hnoflo:
        hmax1=-999
        hmin1=9999
    else:
        hmax1=hmin1=hmeas[0]
    if h_SF[0]==hnoflo:
        hmax2=-999
        hmin2=9999
    else:
        hmax2=hmin2=h_SF[0]

    for i in range(0,len(DateInput)):
        if h_MF[i]>hmax0:
            hmax0=h_MF[i]
        elif h_MF[i]<hmin0:
            hmin0=h_MF[i]
        hmax.append(hmax0)
        hmin.append(hmin0)

        if hmeas[i]>hmax1:
            hmax1=hmeas[i]
        elif hmeas[i]<hmin1:
            hmin1=hmeas[i]
        hmax.append(hmax1)
        hmin.append(hmin1)

        if h_SF[i]>hmax2:
            hmax2=h_SF[i]
        elif h_MF[i]<hmin2:
            hmin2=h_SF[i]
        hmax.append(hmax2)
        hmin.append(hmin2)

    hmax = max(hmax)
    hmin = min(hmin)
    ybuffer=0.1*(hmax-hmin)
    if ybuffer == 0.0:
        ybuffer = 1.0
    plt.ylim((hmin - ybuffer, hmax + ybuffer))
    plt.ylabel('m')
    plt.legend((r'hmeas',r'hsim_MF',r'hsim_SF'), loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.xlabel('Date', fontsize=10)
    plt.grid(True)

    plt.subplots_adjust(left=0.05, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig(plot_export_fn,dpi=150)
#    plt.show()


def plotMBerror(DateInput, MB, plot_export_fn):
    """
    Plot mass balance error
    """
    monthsFmt=DateFormatter('%y-%m-%d')
    lblspc = 0.05
    mkscale = 0.5
    fig = plt.figure(figsize=(11.7, 8.27))

    ax1=fig.add_subplot(1,1,1)
    plt.setp( ax1.get_xticklabels(), visible=False)
    plt.setp( ax1.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,MB,'-', c='r')
    plt.legend(['MB'], loc=0, labelspacing=lblspc, markerscale=mkscale)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()
    plt.setp(ltext, fontsize=8 )
    plt.ylabel('mm', fontsize=10)
    plt.grid(True)

    plt.subplots_adjust(left=0.05, bottom=0.10, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig(plot_export_fn,dpi=150)
#    plt.show()