# -*- coding: cp1252 -*-

from pylab import *
from matplotlib.dates import MonthLocator, DateFormatter
from matplotlib.ticker import FormatStrFormatter
from sys import *
import os
import datetime



##=========================================================================================##
##===========================| PIEZO GRAPHS |========================================##
##=========================================================================================##


def piezocalibGRAPH(DateInput, h, hmeas, P, Pe):
    """
    calibGRAPH: GRAPH the computed data and the piezometric calibration one
    Use Matplotlib
    _______________________________________________________________________________

    INPUTS
            STATE VARIABLES
                TS              Time step
                h               Daily water level
                hmeas           Daily measured water level
    ______________________________________________________________________________
    ______________________________________________________________________________
    """

    months=MonthLocator()
    monthsFmt=DateFormatter('%y-%m')

#__________________Create outputs plots______________________#
#    ioff()
    figCalib=figure()
    figCalib.Title='Calibration graphs'

    ax2=subplot(111)
    setp( ax2.get_xticklabels(), fontsize=8)
    setp( ax2.get_yticklabels(), fontsize=8)
    plot_date(DateInput,hmeas,'ro', markersize=2, markeredgecolor='r')
    plot_date(DateInput,h,'-')
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))
    labels=ax2.get_yticklabels()
    setp(labels, 'rotation', 90)
    ax2.xaxis.set_major_locator(months)
    ax2.xaxis.set_major_formatter(monthsFmt)
    xlim((DateInput[0],DateInput[len(h)-1]))
    if h[0]==-999:
        hmax=-999
        hmin=9999
    else:
        hmax=hmin=h[0]
    labels=ax2.get_xticklabels()
    setp(labels, 'rotation', 90)
    for i in range(0,len(DateInput)):
        if h[i]>hmax:
            hmax=h[i]
        elif h[i]<hmin:
            hmin=h[i]
    if hmeas[0]==-999:
        hmmax=-999
        hmmin=9999
    else:
        hmmax=hmmin=hmeas[0]
    for i in range(0,len(DateInput)):
        if hmeas[i]!=-999:
            if hmeas[i]>hmmax:
                hmmax=hmeas[i]
            elif hmeas[i]<hmmin:
                hmmin=hmeas[i]
    if hmin>hmmin:
        hmin=hmmin
    if hmax<hmmax:
        hmax=hmmax
    ybuffer=0.1*(hmax-hmin)
    ylim((hmin - ybuffer, hmax + ybuffer))
    ylabel('m')
    legend((r'h obs',r'h sim'), loc=0)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
    xlabel(r'Date')
#    ax2.autoscale_view()
    grid(True)

##    ax3=subplot(211, sharex=ax2)
##    setp( ax3.get_xticklabels(), visible=False)
##    setp( ax3.get_yticklabels(), fontsize=8)
###    DateInput1=range(0,len(DateInput))
##    bar(DateInput,P,color='b', linewidth=0, align = 'edge')
##    bar(DateInput,Pe,color='m', linewidth=0, align = 'edge')
##    legend((r'P', r'Pe'), loc=0)
##    leg = gca().get_legend()
##    ltext  = leg.get_texts()  # all the text.Text instance in the legend
##    setp(ltext, fontsize='small')    # the legend text fontsize
##    setp(labels, 'rotation', 90)
##    grid(True)
##    ylabel('mm')

    subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.95, wspace=0.1, hspace=0.05)
#    draw()
#    ion()
    show()
    del DateInput, h, hmeas, P, Pe



##=========================================================================================##
##===========================| CALIBRATION GRAPHS |========================================##
##=========================================================================================##


def calibGRAPH(DateInput, P, PET, Pe, ETa, S, R, h, hmeas, Smeas, Sm, Sr):
    """
    calibGRAPH: GRAPH the computed data and the calibration one, that it h and S
    Use Matplotlib
    _______________________________________________________________________________

    INPUTS
            STATE VARIABLES
                TS              Time step
                S               Daily soil moisture
                Smeas           Daily measured soil moisture
                h               Daily water level
                hmeas           Daily measured water level
    ______________________________________________________________________________
    ______________________________________________________________________________
    """

    months=MonthLocator()
    monthsFmt=DateFormatter('%y-%m')

#__________________Create outputs plots______________________#
#    ioff()
    figCalib=figure()
    figCalib.Title='Calibration graphs'

    ax5=subplot(515)
    setp( ax5.get_xticklabels(), fontsize=8)
    setp( ax5.get_yticklabels(), fontsize=8)
    plot_date(DateInput,hmeas,'ro', markersize=2, markeredgecolor='r')
    plot_date(DateInput,h,'-')
    ax5.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))
    labels=ax5.get_yticklabels()
#    setp(labels, 'rotation', 90)
##    ax5.xaxis.set_major_locator(months)
##    ax5.xaxis.set_major_formatter(monthsFmt)
##    xlim((DateInput[0],DateInput[len(h)-1]))
    if h[0]==-999:
        hmax=-999
        hmin=9999
    else:
        hmax=hmin=h[0]
##    labels=ax5.get_xticklabels()
##    setp(labels, 'rotation', 90)
    for i in range(0,len(DateInput)):
        if h[i]>hmax:
            hmax=h[i]
        elif h[i]<hmin:
            hmin=h[i]
    if hmeas[0]==-999:
        hmmax=-999
        hmmin=9999
    else:
        hmmax=hmmin=hmeas[0]
    for i in range(0,len(DateInput)):
        if hmeas[i]!=-999:
            if hmeas[i]>hmmax:
                hmmax=hmeas[i]
            elif hmeas[i]<hmmin:
                hmmin=hmeas[i]
    if hmin>hmmin:
        hmin=hmmin
    if hmax<hmmax:
        hmax=hmmax
    ybuffer=0.1*(hmax-hmin)
    ylim((hmin - ybuffer, hmax + ybuffer))
    ylabel('m')
    legend((r'h obs',r'h sim'), loc=0)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
    xlabel(r'Date')
#    ax2.autoscale_view()
    grid(True)

    ax4=subplot(514, sharex=ax5)
    setp( ax4.get_xticklabels(), visible=False)
    setp( ax4.get_yticklabels(), fontsize=8)
    grid(True)
    plot_date(DateInput,R,'-')
    legend('R', loc=0)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
    ylabel('mm')

    ax3=subplot(513, sharex=ax5)
    setp( ax3.get_xticklabels(), visible=False)
    setp( ax3.get_yticklabels(), fontsize=8)
    plot_date(DateInput, Smeas, 'ro', markersize=2, markeredgecolor='m')
    plot_date(DateInput,S, '-')
    legend((r'S obs',r'S sim'), loc=0)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
#    ybuffer=0.01*(float(Sm)-float(Sr))
#    ylim((float(Sr) - ybuffer,float(Sm) + ybuffer))
#    ylabel('mm')
#    ylim(0,1)
    ylim(float(Sr),float(Sm))
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
    labels=ax3.get_yticklabels()
 #   setp(labels, 'rotation', 90)
    grid(True)
    ylabel('%')

    ax2=subplot(512, sharex=ax5)
    setp( ax2.get_xticklabels(), visible=False)
    setp( ax2.get_yticklabels(), fontsize=8)
    plot_date(DateInput,PET,'b-')
    plot_date(DateInput,ETa,'r-')
    legend((r'PET',r'ETa'), loc=0)  #,  fontsize=10)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
    grid(True)
    ylabel('mm')

    ax1=subplot(511, sharex=ax5)
    setp( ax1.get_xticklabels(), visible=False)
    setp( ax1.get_yticklabels(), fontsize=8)
#    DateInput1=range(0,len(DateInput))
    bar(DateInput,P,color='b', linewidth=0, align = 'edge')
    bar(DateInput,Pe,color='r', linewidth=0, align = 'edge')
    legend(['P', 'Pe'], loc=0)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
 #   setp(labels, 'rotation', 90)
    grid(True)
    ylabel('mm')
    ax1.xaxis.set_major_locator(months)
    ax1.xaxis.set_major_formatter(monthsFmt)
    xlim((DateInput[0],DateInput[len(h)-1]))
    labels=ax5.get_xticklabels()
    setp(labels, 'rotation', 90)

    subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.95, wspace=0.1, hspace=0.05)
#    draw()
#    ion()
    show()
    del DateInput, P, PET, Pe, ETa, S, R, h, hmeas, Smeas, Sm, Sr


##=========================================================================================##
##===========================| ALL GRAPHS |====================================================##
##=========================================================================================##

def allGRAPH(DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, Sm, Sr):
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
                ETa             Daily evapotranspiration
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

    months=MonthLocator()
    monthsFmt=DateFormatter('%y-%m')

#__________________Create outputs plots______________________#
#    ioff()
    figCalib=figure()
    figCalib.Title='All graphs'

# First column of graphs

    ax7=subplot(427)
    setp(ax7.get_xticklabels(), fontsize=8)
    setp(ax7.get_yticklabels(), fontsize=8)
    plot_date(DateInput,Qs,'r-')
    ax7.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))
    labels=ax7.get_yticklabels()
    setp(labels, 'rotation', 90)
    ax7.xaxis.set_major_locator(months)
    ax7.xaxis.set_major_formatter(monthsFmt)
#    ax2.autoscale_view()
    labels=ax7.get_xticklabels()
    setp(labels, 'rotation', 90)
#   DateInput1=range(0,len(DateInput))
    bar(DateInput, SUST, linewidth=0, align = 'edge')
    xlim((DateInput[0],DateInput[len(S)-1]))
    ylabel('mm')
    legend(['Qs', 'SUST'], loc=0)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
    xlabel(r'Date')
#    xlim((0,len(P)))
    grid(True)

    ax5=subplot(425, sharex=ax7)
    setp( ax5.get_xticklabels(), visible=False)
    setp( ax5.get_yticklabels(), fontsize=8)
    plot_date(DateInput, Smeas, 'mo', markersize=2, markeredgecolor='m')
    plot_date(DateInput,S, '-')
    ybuffer=0.1*(float(Sm)-float(Sr))
    ylim((float(Sr) - ybuffer,float(Sm) + ybuffer))
    ylabel('%')
    ax5.yaxis.set_major_formatter(FormatStrFormatter('%1.1f'))
    labels=ax5.get_yticklabels()
    setp(labels, 'rotation', 90)
    legend((r'S obs',r'S sim'), loc=0)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
    grid(True)

    ax3=subplot(423, sharex=ax7)
    setp( ax3.get_xticklabels(), visible=False)
    setp( ax3.get_yticklabels(), fontsize=8)
    plot_date(DateInput,PET,'b-')
    plot_date(DateInput,ETa,'r-')
    legend((r'PET',r'ETa'), loc=0)  #,  fontsize=10)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
    grid(True)
    ylabel('mm')

    ax1=subplot(421, sharex=ax7)
    setp( ax1.get_xticklabels(), visible=False)
    setp( ax1.get_yticklabels(), fontsize=8)
#    DateInput1=range(0,len(DateInput))
    bar(DateInput,P,color='b', linewidth=0, align = 'edge')
    legend(['P'], loc=0)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
    grid(True)
    ylabel('mm')

# Second column of graphs
    ax2=subplot(422, sharex=ax7)
    setp( ax2.get_xticklabels(), visible=False)
    setp( ax2.get_yticklabels(), fontsize=8)
#    DateInput1=range(0,len(DateInput))
    bar(DateInput,Pe,color='m', linewidth=0, align = 'edge')
    legend(['Pe'], loc=0)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
    grid(True)
    ylabel('mm')
    labels=ax2.get_xticklabels()
    setp(labels, 'horizontalalignment', 'right')
    #setp(gca(), 'horizontalalignment', 'right')

    ax4=subplot(424, sharex=ax7)
    setp( ax4.get_xticklabels(), visible=False)
    setp( ax4.get_yticklabels(), fontsize=8)
#    DateInput1=range(0,len(DateInput))
    bar(DateInput,Rp,linewidth=0, align = 'edge')
    legend(['Rp'], loc=0)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
    grid(True)
    ylabel('mm')

    ax6=subplot(426, sharex=ax7)
    setp( ax6.get_xticklabels(), visible=False)
    setp( ax6.get_yticklabels(), fontsize=8)
    grid(True)
    plot_date(DateInput,R,'-')
    legend('R', loc=0)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
    ylabel('mm')

    ax8=subplot(428, sharex=ax7)
    setp( ax8.get_xticklabels(), fontsize=8)
    setp( ax8.get_yticklabels(), fontsize=8)
    plot_date(DateInput,hmeas,'ro', markersize=2, markeredgecolor='r')
    plot_date(DateInput,h,'-')
    ax8.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    labels=ax8.get_yticklabels()
    setp(labels, 'rotation', 90)
    ax8.xaxis.set_major_locator(months)
    ax8.xaxis.set_major_formatter(monthsFmt)
#    ax2.autoscale_view()
    labels=ax8.get_xticklabels()
    setp(labels, 'rotation', 90)
    xlim((DateInput[0],DateInput[len(S)-1]))
    if h[0]==-999:
        hmax=-999
        hmin=9999
    else:
        hmax=hmin=h[0]
    for i in range(0,len(DateInput)):
        if h[i]>hmax:
            hmax=h[i]
        elif h[i]<hmin:
            hmin=h[i]
    if hmeas[0]==-999:
        hmmax=-999
        hmmin=9999
    else:
        hmmax=hmmin=hmeas[0]
    for i in range(0,len(DateInput)):
        if hmeas[i]!=-999:
            if hmeas[i]>hmmax:
                hmmax=h[i]
            elif hmeas[i]<hmmin:
                hmmin=h[i]
    if hmin>hmmin:
        hmin=hmmin
    if hmax<hmmax:
        hmax=hmmax
    ybuffer=0.1*(hmax-hmin)
    ylim((hmin - ybuffer, hmax + ybuffer))
    ylabel('m')
    legend((r'h obs',r'h sim'), loc=0)
    leg = gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    setp(ltext, fontsize='small')    # the legend text fontsize
    xlabel(r'Date')
    grid(True)

    #figure.title('EARTH',fontsize=10)
    subplots_adjust(left=0.05, bottom=0.07, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
#    draw()
#    ion()
    show()
    del DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, Sm, Sr

##    #__________________Export graphs as pdf______________________#
##        ans=''
##        while ans!='y' or ans!='n':
##            ans = str(raw_input('\nDo U want to export the graph as a pdf file?\n(y or n)'))
##            if ans == 'n':
##                break
##            elif ans == 'y':
##                matplotlib.use('PDF')
##                savefig('C:\_alf\MOD9_GEOPROC\P2\EARTHgraph.pdf',dpi=600)
##                print 'Graph saved as ' + 'C:\_alf\MOD9_GEOPROC\P2\EARTHgraph.pdf'
##                break
