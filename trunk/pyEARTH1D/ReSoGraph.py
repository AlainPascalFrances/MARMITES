# -*- coding: UTF-8 -*-
import matplotlib as mpl
import matplotlib.pyplot as plt

##=========================================================================================##
##===========================| PIEZO GRAPHS |========================================##
##=========================================================================================##

def piezocalibGRAPH(DateInput,h,hmeas,P,Pe):
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
    first_date = mpl.dates.num2date(DateInput[0])
    months=mpl.dates.MonthLocator(bymonth=first_date.month,
    bymonthday=first_date.day, interval=1)
    monthsFmt=mpl.dates.DateFormatter('%y-%m')

#__________________Create outputs plots______________________#
#    ioff()
    figCalib=plt.figure(num=None, figsize=(11.7, 8.27), dpi=150)
    figCalib.Title='Calibration graphs'

    ax2=plt.subplot(111)
    plt.setp(ax2.get_xticklabels(), fontsize=8)
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,hmeas,'o', markersize=5, markerfacecolor = 'lime', markeredgecolor='green')
    plt.plot_date(DateInput,h,'-', color = 'blue')
    ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%1.1f'))
    ax2.xaxis.set_major_locator(months)
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
                hmmax=hmeas[i]
            elif hmeas[i]<hmmin:
                hmmin=hmeas[i]
    hmin=min(hmin,hmmin)
    hmax=max(hmax,hmmax)
    ybuffer=0.1*(hmax-hmin)
    plt.ylim((hmin - ybuffer, hmax + ybuffer))
    plt.ylabel('m')
    plt.legend((r'h obs',r'h sim'), loc=0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.xlabel(r'Date')
#    ax2.autoscale_view()
    plt.grid(True)
    ax2.xaxis.set_major_formatter(monthsFmt)
    plt.xlim((DateInput[0],DateInput[len(h)-1]))
    ax2.tick_params(axis='x', labelrotation=90)
    for label in ax2.get_xticklabels():
        label.set_verticalalignment('top')
    ax2.tick_params(axis='y', labelrotation=90)
    for label in ax2.get_yticklabels():
        label.set_verticalalignment('center')

    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.95, wspace=0.1, hspace=0.05)
#    draw()
#    ion()
    plt.show()
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

    first_date = mpl.dates.num2date(DateInput[0])
    months=mpl.dates.MonthLocator(bymonth=first_date.month,
    bymonthday=first_date.day, interval=1)
    monthsFmt=mpl.dates.DateFormatter('%y-%m')

#__________________Create outputs plots______________________#
#    ioff()
    figCalib=plt.figure(num=None, figsize=(11.7, 8.27), dpi=150)
    figCalib.Title='Calibration graphs'

    ax5=plt.subplot(515)
    plt.setp( ax5.get_xticklabels(), fontsize=8)
    plt.setp( ax5.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,hmeas,'o', markersize=5, markerfacecolor = 'lime', markeredgecolor='green')
    plt.plot_date(DateInput,h,'-', color = 'blue')
    ax5.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%1.1f'))
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
                hmmax=hmeas[i]
            elif hmeas[i]<hmmin:
                hmmin=hmeas[i]
    hmin=min(hmin,hmmin)
    hmax=max(hmax,hmmax)
    ybuffer=0.1*(hmax-hmin)
    plt.ylim((hmin - ybuffer, hmax + ybuffer))
    plt.ylabel('m')
    plt.legend((r'h obs',r'h sim'), loc=0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.xlabel(r'Date')
#    ax2.autoscale_view()
    ax5.tick_params(axis='x', labelrotation=90)
    for label in ax5.get_xticklabels():
        label.set_verticalalignment('top')
    plt.grid(True)
    ax5.tick_params(axis='y', labelrotation=90)
    for label in ax5.get_yticklabels():
        label.set_verticalalignment('center')
    
    ax4=plt.subplot(514, sharex=ax5)
    plt.setp( ax4.get_xticklabels(), visible=False)
    plt.setp( ax4.get_yticklabels(), fontsize=8)
    plt.grid(True)
    plt.plot_date(DateInput,R,'-')
    plt.legend('R', loc=0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.ylabel('mm')
    ax4.tick_params(axis='y', labelrotation=90)
    for label in ax4.get_yticklabels():
        label.set_verticalalignment('center')
    
    ax3=plt.subplot(513, sharex=ax5)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput, Smeas, 'o', markersize=5, markerfacecolor = 'lime', markeredgecolor='green')
    plt.plot_date(DateInput,S, '-', color = 'brown')
    plt.legend((r'S obs',r'S sim'), loc=0)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
#    ybuffer=0.01*(float(Sm)-float(Sr))
#    ylim((float(Sr) - ybuffer,float(Sm) + ybuffer))
#    ylabel('mm')
#    ylim(0,1)
    plt.ylim(float(Sr),float(Sm))
    ax3.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%0.2f'))
    ax3.tick_params(axis='y', labelrotation=90)
    for label in ax3.get_yticklabels():
        label.set_verticalalignment('center')
    plt.grid(True)
    plt.ylabel('%')

    ax2=plt.subplot(512, sharex=ax5)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,PET,'b-')
    plt.plot_date(DateInput,ETa,'r-')
    plt.legend((r'PET',r'ETa'), loc=0)  #,  fontsize=10)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.grid(True)
    plt.ylabel('mm')
    ax2.tick_params(axis='y', labelrotation=90)
    for label in ax2.get_yticklabels():
        label.set_verticalalignment('center')

    ax1=plt.subplot(511, sharex=ax5)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), fontsize=8)
#    DateInput1=range(0,len(DateInput))
    plt.bar(DateInput,P,color='b', linewidth=0, align = 'edge', label = 'P')
    plt.bar(DateInput,Pe,color='deepskyblue', linewidth=0, align = 'edge', label = 'Pe')
    plt.legend(loc=0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.grid(True)
    plt.ylabel('mm')
    ax1.xaxis.set_major_locator(months)
    ax1.xaxis.set_major_formatter(monthsFmt)
    plt.xlim((DateInput[0],DateInput[len(h)-1]))
    ax1.tick_params(axis='y', labelrotation=90)
    for label in ax1.get_yticklabels():
        label.set_verticalalignment('center')

    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.95, wspace=0.1, hspace=0.05)
#    draw()
#    ion()
    plt.show()
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

    first_date = mpl.dates.num2date(DateInput[0])
    months=mpl.dates.MonthLocator(bymonth=first_date.month,
    bymonthday=first_date.day, interval=1)
    monthsFmt=mpl.dates.DateFormatter('%y-%m')

#__________________Create outputs plots______________________#
#    ioff()
    figCalib=plt.figure(num=None, figsize=(11.7, 8.27), dpi=150)
    figCalib.Title='All graphs'

# First column of graphs

    ax7=plt.subplot(427)
    plt.setp(ax7.get_xticklabels(), fontsize=8)
    plt.setp(ax7.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,Qs,'r-')
    ax7.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%1.1f'))
    ax7.xaxis.set_major_locator(months)
    ax7.xaxis.set_major_formatter(monthsFmt)
    plt.bar(DateInput, SUST, linewidth=0, align = 'edge')
    plt.xlim((DateInput[0],DateInput[len(S)-1]))
    plt.ylabel('mm')
    plt.legend(['Qs', 'SUST'], loc=0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.xlabel(r'Date')
#    xlim((0,len(P)))
    plt.grid(True)
    ax7.tick_params(axis='x', labelrotation=90)
    for label in ax7.get_yticklabels():
        label.set_verticalalignment('center')
    ax7.tick_params(axis='y', labelrotation=90)
    for label in ax7.get_yticklabels():
        label.set_verticalalignment('center')
    
    ax5=plt.subplot(425, sharex=ax7)
    plt.setp(ax5.get_xticklabels(), visible=False)
    plt.setp(ax5.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput, Smeas, 'o', markersize=5, markerfacecolor = 'lime', markeredgecolor='green')
    plt.plot_date(DateInput,S, '-', color = 'brown')
    ybuffer=0.1*(float(Sm)-float(Sr))
    plt.ylim((float(Sr) - ybuffer,float(Sm) + ybuffer))
    plt.ylabel('%')
    ax5.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%1.1f'))
    plt.legend((r'S obs',r'S sim'), loc=0)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.grid(True)
    ax5.tick_params(axis='y', labelrotation=90)
    for label in ax5.get_yticklabels():
        label.set_verticalalignment('center')

    ax3=plt.subplot(423, sharex=ax7)
    plt.setp( ax3.get_xticklabels(), visible=False)
    plt.setp( ax3.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,PET,'b-')
    plt.plot_date(DateInput,ETa,'r-')
    plt.legend((r'PET',r'ETa'), loc=0)  #,  fontsize=10)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.grid(True)
    plt.ylabel('mm')
    ax3.tick_params(axis='y', labelrotation=90)
    for label in ax3.get_yticklabels():
        label.set_verticalalignment('center')

    ax1=plt.subplot(421, sharex=ax7)
    plt.setp( ax1.get_xticklabels(), visible=False)
    plt.setp( ax1.get_yticklabels(), fontsize=8)
#    DateInput1=range(0,len(DateInput))
    plt.bar(DateInput,P,color='b', linewidth=0, align = 'edge')
    plt.legend(['P'], loc=0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.grid(True)
    plt.ylabel('mm')
    ax1.tick_params(axis='y', labelrotation=90)
    for label in ax1.get_yticklabels():
        label.set_verticalalignment('center')

# Second column of graphs
    ax2=plt.subplot(422, sharex=ax7)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), fontsize=8)
#    DateInput1=range(0,len(DateInput))
    plt.bar(DateInput,Pe,color='deepskyblue', linewidth=0, align = 'edge')
    plt.legend(['Pe'], loc=0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.grid(True)
    plt.ylabel('mm')
    ax2.tick_params(axis='y', labelrotation=90)
    for label in ax2.get_yticklabels():
        label.set_verticalalignment('center')

    ax4=plt.subplot(424, sharex=ax7)
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), fontsize=8)
#    DateInput1=range(0,len(DateInput))
    plt.bar(DateInput,Rp,linewidth=0, align = 'edge')
    plt.legend(['Rp'], loc=0)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.grid(True)
    plt.ylabel('mm')
    ax4.tick_params(axis='y', labelrotation=90)
    for label in ax4.get_yticklabels():
        label.set_verticalalignment('center')

    ax6=plt.subplot(426, sharex=ax7)
    plt.setp(ax6.get_xticklabels(), visible=False)
    plt.setp(ax6.get_yticklabels(), fontsize=8)
    plt.grid(True)
    plt.plot_date(DateInput,R,'-')
    plt.legend('R', loc=0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.ylabel('mm')
    ax6.tick_params(axis='y', labelrotation=90)
    for label in ax6.get_yticklabels():
        label.set_verticalalignment('center')

    ax8=plt.subplot(428, sharex=ax7)
    plt.setp(ax8.get_xticklabels(), fontsize=8)
    plt.setp(ax8.get_yticklabels(), fontsize=8)
    plt.plot_date(DateInput,hmeas,'o', markersize=5, markerfacecolor = 'lime', markeredgecolor='green')
    plt.plot_date(DateInput,h,'-', color = 'blue')
    ax8.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
    ax8.xaxis.set_major_locator(months)
    ax8.xaxis.set_major_formatter(monthsFmt)
#    ax2.autoscale_view()
    plt.xlim((DateInput[0],DateInput[len(S)-1]))
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
    hmin=min(hmin,hmmin)
    hmax=max(hmax,hmmax)
    ybuffer=0.1*(hmax-hmin)
    plt.ylim((hmin - ybuffer, hmax + ybuffer))
    plt.ylabel('m')
    plt.legend((r'h obs',r'h sim'), loc=0)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.xlabel(r'Date')
    plt.grid(True)
    ax8.tick_params(axis='x', labelrotation=90)
    for label in ax8.get_xticklabels():
        label.set_verticalalignment('top')
    ax8.tick_params(axis='y', labelrotation=90)
    for label in ax8.get_yticklabels():
        label.set_verticalalignment('center')

    #figure.title('EARTH',fontsize=10)
    plt.subplots_adjust(left=0.05, bottom=0.07, right=0.95, top=0.95, wspace=0.1, hspace=0.1)
#    draw()
#    ion()
    plt.show()
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
