#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      alf
#
# Created:     25-11-2010
# Copyright:   (c) alf 2010
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
import matplotlib as mpl
if mpl.get_backend()<>'agg':
    mpl.use('agg')
import matplotlib.pyplot as plt
import CreateColors

def plot(strTitle = 'Title', x = [], \
        y1 = [], y2 = [] , y3 = [], lbl_y1 = '', lbl_y2 = 'E0', lbl_y3 = ''
        ,plot_exportPET_fn= ''):
##    if len(x)>35:
##        locX = mpl.dates.MonthLocator()
##    else:
##        locX = mpl.dates.HourLocator(byhour=range(24), interval = 24)
    monthsFmt = mpl.dates.DateFormatter('%y-%m-%d %H:%M')

    fig=plt.figure(figsize=(2*11.7, 2*8.27))
    colors_y1 = CreateColors.main(hi=100, hf=110, numbcolors = (len(y1)-1))
    colors_y3 = CreateColors.main(hi=70, hf=80, numbcolors = (len(y3)))
    lbl_y = []
    lbl_y.append(lbl_y2)
    for i in range(len(lbl_y1)):
        lbl_y.append(lbl_y1[i])
    for i in range(len(lbl_y3)):
        lbl_y.append(lbl_y3[i])

    ax2=plt.subplot(111)
    ax2.set_title(strTitle)
    plt.setp( ax2.get_xticklabels(), fontsize=8)
    plt.setp( ax2.get_yticklabels(), fontsize=8)
    plt.plot_date(x,y2,'-', color='blue')
    plt.plot_date(x,y1[0],'-', color='green')
    for i in range(1,len(y1)):
        plt.plot_date(x,y1[i],'-.', color=colors_y1[i-1], linewidth = 2)
    for i in range(len(y3)):
        plt.plot_date(x,y3[i],'-.', color=colors_y3[i], linewidth = 2)
    ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%1.2f'))
    labels=ax2.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
    labels=ax2.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    #ax2.xaxis.set_major_locator(locX)
    ax2.xaxis.set_major_formatter(monthsFmt)
    plt.ylabel('mm/day')
    plt.legend(lbl_y, loc=0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.xlabel(r'Date')
    plt.grid(True)

    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.95, wspace=0.1, hspace=0.05)
    #show()
    plt.savefig(plot_exportPET_fn,dpi=150)
    plt.close()

def plotVAR(strTitle = 'Title', x = [], \
        y1 = [], y2 = [], y3 =[], y4=[], y5=[], y6=[]\
        ,lbl_y1 = '', lbl_y2 = '', lbl_y3 = '', lbl_y4 = '', lbl_y5 = '', lbl_y6 = ''
        , plot_exportVAR_fn = ''
        ):

##    if len(x)>365:
##        locX = mpl.dates.MonthLocator()
##    else:
##        locX = mpl.dates.HourLocator(byhour=range(24), interval = 24)
    DateFmt = mpl.dates.DateFormatter('%y-%m-%d %H:%M')

    fig = plt.figure(figsize=(2*11.7, 2*8.27))

    ax2=plt.subplot(111)
    ax2.set_title(strTitle)
    plt.setp( ax2.get_xticklabels(), fontsize=8)
    plt.setp( ax2.get_yticklabels(), fontsize=8)
    plt.plot_date(x,y1,'-', color='red')
    if y2<>[]:
        plt.plot_date(x,y2,'-', color='yellow')
    if y3<>[]:
        plt.plot_date(x,y3,'-.', color='orange', linewidth = 2)
    if y4<>[]:
        plt.plot_date(x,y4,'-.', color='blue')
    if y5<>[]:
        plt.plot_date(x,y5,'--', color='green')
    if y6<>[]:
        plt.plot_date(x,y6,'--', color='black')
    ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%1.2f'))
    labels=ax2.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
    labels=ax2.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    #ax2.xaxis.set_major_locator(locX)
    ax2.xaxis.set_major_formatter(DateFmt)
    plt.ylabel('mm/day')
    plt.legend((lbl_y1,lbl_y2,lbl_y3,lbl_y4,lbl_y5, lbl_y6), loc=0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.xlabel(r'Date')
    plt.grid(True)

    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.95, wspace=0.1, hspace=0.05)
#    plt.show()
    plt.savefig(plot_exportVAR_fn,dpi=150)
    plt.close()
# EOF
##############################
