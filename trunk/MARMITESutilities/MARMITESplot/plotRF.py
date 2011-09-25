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
import sys
import CreateColors

def plot(strTitle = 'Title', x = [], \
        y1 = [], y2 = [], y3 =[], y4=[],\
        lbl_y1 = '', lbl_y2 = '', lbl_y3 = '', lbl_y4 = '', lbl_veg = ''
        ,plot_exportRF_fn = ''
        ):


##    if len(x)>365:
##        locX = mpl.dates.MonthLocator()
##    else:
##        locX = mpl.dates.HourLocator(byhour=range(24), interval = 24)
    monthsFmt = mpl.dates.DateFormatter('%y-%m-%d %H:%M')

    colors_y3 = CreateColors.main(hi=120, hf=130, numbcolors = (len(y3)))
    colors_y4 = CreateColors.main(hi=160, hf=180, numbcolors = (len(y4)))
    lbls_y3 = [lbl_y3 + " " +  lbl for lbl in lbl_veg]
    lbls_y4 = [lbl_y4 + " " +  lbl for lbl in lbl_veg]

    fig = plt.figure(figsize=(2*11.7, 2*8.27))

    ax1=fig.add_subplot(4,1,4)
    plt.setp( ax1.get_xticklabels(), fontsize=8)
    plt.setp( ax1.get_yticklabels(), fontsize=8)
    plt.bar( x, y1, color='b', linewidth = 0, align='edge', width = 0.8, label=lbl_y1)
#    ax1.yaxis.set_major_formatter(plt.FormatStrFormatter('%1.2f'))
    labels=ax1.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
    labels=ax1.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    #ax1.xaxis.set_major_locator(locX)
    ax1.xaxis.set_major_formatter(monthsFmt)
    ax1.yaxis.set_major_formatter(plt.FormatStrFormatter('%0.1f'))
    lbls_y1 = []
    lbls_y1.append(lbl_y1)
    plt.legend(loc = 0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.xlabel(r'Date')
    plt.grid(True)

    ax2=fig.add_subplot(4,1,3, sharex=ax1)
    plt.setp( ax2.get_xticklabels(), visible=False)
    plt.setp( ax2.get_yticklabels(), fontsize=8)
    for i, (y, color, lbl) in enumerate(zip(y4, colors_y4, lbls_y4)) :
        ax2.bar(x+(0.8*float(i)/len(y4)), y, color=color, label=lbl, linewidth=0, align='edge', width=(0.8/len(y4)))
#    ax2.yaxis.set_major_formatter(plt.FormatStrFormatter('%1.2f'))
    labels=ax2.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
#    ax2.set_xticks(x+(0.4/len(y4)))
    labels=ax2.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    #ax2.xaxis.set_major_locator(locX)
    ax2.xaxis.set_major_formatter(monthsFmt)
    ax2.yaxis.set_major_formatter(plt.FormatStrFormatter('%0.1f'))
    plt.legend(loc = 0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.grid(True)

    ax3=fig.add_subplot(4,1,2, sharex=ax1)
    plt.setp( ax3.get_xticklabels(), visible=False)
    plt.setp( ax3.get_yticklabels(), fontsize=8)
    plt.bar(x,y2,color='aqua', linewidth=0, align='edge', width = 0.8, label=lbl_y2)
#    ax3.yaxis.set_major_formatter(plt.FormatStrFormatter('%1.2f'))
    labels=ax3.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
    labels=ax3.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    #ax3.xaxis.set_major_locator(locX)
    ax3.xaxis.set_major_formatter(monthsFmt)
    ax3.yaxis.set_major_formatter(plt.FormatStrFormatter('%0.1f'))
    lbls_y2 = []
    lbls_y2.append(lbl_y2)
    plt.legend(loc = 0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.grid(True)

    ax4=fig.add_subplot(4,1,1, sharex=ax1)
    ax4.set_title(strTitle)
    plt.setp( ax4.get_xticklabels(), visible=False)
    plt.setp( ax4.get_yticklabels(), fontsize=8)

    for i, (y, color, lbl) in enumerate(zip(y3, colors_y3, lbls_y3)) :
        ax4.bar(x+(0.8*float(i)/len(y3)), y, color=color, label=lbl, linewidth=0, align='edge', width=(0.8/len(y3)))
#    ax1.yaxis.set_major_formatter(plt.FormatStrFormatter('%1.2f'))
    labels=ax4.get_yticklabels()
    plt.setp(labels, 'rotation', 90)
#    ax4.set_xticks(x+(0.4/len(y3)))
    labels=ax4.get_xticklabels()
    plt.setp(labels, 'rotation', 90)
    #ax4.xaxis.set_major_locator(locX)
    ax4.xaxis.set_major_formatter(monthsFmt)
    ax4.yaxis.set_major_formatter(plt.FormatStrFormatter('%0.1f'))
    plt.legend(loc = 0)
    leg = plt.gca().get_legend()
    ltext  = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize='small')    # the legend text fontsize
    plt.grid(True)

    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.95, wspace=0.1, hspace=0.15)
    #plt.show()
    plt.savefig(plot_exportRF_fn,dpi=150)
    plt.close()
# EOF
#############################