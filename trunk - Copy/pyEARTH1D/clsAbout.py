# -*- coding: iso-8859-1 -*-
# generated by wxGlade 0.6.3 on Fri Nov 19 00:38:46 2010

import wx
import os

# begin wxGlade: dependencies
# end wxGlade

# begin wxGlade: extracode

# end wxGlade

class clsAbout(wx.Dialog):
    def __init__(self, *args, **kwds):
        # begin wxGlade: clsAbout.__init__
        kwds["style"] = wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER|wx.THICK_FRAME|wx.STAY_ON_TOP
        wx.Dialog.__init__(self, *args, **kwds)
        self.lbl_About = wx.StaticText(self, -1, u"\n-x-x-x-x-x-x-x-x-\n\npyEARTH1D\n\nBased on:\nVan der Lee J. and Gehrels, J. (1990)\nModelling Aquifer Recharge - Introduction to the Lumped Parameter Model EARTH\nFree University of Amsterdam, The Netherlands\n\nSee also:\nFranc�s, A. P. (2008)\nSpatio - temporal groundwater recharge assessment : a data - integration and modelling approach\nM.Sc. thesis, ITC-WRS, Enschede\n\nContact:\na.p.frances@utwente.nl\n\n-x-x-x-x-x-x-x-x-\n", style=wx.ALIGN_CENTRE)

        self.__set_properties()
        self.__do_layout()
        # end wxGlade

    def __set_properties(self):
        _icon = wx.EmptyIcon()
        _icon.CopyFromBitmap(wx.Bitmap(os.getcwd() + '\\aquaflagged.ico', wx.BITMAP_TYPE_ANY))
        self.SetIcon(_icon)
        # begin wxGlade: clsAbout.__set_properties
        self.SetTitle("pyEARTH1D - Recharge and Soil moisture modelling")
        self.SetFocus()
        # end wxGlade

    def __do_layout(self):
        # begin wxGlade: clsAbout.__do_layout
        sz_About = wx.BoxSizer(wx.VERTICAL)
        sz_About.Add(self.lbl_About, 0, wx.ALIGN_CENTER_HORIZONTAL|wx.ALIGN_CENTER_VERTICAL, 0)
        self.SetSizer(sz_About)
        sz_About.Fit(self)
        self.Layout()
        # end wxGlade

# end of class clsAbout

