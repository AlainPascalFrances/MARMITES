# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:        startMARMITES
# Purpose:
#
# Author:      frances08512
#
# Created:     25-11-2010
# Copyright:   (c) frances08512 2010
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

""" See info in MARMITESsoil_v3.py"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.3"
__date__ = "2012"

import sys, os, h5py, shutil, itertools
from win32com.shell import shell, shellcon
import matplotlib as mpl
if mpl.get_backend!='agg':
    mpl.use('agg')
import matplotlib.pyplot as plt
plt.ioff()
import numpy as np
import startMARMITESsurface as startMMsurf
import MARMITESsoil_v3 as MMsoil
import ppMODFLOW_flopy_v3 as ppMF
import MARMITESplot_v3 as MMplot
import MARMITESutilities as MMutils
#import MARMITESprocess_v3 as MMproc

# TODO verify if thickness of MF layer 1 > thickness of soil

#####################################

timestart = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
fmt_DH = mpl.dates.DateFormatter('%Y-%m-%d %H:%M')
print '\n##############\nMARMITES v0.3 started!\n%s\n##############' % mpl.dates.DateFormatter.format_data(fmt_DH, timestart)

cUTIL = MMutils.clsUTILITIES(fmt = fmt_DH)

# workspace (ws) definition
# read input file (called _input.ini in the MARMITES workspace
# the first character on the first line has to be the character used to comment
# the file can contain any comments as the user wish, but the sequence of the input has to be respected
# 00_TESTS\MARMITESv3_r13c6l2  00_TESTS\r40c20  00_TESTS\r20c40  r130c60l2   r130c60l2new
# SARDON2013  CARRIZAL3 CARRIZAL3newera LAMATA LaMata_new
startMM_fn = shell.SHGetFolderPath (0, shellcon.CSIDL_DESKTOP, 0, 0)
inputFile = cUTIL.readFile(startMM_fn, 'startMM_fn.txt')

try:
    MM_ws = str(inputFile[0].strip())
    MM_ini_fn = str(inputFile[1].strip())
except:
    cUTIL.ErrorExit(msg = '\nFATAL ERROR!\nInput ini MM file [%s] is not correct.' % inputFile)

inputFile = cUTIL.readFile(MM_ws,MM_ini_fn)

l=0
try:
    # run_name
    run_name = inputFile[l].strip()
    l += 1
    fmt_DHshort = mpl.dates.DateFormatter('%Y%m%d%H%M')
    MM_ws_out = os.path.join(MM_ws,'out_%s_%s'% (mpl.dates.DateFormatter.format_data(fmt_DHshort, timestart), run_name))
    f = 1
    if os.path.exists(MM_ws_out):
        MM_ws_out1 = MM_ws_out
        while os.path.exists(MM_ws_out1):
            MM_ws_out1 = os.path.join(MM_ws,'out_%s_%s_%s' % (mpl.dates.DateFormatter.format_data(fmt_DHshort, timestart), run_name, f))
            f +=1
        os.makedirs(MM_ws_out1)
        MM_ws_out = MM_ws_out1
        del MM_ws_out1
    else:
        os.makedirs(MM_ws_out)
    # # ECHO ON/OFF (1 - interpreter verbose BUT NO report, 0 - NO interpreter verbose BUT report)
    verbose = int(inputFile[l].strip())
    l += 1
    # output plot (1 is YES, 0 is NO)
    plt_out  = int(inputFile[l].strip())
    l += 1
    plt_freq =  int(inputFile[l].strip())
    l += 1
    nrangeMM =  int(inputFile[l].strip())
    l += 1
    nrangeMF =  int(inputFile[l].strip())
    l += 1
    ctrsMM =  int(inputFile[l].strip())
    if ctrsMM == 1:
        ctrsMM = True
    else:
        ctrsMM = False
    l += 1
    ctrsMF =  int(inputFile[l].strip())
    if ctrsMF == 1:
        ctrsMF = True
    else:
        ctrsMF = False
    l += 1
    ntick =  int(inputFile[l].strip())
    l += 1
    animation =  int(inputFile[l].strip())
    # read observations?
    l += 1
    plt_out_obs = int(inputFile[l].strip())
    l += 1
    plt_WB_unit = inputFile[l].strip()
    if plt_WB_unit == 'year':
        facTim = 365.0
    else:
        facTim = 1.0
    l += 1
    iniMonthHydroYear = int(inputFile[l].strip())
    if iniMonthHydroYear < 1 or iniMonthHydroYear > 12:
        print '\nWARNING!\nInvalid starting month of the hydrologic year. Please correct your input (currently %d) to a value between 1 and 12 inclusive. The starting month was now defined as October (10).' % iniMonthHydroYear
        iniMonthHydroYear = 10
    l += 1
    # map of inputs (1 is YES, 0 is NO)
    plt_input = int(inputFile[l].strip())
    l += 1
    #run MARMITESsurface  (1 is YES, 0 is NO)
    MMsurf_yn = int(inputFile[l].strip())
    l += 1
    #plot MARMITESsurface  (1 is YES, 0 is NO)
    MMsurf_plot = int(inputFile[l].strip())
    l += 1
    #run MARMITESsoil  (1 is YES, 0 is NO)
    MMsoil_yn = int(inputFile[l].strip())
    l += 1
    #run MODFLOW  (1 is YES, 0 is NO)
    MF_yn = int(inputFile[l].strip())
    l += 1
    convcrit = float(inputFile[l].strip())
    l += 1
    convcritmax = float(inputFile[l].strip())
    l += 1
    ccnum = int(inputFile[l].strip())
    l += 1
    # Define MARMITESsurface folder
    MMsurf_ws = inputFile[l].strip()
    l += 1
    # METEO/VEGETATION/SOIL/WATER PARAMETERS file name
    inputFile_PAR_fn = inputFile[l].strip()
    l += 1
    # METEO TIME SERIES file name
    inputFile_TS_fn = inputFile[l].strip()
    l += 1
    # OPTIONNAL IRRIGATION FILES
    irr_yn = int(inputFile[l].strip())
    if irr_yn == 1 :
        l += 1
        inputFile_TSirr_fn = inputFile[l].strip()
        l += 1
        gridIRR_fn = inputFile[l].strip()
        l += 1
    else:
        l += 3
    # ouputprefix
    outputFILE_fn = inputFile[l].strip()
    l += 1
    # ZONEVEGSOILfile
    outMMsurf_fn = inputFile[l].strip()
    l += 1
    # Define MODFLOW ws folders
    MF_ws = inputFile[l].strip()
    l += 1
    MF_ini_fn = inputFile[l].strip()
    l += 1
    #GRID (ll means lower left)
    xllcorner = float(inputFile[l].strip())
    l += 1
    yllcorner = float(inputFile[l].strip())
    l += 1
    gridMETEO_fn = inputFile[l].strip()
    l += 1
    gridSOIL_fn = inputFile[l].strip()
    l += 1
    gridSOILthick_fn = inputFile[l].strip()
    l += 1
    gridSsurfhmax_fn =  inputFile[l].strip()
    l += 1
    gridSsurfw_fn =  inputFile[l].strip()
    l += 1
    SOILparam_fn = inputFile[l].strip()
    l += 1
    inputObs_fn = inputFile[l].strip()
    l += 1
    inputObsHEADS_fn = inputFile[l].strip()
    l += 1
    inputObsSM_fn = inputFile[l].strip()
    l += 1
    inputObsRo_fn = inputFile[l].strip()
    l += 1
    rmseHEADSmax  = float(inputFile[l].strip())
    l += 1
    rmseSMmax  = float(inputFile[l].strip())
    l += 1
    chunks = int(inputFile[l].strip())
    if MMsoil_yn == 1:
        MF_yn = 1
    if MMsurf_plot == 1:
        plt_out = 0
        plt_out_obs = 0
        print "\nYou required the MMsurf plots to appear on the screen. Due to backends limitation, MM and MF plots were disabled. Run again MM with MMsurf_plot = 0 to obtain the MM and MF plots."
    if animation == 1:
        test_ffmpeg = cUTIL.which(program = 'ffmpeg.exe')
        if test_ffmpeg == None:
            print '\nWARNING!\nFFmpeg was not found on your computer!\nSpatial animations will not be produced.\nYou can solve this problem downloading FFmpeg at www.ffmpeg.org and install it.'
            animation = 0
        del test_ffmpeg
except:
    cUTIL.ErrorExit(msg = '\nFATAL ERROR!\nType error in the input file %s' % MM_ini_fn)

del inputFile

shutil.copy2(os.path.join(MM_ws,MM_ini_fn), os.path.join(MM_ws_out,'__%s%s'% (mpl.dates.DateFormatter.format_data(fmt_DHshort, timestart), MM_ini_fn)))
shutil.copy2(os.path.join(MM_ws, MF_ws,MF_ini_fn), os.path.join(MM_ws_out,'__%s%s'% (mpl.dates.DateFormatter.format_data(fmt_DHshort, timestart), MF_ini_fn)))
shutil.copy2(os.path.join(MM_ws, MMsurf_ws,inputFile_PAR_fn), os.path.join(MM_ws_out,'__%s%s'% (mpl.dates.DateFormatter.format_data(fmt_DHshort, timestart), inputFile_PAR_fn)))
shutil.copy2(os.path.join(MM_ws, SOILparam_fn), os.path.join(MM_ws_out,'__%s__%s'% (mpl.dates.DateFormatter.format_data(fmt_DHshort, timestart), SOILparam_fn.split('\\')[-1])))

if verbose == 0:
# capture interpreter output to be written in to a report file
    report_fn = os.path.join(MM_ws_out,'__%s_MMMFrun_report.txt' % (mpl.dates.DateFormatter.format_data(fmt_DHshort, timestart)))
    print '\nECHO OFF (no screen output).\nSee the report of the MM-MF run in file:\n%s\n' % report_fn
    s = sys.stdout
    report = open(report_fn, 'w')
    sys.stdout = report
    print '\n##############\nMARMITES v0.3 started!\n%s\n##############' % mpl.dates.DateFormatter.format_data(fmt_DH, timestart)

    cUTIL = MMutils.clsUTILITIES(fmt = fmt_DH, verbose = verbose, s = s, report = report, report_fn = report_fn)
else:
    report = None

MMsurf_ws = os.path.join(MM_ws,MMsurf_ws)
MF_ws = os.path.join(MM_ws,MF_ws)
if os.path.exists(MM_ws):
    if os.path.exists(MMsurf_ws):
        if os.path.exists(MF_ws):
            print ('\nMARMITES workspace:\n%s\n\nMARMITESsurf workspace:\n%s\n\nMODFLOW workspace:\n%s' % (MM_ws, MMsurf_ws, MF_ws))
        else:
            print "The folder %s doesn't exist!\nPlease create it and run the model again." % (MF_ws)
    else:
        print "The folder %s doesn't exist!\nPlease create it and run the model again." % (MMsurf_ws)
else:
    print "The folder %s doesn't exist!\nPlease create it and run the model again." % (MM_ws)

# #############################
# 1st phase: INITIALIZATION #####
# #############################

# #############################
# ###  MARMITES surface  ######
# #############################
print'\n##############'
print 'MARMITESsurf RUN'
if MMsurf_yn>0:
    durationMMsurf = 0.0
    timestartMMsurf = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    if irr_yn == 1 :
        outMMsurf_fn = startMMsurf.MMsurf(cUTIL, MMsurf_ws, inputFile_TS_fn, inputFile_PAR_fn, outputFILE_fn, MM_ws, outMMsurf_fn, MMsurf_plot, inputFile_TSirr_fn)
    else:
        outMMsurf_fn = startMMsurf.MMsurf(cUTIL, MMsurf_ws, inputFile_TS_fn, inputFile_PAR_fn, outputFILE_fn, MM_ws, outMMsurf_fn, MMsurf_plot)
    timeendMMsurf = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    durationMMsurf=(timeendMMsurf-timestartMMsurf)
inputFile = cUTIL.readFile(MM_ws,outMMsurf_fn)
l=0
VegName = []
Zr = []
kT_min = []
kT_max = []
kT_n = []
try:
    NMETEO = int(inputFile[l].strip())
    l += 1
    NVEG = int(inputFile[l].strip())
    l += 1
    NSOIL = int(inputFile[l].strip())
    l += 1
    inputDate_fn = str(inputFile[l].strip())
    l += 1
    inputZON_dSP_RF_veg_fn = str(inputFile[l].strip())
    l += 1
    inputZON_dSP_RFe_veg_fn = str(inputFile[l].strip())
    l += 1
    inputZON_dSP_PT_fn = str(inputFile[l].strip())
    l += 1
    input_dSP_LAI_veg_fn = str(inputFile[l].strip())
    l += 1
    inputZON_dSP_PE_fn = str(inputFile[l].strip())
    l += 1
    inputZON_dSP_E0_fn = str(inputFile[l].strip())
    l += 1
    line = inputFile[l].split()
    for v in range(NVEG):
        VegName.append((line[v]))
    l += 1
    line = inputFile[l].split()
    for v in range(NVEG):
        Zr.append(float(line[v]))
    l += 1
    line = inputFile[l].split()
    for v in range(NVEG):
        kT_min.append(float(line[v]))
    l += 1
    line = inputFile[l].split()
    for v in range(NVEG):
        kT_max.append(float(line[v]))
    l += 1
    line = inputFile[l].split()
    for v in range(NVEG):
        kT_n.append(float(line[v]))
    if irr_yn == 1:
        l += 1
        NCROP = int(inputFile[l].strip())
        l += 1
        NFIELD = int(inputFile[l].strip())
        l += 1
        inputZON_dSP_RF_irr_fn = str(inputFile[l].strip())
        l += 1
        inputZON_dSP_RFe_irr_fn = str(inputFile[l].strip())
        l += 1
        inputZON_dSP_PT_irr_fn = str(inputFile[l].strip())
        l += 1
        Zr_c = []
        kT_min_c = []
        kT_max_c = []
        kT_n_c = []
        line = inputFile[l].split()
        for c in range(NCROP):
            Zr_c.append(float(line[c]))
        Zr_c = np.asarray(Zr_c)
        l += 1
        line = inputFile[l].split()
        for c in range(NCROP):
            kT_min_c.append(float(line[c]))
        kT_min_c = np.asarray(kT_min_c)
        l += 1
        line = inputFile[l].split()
        for c in range(NCROP):
            kT_max_c.append(float(line[c]))
        kT_max_c = np.asarray(kT_max_c)
        l += 1
        line = inputFile[l].split()
        for c in range(NCROP):
            kT_n_c.append(float(line[c]))
        kT_n_c = np.asarray(kT_n_c)
        l += 1
        input_dSP_crop_irr_fn = str(inputFile[l].strip())
except:
    cUTIL.ErrorExit('\nFATAL ERROR!\Error in reading file [%s].' % outMMsurf_fn)

del inputFile
numDays = len(cUTIL.readFile(MM_ws, inputDate_fn))

# Compute hydrologic year start and end dates

inputFile_fn = os.path.join(MM_ws, inputDate_fn)    
if os.path.exists(inputFile_fn):
    date = np.loadtxt(inputFile_fn, skiprows = 1, dtype = str, delimiter = ',')[:,0]
else:
    cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe input file [" + inputFile_fn + "] doesn't exist, verify name and path!")
DATE = np.zeros(len(date), dtype = float)
for i, d in enumerate(date):
    DATE[i] = mpl.dates.datestr2num(d)
del date

year_lst = []
HYindex = []
if mpl.dates.num2date(DATE[0]).month<iniMonthHydroYear or (mpl.dates.num2date(DATE[0]).month == iniMonthHydroYear and mpl.dates.num2date(DATE[0]).day == 1):
    year_lst.append(mpl.dates.num2date(DATE[0]).year)
else:
    year_lst.append(mpl.dates.num2date(DATE[0]).year + 1)
HYindex.append(np.argmax(DATE[:] == mpl.dates.datestr2num('%d-%d-01' % (year_lst[0], iniMonthHydroYear))))

if sum(DATE[:] == mpl.dates.datestr2num('%d-%d-01' % (year_lst[0]+1, iniMonthHydroYear)))==0:
    cUTIL.ErrorExit(msg = '\nThe data file does not contain a full hydrological year starting at date 01/%d' % iniMonthHydroYear)

y = 0    
while DATE[-1] >= mpl.dates.datestr2num('%d-%d-01' % (year_lst[y]+1, iniMonthHydroYear)):
    if iniMonthHydroYear == 1:
        iniMonth_prev = 12
        add = 0
    else:
        iniMonth_prev = iniMonthHydroYear - 1
        add = 1
    indexend = np.argmax(DATE[:] == mpl.dates.datestr2num('%d-%d-30' % (year_lst[y]+add, iniMonth_prev)))
    if DATE[-1] >= mpl.dates.datestr2num('%d-%d-30' % (year_lst[y]+2, iniMonth_prev)):
        year_lst.append(year_lst[y]+1)
        HYindex.append(np.argmax(DATE[:] == mpl.dates.datestr2num('%d-%d-01' % (year_lst[y]+1, iniMonthHydroYear))))
        indexend = np.argmax(DATE[:] == mpl.dates.datestr2num('%d-%d-30' % (year_lst[y]+2, iniMonth_prev)))
        y += 1
    else:
        break
HYindex.append(indexend)
HYindex.insert(0, 0)
del indexend
StartDate = DATE[HYindex[1]]
EndDate   = DATE[HYindex[-1]]
        
#    print year_lst
#    print index
print '-------\nStarting date of time serie:\n%s' % (mpl.dates.DateFormatter.format_data(fmt_DH, DATE[0]))
print '-------\nStarting date of hydrological year(s):'
for j in HYindex[1:-1]:
    print mpl.dates.DateFormatter.format_data(fmt_DH, DATE[j])
print 'End date of last hydrological year:'
print mpl.dates.DateFormatter.format_data(fmt_DH, DATE[HYindex[-1]])
print '-------\nEnd date of time serie:\n%s' % (mpl.dates.DateFormatter.format_data(fmt_DH, DATE[-1]))

# #############################
# ###  READ MODFLOW CONFIG ####
# #############################
print'\n##############'
print 'MODFLOW initialization'
durationMF = 0.0
timestartMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
cMF = ppMF.clsMF(cUTIL, MM_ws, MM_ws_out, MF_ws, MF_ini_fn, xllcorner, yllcorner, numDays = numDays)
# compute cbc conversion factor from volume to mm
if cMF.lenuni == 1:
    conv_fact = 304.8
    lenuni_str = 'feet'
elif cMF.lenuni == 2:
    conv_fact = 1000.0
    lenuni_str = 'm'
elif cMF.lenuni == 3:
    conv_fact = 10.0
    lenuni_str = 'cm'
else:
    cUTIL.ErrorExit('\nFATAL ERROR!\nDefine the length unit in the MODFLOW ini file!\n (see USGS Open-File Report 00-92)')
    # TODO if lenuni!=2 apply conversion factor to delr, delc, etc...
if cMF.laytyp[0]==0:
    cUTIL.ErrorExit('\nFATAL ERROR!\nThe first layer cannot be confined type!\nChange your parameter laytyp in the MODFLOW lpf package.\n(see USGS Open-File Report 00-92)')
if cMF.itmuni != 4:
    cUTIL.ErrorExit('\nFATAL ERROR!\nTime unit is not in days!')
else:
    itmuni_str = 'day'
ncell_MF = []
ncell_MM = []
iboundBOL = np.ones(np.array(cMF.ibound).shape, dtype = bool)
mask_tmp = np.zeros((cMF.nrow, cMF.ncol), dtype = int)
mask = []
mask_Lsup = np.zeros((cMF.nrow, cMF.ncol), dtype = int)
for l in range(cMF.nlay):
    mask_Lsup += np.asarray(cMF.ibound)[l,:,:]
    ncell_MF.append((np.asarray(cMF.ibound)[l,:,:] != 0).sum())
    ncell_MM.append((np.asarray(cMF.iuzfbnd) == l+1).sum())
    iboundBOL[l,:,:] = (np.asarray(cMF.ibound)[l,:,:] != 0)
    mask.append(np.ma.make_mask(iboundBOL[l,:,:]-1))
    mask_tmp += (np.asarray(cMF.ibound)[l,:,:] != 0)
maskAllL = (mask_tmp == 0)
del iboundBOL, mask_tmp
timeendMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
durationMF +=  timeendMF-timestartMF

# #############################
# ### MF time processing
# #############################
# if required by user, compute nper, perlen,etc based on RF analysis in the METEO zones
if isinstance(cMF.nper, str):
    try:
        perlenmax = int(cMF.nper.split()[1].strip())
    except:
        cUTIL.ErrorExit('\nFATAL ERROR!\nError in nper format of the MODFLOW ini file!')
if irr_yn == 0:
    cMF.ppMFtime(inputDate_fn, inputZON_dSP_RF_veg_fn, inputZON_dSP_RFe_veg_fn, inputZON_dSP_PT_fn,input_dSP_LAI_veg_fn, inputZON_dSP_PE_fn, inputZON_dSP_E0_fn, NMETEO, NVEG, NSOIL)
else:
    cMF.ppMFtime(inputDate_fn, inputZON_dSP_RF_veg_fn, inputZON_dSP_RFe_veg_fn, inputZON_dSP_PT_fn, input_dSP_LAI_veg_fn, inputZON_dSP_PE_fn, inputZON_dSP_E0_fn, NMETEO, NVEG, NSOIL, inputZON_dSP_RF_irr_fn, inputZON_dSP_RFe_irr_fn, inputZON_dSP_PT_irr_fn, input_dSP_crop_irr_fn, NFIELD)

print'\n##############'
print 'MARMITESsoil initialization'
MM_SOIL = MMsoil.clsMMsoil(hnoflo = cMF.hnoflo)
MM_SATFLOW = MMsoil.SATFLOW()

# READ input ESRI ASCII rasters
print "\nImporting ESRI ASCII files to initialize MARMITES..."
gridMETEO = cMF.cPROCESS.inputEsriAscii(grid_fn = gridMETEO_fn, datatype = int)
gridSOIL = cMF.cPROCESS.inputEsriAscii(grid_fn = gridSOIL_fn, datatype = int)
gridSOILthick = cMF.cPROCESS.inputEsriAscii(grid_fn = gridSOILthick_fn, datatype = float)
gridSsurfhmax = cMF.cPROCESS.inputEsriAscii(grid_fn = gridSsurfhmax_fn, datatype = float)
gridSsurfw = cMF.cPROCESS.inputEsriAscii(grid_fn = gridSsurfw_fn, datatype = float)
if irr_yn == 1:
    gridIRR = cMF.cPROCESS.inputEsriAscii(grid_fn = gridIRR_fn, datatype = int)
    if gridIRR.max() > NFIELD:
        cUTIL.ErrorExit('\nFATAL ERROR!\nThere is more fields in the asc file than in the MMsurf file!')

# READ input time series and parameters
if irr_yn == 0:
    gridVEGarea, RF_veg_zoneSP, E0_zonesSP, PT_veg_zonesSP, RFe_veg_zonesSP, LAI_veg_zonesSP, PE_zonesSP = cMF.cPROCESS.inputSP(                       NMETEO                   = NMETEO,
                                NVEG                     = NVEG,
                                NSOIL                    = NSOIL,
                                nper                     = cMF.nper,
                                inputZON_SP_RF_veg_fn    = cMF.inputZON_SP_RF_veg_fn,
                                inputZON_SP_RFe_veg_fn   = cMF.inputZON_SP_RFe_veg_fn,
                                inputZON_SP_LAI_veg_fn   = cMF.inputZON_SP_LAI_veg_fn,
                                inputZON_SP_PT_fn        = cMF.inputZON_SP_PT_fn,
                                inputZON_SP_PE_fn        = cMF.inputZON_SP_PE_fn,
                                inputZON_SP_E0_fn        = cMF.inputZON_SP_E0_fn,
                                )
else:
    gridVEGarea, RF_veg_zoneSP, E0_zonesSP, PT_veg_zonesSP, RFe_veg_zonesSP, LAI_veg_zonesSP, PE_zonesSP, RF_irr_zoneSP, RFe_irr_zoneSP, PT_irr_zonesSP, crop_irr_SP = cMF.cPROCESS.inputSP(
                                NMETEO                   = NMETEO,
                                NVEG                     = NVEG,
                                NSOIL                    = NSOIL,
                                nper                     = cMF.nper,
                                inputZON_SP_RF_veg_fn    = cMF.inputZON_SP_RF_veg_fn,
                                inputZON_SP_RFe_veg_fn   = cMF.inputZON_SP_RFe_veg_fn,
                                inputZON_SP_LAI_veg_fn   = cMF.inputZON_SP_LAI_veg_fn,
                                inputZON_SP_PT_fn        = cMF.inputZON_SP_PT_fn,
                                inputZON_SP_PE_fn        = cMF.inputZON_SP_PE_fn,
                                inputZON_SP_E0_fn        = cMF.inputZON_SP_E0_fn,
                                NFIELD                   = NFIELD,
                                inputZON_SP_RF_irr_fn    = cMF.inputZON_SP_RF_irr_fn,
                                inputZON_SP_RFe_irr_fn   = cMF.inputZON_SP_RFe_irr_fn,
                                inputZON_SP_PT_irr_fn    = cMF.inputZON_SP_PT_irr_fn,
                                input_SP_crop_irr_fn     = cMF.input_SP_crop_irr_fn
                                )

# SOIL PARAMETERS
_nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _S_ini, _Ks = cMF.cPROCESS.inputSoilParam(SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
_nslmax = max(_nsl)
for l in range(NSOIL):
    _slprop[l] = np.asarray(_slprop[l])

# compute thickness, top and bottom elevation of each soil layer
cMF.elev = np.ma.masked_values(np.asarray(cMF.elev), cMF.hnoflo, atol = 0.09)
cMF.top = np.ma.masked_values(cMF.elev, cMF.hnoflo, atol = 0.09) - gridSOILthick
cMF.botm = np.asarray(cMF.botm)
for l in range(cMF.nlay):
    cMF.botm[l,:,:] = np.ma.masked_values(cMF.botm[l,:,:], cMF.hnoflo, atol = 0.09) - gridSOILthick
cMF.botm = np.ma.masked_values(cMF.botm, cMF.hnoflo, atol = 0.09)
botm_l0 = np.asarray(cMF.botm)[:,:,0]

# create MM array
h5_MM_fn = os.path.join(MM_ws,'_h5_MM.h5')
# indexes of the HDF5 output arrays
index = {'iRF':0, 'iPT':1, 'iPE':2, 'iRFe':3, 'iSsurf':4, 'iRo':5, 'iEXF':6, 'iEsurf':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSsurf':13, 'iETg':14, 'iETsoil':15, 'iSsoil_pc':16, 'idSsoil':17, 'iinf':18, 'iHEADScorr':19, 'idgwt':20, 'iuzthick':21}
index_S = {'iEsoil':0, 'iTsoil':1,'iSsoil_pc':2, 'iRp':3, 'iRexf':4, 'idSsoil':5, 'iSsoil':6, 'iSAT':7, 'iMB_l':8}

# READ observations time series (heads and soil moisture)
print "\nReading observations time series (hydraulic heads and soil moisture)..."
obs, obs_list, obs_catch, obs_catch_list = cMF.cPROCESS.inputObs(
                              inputObs_fn      = inputObs_fn,
                              inputObsHEADS_fn = inputObsHEADS_fn,
                              inputObsSM_fn    = inputObsSM_fn,
                              inputObsRo_fn    = inputObsRo_fn,
                              inputDate        = cMF.inputDate,
                              _nslmax          = _nslmax,
                              nlay             = cMF.nlay
                              )
i = []
j = []
lay = []
lbl = []
for o_ref in obs_list:
    for o in obs.keys():
        if o == o_ref:
            i.append(obs.get(o)['i']+1)
            j.append(obs.get(o)['j']+1)
            lay.append(obs.get(o)['lay'])
            lbl.append(obs.get(o)['lbl'])
obs4map = [lbl, i, j, lay]
# To write MM output in a txt file
Ssoil_str   = ''
Ssoilpc_str = ''
dSsoil_str  = ''
Rp_str   = ''
Rexf_str = ''
Esoil_str   = ''
Tsoil_str   = ''
Smeasout = ''
MB_str   = ''
for l in range(_nslmax):
    Ssoil_str = Ssoil_str + 'Ssoil_l' + str(l+1) + ','
    Ssoilpc_str = Ssoilpc_str + 'Ssoilpc_l' + str(l+1) + ','
    dSsoil_str = dSsoil_str + 'dSsoil_l' + str(l+1) + ','
    Esoil_str = Esoil_str + 'Esoil_l' + str(l+1) + ','
    Tsoil_str = Tsoil_str + 'Tsoil_l' + str(l+1) + ','
    Rp_str = Rp_str + 'Rp_l' + str(l+1) + ','
    Rexf_str = Rexf_str + 'Rexf_l' + str(l+1) + ','
    MB_str = MB_str + 'MB_l' + str(l+1) + ','
    Smeasout = Smeasout + 'Smeas_' + str(l+1) + ','
header='Date,MF_SP,veg_crop,RF,E0,PT,PE,RFe,I,' + Esoil_str + Tsoil_str + 'Eg,Tg,ETg,WEL_MF,Esurf,' + Ssoil_str + Ssoilpc_str + dSsoil_str + 'dSsurf,Ssurf,Ro,Romeas,GW_EXF,' + Rp_str + Rexf_str + 'R_MF,hSATFLOW,hMF,hMFcorr,hmeas,dgwt,' + Smeasout + MB_str + 'MB\n'
if cMF.uzf_yn == 1:
    cMF.uzf_obs(obs = obs)

# EXPORT INPUT MAPS
if plt_input == 1:
    ibound = np.asarray(cMF.ibound)
    print'\n##############'
    print 'Exporting input maps...'
    i_lbl = 1
    T = np.zeros((cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
    if cMF.nlay > 1:
        T = np.asarray(cMF.hk_actual) * np.asarray(cMF.thick)
    else:
        T[0,:,:] = np.asarray(cMF.hk_actual) * np.asarray(cMF.thick[:,:,0])
    top_tmp = np.zeros((cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
    top_tmp[0,:,:] = cMF.top
    for l in range(1,cMF.nlay):
        top_tmp[l,:,:] = cMF.botm[l-1,:,:]
    lst = [cMF.elev, top_tmp, cMF.botm, cMF.thick, np.asarray(cMF.strt), gridSOILthick, gridSsurfhmax, gridSsurfw, np.asarray(cMF.hk_actual), T, np.asarray(cMF.ss_actual), np.asarray(cMF.sy_actual), np.asarray(cMF.vka_actual)]
    lst_lbl = ['elev', 'top', 'botm', 'thick', 'strt', 'gridSOILthick', 'gridSsurfhmax', 'gridSsurfw', 'hk', 'T', 'Ss', 'Sy', 'vka']
    lst_lblCB = ['Elev.', 'Aq. top - $top$', 'Aq. bot. - $botm$', 'Aq. thick.', 'Init. heads - $strt$', 'Soil thick.', 'Max. stream heigth', 'Stream width', 'Horizontal hydraulic cond. - $hk$', 'Transmissivity - $T$','Specific storage - $S_s$', 'Specific yield - $S_y$', 'Vertical hydraulic cond. - $vka$']
    if cMF.drn_yn == 1:
        lst.append(cMF.drn_cond_array)
        lst_lbl.append('drn_cond')
        lst_lblCB.append('Drain cond.')
        lst.append(cMF.drn_elev_array)
        lst_lbl.append('drn_elev')
        lst_lblCB.append('Drain elev.')
    if cMF.ghb_yn == 1:
        lst.append(cMF.ghb_cond_array)
        lst_lbl.append('ghb_cond')
        lst_lblCB.append('GHB cond.')
        lst.append(cMF.ghb_head_array)
        lst_lbl.append('ghb_head')
        lst_lblCB.append('GHB head')
    if cMF.uzf_yn == 1:
        lst.append(np.asarray(cMF.eps_actual))
        lst_lbl.append('eps')
        lst_lblCB.append('Epsilon - $eps$')
        lst.append(np.asarray(cMF.thts_actual))
        lst_lbl.append('thts')
        lst_lblCB.append('Sat. water content - $thts$')
        if cMF.vks_actual != 0.0:
            lst.append(np.asarray(cMF.vks_actual))
            lst_lbl.append('vks')
            lst_lblCB.append('Sat. vert. hydraulic conductivity - $vks$')
        try:
            lst.append(np.asarray(cMF.thti_actual))
            lst_lbl.append('thti')
            lst_lblCB.append('Initial water content - $thti$')
        except:
            pass
        try:
            lst.append(np.asarray(cMF.thtr_actual))
            lst_lbl.append('thtr')
            lst_lblCB.append('Residual water content - $thtr$')
        except:
            pass
    elev_max = []
    elev_min = []
    for e in [cMF.elev, cMF.top, cMF.botm]:
        elev_max.append(np.max(e))
        elev_min.append(np.min(e))
    elev_max = max(elev_max)
    elev_min = min(elev_min)
    for i, l in enumerate(lst):
        V = np.zeros((1, cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        Vmax = []
        Vmin = []
        mask_tmp = np.zeros((cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        for L in range(cMF.nlay):
            if l.shape == (cMF.nlay, cMF.nrow, cMF.ncol) or l.shape == (cMF.nrow, cMF.ncol):
                try:
                    V[0,L,:,:] = l[L,:,:]
                    mask_tmp[L,:,:] = mask[L]
                    Vmax.append(np.ma.max(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], mask[L]), cMF.hnoflo, atol = 0.09)))
                    Vmin.append(np.ma.min(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], mask[L]), cMF.hnoflo, atol = 0.09)))
                    nplot = cMF.nlay
                except:
                    V[0,L,:,:] = l
                    mask_tmp[L,:,:] = maskAllL
                    Vmax.append(np.ma.max(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], maskAllL), cMF.hnoflo, atol = 0.09)))
                    Vmin.append(np.ma.min(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], maskAllL), cMF.hnoflo, atol = 0.09)))
                    nplot = 1
            elif l.shape != ():
                V[0,L,:,:] = l[L]*ibound[L,:,:]
                mask_tmp[L,:,:] = mask[L]
                Vmax.append(np.ma.max(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], mask[L]), cMF.hnoflo, atol = 0.09)))
                Vmin.append(np.ma.min(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], mask[L]), cMF.hnoflo, atol = 0.09)))
                nplot = cMF.nlay
            else:
                V[0,L,:,:] = l*np.invert(maskAllL)
                mask_tmp[L,:,:] = maskAllL
                Vmax.append(np.ma.max(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], maskAllL), cMF.hnoflo, atol = 0.09)))
                Vmin.append(np.ma.min(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], maskAllL), cMF.hnoflo, atol = 0.09)))
                nplot = 1
        if lst_lbl[i] == 'Ss' or lst_lbl[i] == 'hk' or lst_lbl[i] == 'T' or lst_lbl[i] == 'drn_cond' or lst_lbl[i] == 'ghb_cond':
            fmt = '%5.e'
        elif lst_lbl[i] == 'Sy' or lst_lbl[i] == 'thts' or lst_lbl[i] == 'thti' or lst_lbl[i] == 'thtr' or lst_lbl[i] == 'gridSOILthick' or lst_lbl[i] == 'gridSsurfhmax' or lst_lbl[i] == 'gridSsurfw':
            fmt = '%5.3f'
        else:
            fmt = '%5.1f'
        if lst_lbl[i] == 'Ss' or lst_lbl[i] == 'eps':
           CBlabel = lst_lblCB[i] + ' $([-])$'
        elif lst_lbl[i] == 'hk' or lst_lbl[i] == 'drn_cond' or lst_lbl[i] == 'ghb_cond':
            CBlabel = lst_lblCB[i] + ' $(%s/%s)$' % (lenuni_str, itmuni_str)
        elif lst_lbl[i] == 'T':
            CBlabel = lst_lblCB[i] + ' $(%s^{2}/%s)$' % (lenuni_str, itmuni_str)
        elif lst_lbl[i] == 'Sy' or lst_lbl[i] == 'thts' or lst_lbl[i] == 'thti' or lst_lbl[i] == 'thtr':
            CBlabel = lst_lblCB[i] + ' $(L^3/L^{-3})$'
        elif lst_lbl[i] == 'vka':
            CBlabel = lst_lblCB[i] + ' $([-]$ $or$ $%s/%s,$ $layvka = %s)$' % (lenuni_str, itmuni_str, cMF.layvka)
        else:
            CBlabel = lst_lblCB[i] + ' $(m)$'
        if lst_lbl[i] == 'elev' or lst_lbl[i] == 'top' or lst_lbl[i] == 'botm':
            Vmax = elev_max
            Vmin = elev_min
        else:
            Vmax = np.ma.max(Vmax)
            Vmin = np.ma.min(Vmin)
        MMplot.plotLAYER(timesteps = [0], Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = nplot, V = V,  cmap = plt.cm.gist_rainbow_r, CBlabel = CBlabel, msg = '', plt_title = 'IN_%03d_%s' % (i_lbl,lst_lbl[i]), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsMM, Vmax = [Vmax], Vmin = [Vmin], ntick = ntick, fmt = fmt, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)
        i_lbl += 1
    del V, lst, lst_lbl, nplot, Vmax, Vmin

    Vmax = 100.0
    Vmin = 0.0
    V = np.zeros((1, cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
    for v in range(NVEG):
        V[0,0,:,:] = gridVEGarea[v,:,:]
        V_lbl = 'veg%02d_%s' %(v+1, VegName[v])
        V_lblCB = 'Frac. area of veg. #%02d - $%s$ $(\%%)$' %(v+1, VegName[v])
        for L in range(cMF.nlay):
            mask_tmp[L,:,:] = maskAllL
        #print V_lbl
        MMplot.plotLAYER(timesteps = [0], Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.gist_rainbow_r, CBlabel = V_lblCB, msg = '', plt_title = 'IN_%03d_%s'% (i_lbl,V_lbl), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = False, Vmax = [Vmax], Vmin = [Vmin], ntick = ntick, fmt = '%5.1f', points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)
        i_lbl += 1
    del V, V_lbl, Vmax, Vmin

    lst = [ibound, gridSOIL, gridMETEO]
    lst_lbl = ['ibound', 'gridSOIL', 'gridMETEO']
    lst_lblCB = ['MF cell type - $ibound$', 'Soil type', 'Meteo. zone']
    if irr_yn == 1:
        lst.append(gridIRR)
        lst_lbl.append('gridIRR')
        lst_lblCB.append('Irrigation plots')
    if cMF.uzf_yn == 1:
        lst.append(np.asarray(cMF.iuzfbnd))
        lst_lbl.append('iuzfbnd')
        lst_lblCB.append('Recharge layer - $iuzfbnd$')
    for i, l in enumerate(lst):
        V = np.zeros((1, cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        Vmax = []
        Vmin = []
        for L in range(cMF.nlay):
            try:
                V[0,L,:,:] = l[L,:,:]
                mask_tmp[L,:,:] = mask[L]
                Vmax.append(np.ma.max(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], maskAllL), cMF.hnoflo, atol = 0.09)))
                Vmin.append(np.ma.min(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], maskAllL), cMF.hnoflo, atol = 0.09)))
                nplot = cMF.nlay
            except:
                V[0,L,:,:] = l
                mask_tmp[L,:,:] = maskAllL
                Vmax.append(np.ma.max(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], maskAllL), cMF.hnoflo, atol = 0.09)))
                Vmin.append(np.ma.min(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], maskAllL), cMF.hnoflo, atol = 0.09)))
                nplot = 1
        Vmax = np.ma.max(Vmax)
        Vmin = np.ma.min(Vmin)
        MMplot.plotLAYER(timesteps = [0], Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = nplot, V = V,  cmap = plt.cm.gist_rainbow_r, CBlabel = lst_lblCB[i], msg = '', plt_title = 'IN_%03d_%s'% (i_lbl,lst_lbl[i]), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = False, Vmax = [Vmax], Vmin = [Vmin], ntick = ntick, fmt = '%2.2f', points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)
        i_lbl += 1
    del V, lst, lst_lbl, lst_lblCB, nplot, Vmax, Vmin

# #############################
# ### 1st MODFLOW RUN with initial user-input recharge
# #############################
print'\n##############'
if MF_yn == 1 :
    timestartMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    print 'MODFLOW RUN (initial user-input fluxes)\n'
    if verbose == 0:
        print '\n--------------'
        sys.stdout = s
        report.close()
        s = sys.stdout
        report = open(report_fn, 'a')
        sys.stdout = report
    cMF.runMF(finf_MM = (h5_MM_fn, 'finf'), wel_MM = (h5_MM_fn, 'ETg'), report = report, verbose = verbose, chunks = chunks, numDays = numDays)
    timeendMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    durationMFtmp =  timeendMF - timestartMF
    durationMF +=  durationMFtmp
    print 'MF run time: %02.fmn%02.fs' % (int(durationMFtmp*24.0*60.0), (durationMFtmp*24.0*60.0-int(durationMFtmp*24.0*60.0))*60)
    del durationMFtmp

if os.path.exists(cMF.h5_MF_fn):
    print 'Reading MF fluxes...'
    try:
        h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
        # heads format is : timestep, nlay, nrow, ncol
        # cbc format is: (kstp), kper, textprocess, nlay, nrow, ncol
        cbc_nam = []
        cbc_uzf_nam = []
        for c in h5_MF['cbc_nam']:
            cbc_nam.append(c.strip())
        if cMF.uzf_yn == 1:
            for c in h5_MF['cbc_uzf_nam']:
                cbc_uzf_nam.append(c.strip())
        elif cMF.rch_yn == 1:
            #imfRCH = cbc_nam.index('RECHARGE')
            cUTIL.ErrorExit('\nFATAL ERROR!\nMM has to be run together with the UZF1 package of MODFLOW-NWT, thus the RCH package should be desactivacted!\nExisting MM.')
        cbc_uzf_nam_tex = [0]*len(cbc_uzf_nam)
        imfSTO = cbc_nam.index('STORAGE')
        imfFLF = cbc_nam.index('FLOW LOWER FACE')
        if cMF.ghb_yn == 1:
            imfGHB = cbc_nam.index('HEAD DEP BOUNDS')
        if cMF.drn_yn == 1:
            imfDRN = cbc_nam.index('DRAINS')
        if cMF.wel_yn == 1:
            imfWEL = cbc_nam.index('WELLS')
        else:
            print '\nWARNING!\nThe WEL package should be active to take into account ETg!'
        if cMF.uzf_yn == 1:
            imfEXF   = cbc_uzf_nam.index('SURFACE LEAKAGE')
            imfRCH   = cbc_uzf_nam.index('UZF RECHARGE')
        if MMsoil_yn == 0:
            h5_MF.close()
    except:
        cUTIL.ErrorExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.')

# #############################
# 2nd phase : MM/MF loop #####
# #############################
h_diff_surf = None
if MMsoil_yn > 0:
    durationMMsoil = 0.0
    h_pSP = 0
    h_pSP_all = 0
    LOOP = 0
    endloop = 0
    LOOPlst = [LOOP]
    h_diff = [1000]
    h_diff_log = [1]
    h_diff_all = [1000]
    h_diff_all_log = [1]
    plt_ConvLoop_fn = os.path.join(MM_ws_out, '__plt_MM_MF_ConvLoop.png')
    # Create HDF5 arrays to store MARMITES output
    try:
        h5_MM = h5py.File(h5_MM_fn, 'w')
    except:
        os.remove(h5_MM_fn)
        h5_MM = h5py.File(h5_MM_fn, 'w')
        print "WARNING! Previous h5_MM file corrupted!\nIt was deleted and a new one was created."
    # arrays for fluxes independent of the soil layering
    h5_MM.create_dataset(name = 'iMM', data = np.asarray(index.items()))
    if chunks == 1:
        h5_MM.create_dataset(name = 'MM', shape = (sum(cMF.perlen),cMF.nrow,cMF.ncol,len(index)), dtype = np.float, chunks = (1,cMF.nrow,cMF.ncol,len(index)),  compression = 'gzip', compression_opts = 5, shuffle = True)
    else:
        h5_MM.create_dataset(name = 'MM', shape = (sum(cMF.perlen),cMF.nrow,cMF.ncol,len(index)), dtype = np.float)
    # arrays for fluxes in each soil layer
    h5_MM.create_dataset(name = 'iMM_S', data = np.asarray(index_S.items()))
    if chunks == 1:
        h5_MM.create_dataset(name = 'MM_S', shape = (sum(cMF.perlen),cMF.nrow,cMF.ncol,_nslmax,len(index_S)), dtype = np.float, chunks = (1,cMF.nrow,cMF.ncol,_nslmax,len(index_S)),  compression = 'gzip', compression_opts = 5, shuffle = True)
    else:
        h5_MM.create_dataset(name = 'MM_S', shape = (sum(cMF.perlen),cMF.nrow,cMF.ncol,_nslmax,len(index_S)), dtype = np.float)
    # arrays to compute net recharge to be exported to MF
    if chunks == 1:
        h5_MM.create_dataset(name = 'finf', shape = (cMF.nper,cMF.nrow,cMF.ncol), dtype = np.float, chunks = (1,cMF.nrow,cMF.ncol),  compression = 'gzip', compression_opts = 5, shuffle = True)
        h5_MM.create_dataset(name = 'ETg', shape = (cMF.nper,cMF.nrow,cMF.ncol), dtype = np.float, chunks = (1,cMF.nrow,cMF.ncol),  compression = 'gzip', compression_opts = 5, shuffle = True)
    else:
        h5_MM.create_dataset(name = 'finf', shape = (cMF.nper,cMF.nrow,cMF.ncol), dtype = np.float)
        h5_MM.create_dataset(name = 'ETg', shape = (cMF.nper,cMF.nrow,cMF.ncol), dtype = np.float)

    # #############################
    # ###  CONVERGENCE LOOP   #####
    # #############################

    while (abs(h_diff[LOOP]) > convcrit or abs(h_diff_all[LOOP]) > convcritmax) and LOOP <= ccnum :
        if LOOP == 0:
            print '\n##############\nCONVERGENCE LOOP %d (initialization)\n##############' % (LOOP)
        else:
            print '\n##############\nCONVERGENCE LOOP %d/%d\n##############' % (LOOP, ccnum)
        h_MF_average = 0.0
        timestartMMloop = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        # ###########################
        # ###  MARMITES INPUT #######
        # ###########################
        print'\n##############'
        print 'MARMITESsoil RUN'
        # SOIL PARAMETERS
        _nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _S_ini, _Ks = cMF.cPROCESS.inputSoilParam(SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
        _nslmax = max(_nsl)
        for l in range(NSOIL):
            _slprop[l] = np.asarray(_slprop[l])
        if LOOP > 0:
            h5_MM = h5py.File(h5_MM_fn)
        # ###############
        # # main loop: calculation of soil water balance in each cell-grid for each time step inside each stress period
#        t0=0
        print '\nComputing...'
        if irr_yn == 0:
            MM_SOIL.run(_nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _S_ini, botm_l0, _Ks,
                          gridSOIL, gridSOILthick, cMF.elev*1000.0, gridMETEO,
                          index, index_S, gridSsurfhmax, gridSsurfw,
                          RF_veg_zoneSP, E0_zonesSP, PT_veg_zonesSP, RFe_veg_zonesSP, PE_zonesSP, gridVEGarea,
                          LAI_veg_zonesSP, Zr, kT_min, kT_max, kT_n, NVEG,
                          cMF, conv_fact, h5_MF, h5_MM, irr_yn
                          )
        else:
            MM_SOIL.run(_nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _S_ini, botm_l0, _Ks,
                          gridSOIL, gridSOILthick, cMF.elev*1000.0, gridMETEO,
                          index, index_S, gridSsurfhmax, gridSsurfw,
                          RF_veg_zoneSP, E0_zonesSP, PT_veg_zonesSP, RFe_veg_zonesSP, PE_zonesSP, gridVEGarea,
                          LAI_veg_zonesSP, Zr, kT_min, kT_max, kT_n, NVEG,
                          cMF, conv_fact, h5_MF, h5_MM, irr_yn,
                          RF_irr_zoneSP, PT_irr_zonesSP, RFe_irr_zoneSP,
                          crop_irr_SP, gridIRR,
                          Zr_c, kT_min_c, kT_max_c, kT_n_c, NCROP
                          )

        # CHECK MM amd MF CONVERG.
        h_MF_m = np.ma.masked_values(np.ma.masked_values(h5_MF['heads'], cMF.hdry, atol = 1E+25), cMF.hnoflo, atol = 0.09)
        h5_MF.close()
        h_MF_average = np.ma.average(h_MF_m)
        h_diff.append(h_MF_average - h_pSP)
        h_diff_surf = h_MF_m - h_pSP_all
        h_diff_all_max = np.ma.max(h_diff_surf)
        h_diff_all_min = np.ma.min(h_diff_surf)
        if abs(h_diff_all_max)>abs(h_diff_all_min):
            h_diff_all.append(h_diff_all_max)
        else:
            h_diff_all.append(h_diff_all_min)
        del h_diff_all_max, h_diff_all_min
        LOOPlst.append(LOOP)
        LOOP += 1
        h_pSP = h_MF_average
        h_pSP_all = h_MF_m
        del h_MF_m
        if np.absolute(h_diff[LOOP])>0.0:
            h_diff_log.append(np.log10(np.absolute(h_diff[LOOP])))
            h_diff_all_log.append(np.log10(np.absolute(h_diff_all[LOOP])))
        else:
            h_diff_log.append(np.log10(convcrit))
            h_diff_all_log.append(np.log10(convcritmax))

        msg_end_loop = []
        msg_end_loop.append('Average heads:\n%.3f m' % h_MF_average)
        if LOOP > 1:
            msg_end_loop.append('Heads diff. from previous conv. loop: %.3f m' % h_diff[LOOP])
            msg_end_loop.append('Maximum heads difference:             %.3f m' % h_diff_all[LOOP])
        if h_MF_average == 0.0 or str(h_diff[LOOP]) == 'nan':
            print '\Warning!\nModel with DRY cells or NaN values!'
        elif abs(h_diff[LOOP]) < convcrit and abs(h_diff_all[LOOP]) < convcritmax:
            msg_end_loop.append('Successful convergence between MARMITES and MODFLOW!\n(Conv. criterion = %.4f and conv. crit. max. = %.4f)' % (convcrit, convcritmax))
            endloop += 1
        elif LOOP>ccnum:
            msg_end_loop.append('No convergence between MARMITES and MODFLOW!\n(Conv. criterion = %.4f and conv. crit. max. = %.4f)' % (convcrit, convcritmax))
            endloop += 1
        for txt in msg_end_loop:
            print txt
        del h_MF_average

        timeendMMloop = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        durationMMloop = timeendMMloop-timestartMMloop
        print '\nMM run time: %02.fmn%02.fs' % (int(durationMMloop*24.0*60.0), (durationMMloop*24.0*60.0-int(durationMMloop*24.0*60.0))*60)
        durationMMsoil += durationMMloop
        print '%s'% mpl.dates.DateFormatter.format_data(fmt_DH, mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat()))

        # MODFLOW RUN with MM-computed recharge
        timestartMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        print'\n##############'
        if endloop < 1:
            print 'MODFLOW RUN (MARMITES fluxes)'
        else:
            print 'MODFLOW RUN (MARMITES fluxes after conv. loop)'
        if verbose == 0:
            print '\n--------------'
            sys.stdout = s
            report.close()
            s = sys.stdout
            report = open(report_fn, 'a')
            sys.stdout = report
        cMF.runMF(finf_MM = (h5_MM_fn, 'finf'), wel_MM = (h5_MM_fn, 'ETg'), report = report, verbose = verbose, chunks = chunks, numDays = numDays)
        try:
            h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
        except:
            cUTIL.ErrorExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.')
        timeendMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        durationMFtmp =  timeendMF-timestartMF
        durationMF +=  durationMFtmp
        print '\nMF run time: %02.fmn%02.fs' % (int(durationMFtmp*24.0*60.0), (durationMFtmp*24.0*60.0-int(durationMFtmp*24.0*60.0))*60)
        del durationMFtmp
    h5_MF.close()
    # #############################
    # ###  END CONVERGENCE LOOP ###
    # #############################

    # export loop plot
    print'\n##############'
    print 'Exporting plot of the convergence loop...'
    fig = plt.figure()
    fig.suptitle('Convergence loop plot between MM and MF based on heads differences.\nOrange: average heads for the whole model.\nGreen: maximun heads difference observed in the model (one cell)', fontsize=10)
    if LOOP>0:
        ax1=fig.add_subplot(3,1,1)
        plt.setp(ax1.get_xticklabels(), fontsize=8)
        plt.setp(ax1.get_yticklabels(), fontsize=8)
        ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
        plt.ylabel('h_diff [m]', fontsize=10, horizontalalignment = 'center')
        plt.grid(True)
        plt.plot(LOOPlst[1:], h_diff[1:], linestyle='-', marker='o', markersize=5, c = 'orange', markerfacecolor='orange', markeredgecolor='red')
        plt.plot(LOOPlst[1:], h_diff_all[1:], linestyle='-', marker='o', markersize=5, c = 'green', markerfacecolor='green', markeredgecolor='blue')

    if LOOP>1:
        ax2=fig.add_subplot(3,1,2, sharex = ax1)
        plt.setp(ax2.get_xticklabels(), fontsize=8)
        plt.setp(ax2.get_yticklabels(), fontsize=8)
        plt.plot(LOOPlst[2:], h_diff[2:], linestyle='-', marker='o', markersize=5, c = 'orange', markerfacecolor='orange', markeredgecolor='red')
        plt.plot(LOOPlst[2:], h_diff_all[2:], linestyle='-', marker='o', markersize=5, c = 'green', markerfacecolor='green', markeredgecolor='blue')
        ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
        plt.ylabel('h_diff [m]', fontsize=10, horizontalalignment = 'center')
    #        plt.ylabel.Text.position(0.5, -0.5)
        plt.grid(True)

        ax3=fig.add_subplot(3,1,3, sharex = ax1)
        plt.setp(ax3.get_xticklabels(), fontsize=8)
        plt.setp(ax3.get_yticklabels(), fontsize=8)
        ax3.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
        plt.ylabel('log(abs(h_diff)) [log(m)]', fontsize=10, horizontalalignment = 'center')
        plt.grid(True)
        plt.xlabel('loop', fontsize=10)
        plt.plot(LOOPlst[2:], h_diff_log[2:], linestyle='-', marker='o', markersize=5, c = 'orange', markerfacecolor='orange', markeredgecolor='red')
        plt.plot(LOOPlst[2:], h_diff_all_log[2:], linestyle='-', marker='o', markersize=5, c = 'green', markerfacecolor='green', markeredgecolor='blue')

        ax2.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%2d'))
        ax3.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%2d'))

        plt.xlim(0,LOOP-1)
        ax1.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%2d'))
        ax1.xaxis.set_ticks(LOOPlst[1:])

    plt.savefig(plt_ConvLoop_fn)
    plt.cla()
    plt.clf()
    plt.close('all')
    del fig, LOOPlst, h_diff, h_diff_log, h_pSP_all

# #############################
# 3rd phase : export results #####
# #############################

timestartExport = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())

print '\n##############\nMARMITES exporting...'

del RF_veg_zoneSP
del E0_zonesSP
del PT_veg_zonesSP
del RFe_veg_zonesSP
del PE_zonesSP
del gridSsurfhmax
del gridSsurfw

# reorganizing MF output in daily data
if MF_yn == 1 and isinstance(cMF.h5_MF_fn, str):
    print '\nConverting MODFLOW output into daily time step...'
    try:
        h5_MF = h5py.File(cMF.h5_MF_fn)
    except:
        cUTIL.ErrorExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.')
    cMF.cPROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'STO_d', conv_fact = conv_fact, index = imfSTO)
    cMF.cPROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'FLF_d', conv_fact = conv_fact, index = imfFLF)
    if cMF.drn_yn == 1:
        cMF.cPROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'DRN_d', conv_fact = conv_fact, index = imfDRN)
    if cMF.wel_yn == 1:
        cMF.cPROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'WEL_d', conv_fact = conv_fact, index = imfWEL)
    if cMF.ghb_yn == 1:
        cMF.cPROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'GHB_d', conv_fact = conv_fact, index = imfGHB)
    cMF.cPROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'heads', ds_name_new = 'heads_d', conv_fact = conv_fact)
    if cMF.uzf_yn == 1:
        cMF.cPROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc_uzf', ds_name_new = 'RCH_d', conv_fact = conv_fact, index = imfRCH)
        cMF.cPROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc_uzf', ds_name_new = 'EXF_d', conv_fact = conv_fact, index = imfEXF)
    elif cMF.rch_yn == 1:
        cMF.cPROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'RCH_d', conv_fact = conv_fact, index = imfRCH)
    h5_MF.close()

if MMsoil_yn == 1 and isinstance(h5_MM_fn, str):
    try:
        h5_MM = h5py.File(h5_MM_fn)
    except:
        cUTIL.ErrorExit('\nFATAL ERROR!\nInvalid MARMITES HDF5 file. Run MARMITES and/or MODFLOW again.')
    cMF.cPROCESS.procMM(cMF = cMF, h5_MM = h5_MM, ds_name = 'finf', ds_name_new = 'finf_d')
    cMF.cPROCESS.procMM(cMF = cMF, h5_MM = h5_MM, ds_name = 'ETg', ds_name_new = 'ETg_d')
    h5_MM.close()

SP_d = np.ones(sum(cMF.perlen), dtype = int)
t = 0
for n in range(cMF.nper):
    for x in range(cMF.perlen[n]):
        SP_d[t] = n + 1
        t += 1

if os.path.exists(cMF.h5_MF_fn):
    index_cbc = [imfSTO]
    index_cbc.append(imfFLF)
    if cMF.rch_yn == 1:
       index_cbc.append(imfRCH)
    if cMF.drn_yn == 1:
        index_cbc.append(imfDRN)
    if cMF.wel_yn == 1:
        index_cbc.append(imfWEL)
    if cMF.ghb_yn == 1:
        index_cbc.append(imfGHB)
    if cMF.uzf_yn == 1:
        index_cbc_uzf = [imfRCH, imfEXF]
else:
    cbc_DRN = cbc_STO = cbc_RCH = cbc_WEL = np.zeros((sum(cMF.perlen), cMF.nlay, cMF.nrow, cMF.ncol))
    imfDRN = imfSTO = imfFLF = imfRCH = imfWEL = 0

if h_diff_surf != None:
    h_diff_n = 0
    for n in range(cMF.nper):
        for l in range(cMF.nlay):
            for r, c in enumerate(h_diff_surf[n,:,:,l]):
                try:
                    list(c).index(h_diff_all[LOOP])
                    h_diff_n = n
                    break
                except:
                    pass
    del n, l, r, c
    V = []
    Vmax = []
    Vmin = []
    V = np.zeros((1, cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
    mask_tmp = np.zeros((cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
    for L in range(cMF.nlay):
        V[0,L,:,:] = h_diff_surf[h_diff_n,L,:,:]
        mask_tmp[L,:,:] = mask[L]
        Vmax.append(np.ma.max(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], maskAllL), cMF.hnoflo, atol = 0.09)))
        Vmin.append(np.ma.min(np.ma.masked_values(np.ma.masked_array(V[0,L,:,:], maskAllL), cMF.hnoflo, atol = 0.09)))
    Vmax = np.ma.max(Vmax) #float(np.ceil(max(Vmax)))
    Vmin = np.ma.min(Vmin) #float(np.floor(min(Vmin)))
    # TODO JD and Date are not correct since h_diff_n is # stress periods and not # of days (same in the plots of MF and MM)
    MMplot.plotLAYER(timesteps = [h_diff_n], Date = [cMF.inputDate[h_diff_n]], JD = [cMF.JD[h_diff_n]], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = ('(m)'), msg = 'no value', plt_title = 'HEADSmaxdiff_ConvLoop', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = [Vmax], Vmin = [Vmin], contours = ctrsMF, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo, pref_plt_title = '__sp_plt')
    del h_diff_n, V, Vmin, Vmax

# exporting sm computed by MM for PEST (smp format)
if os.path.exists(h5_MM_fn):
    try:
        h5_MM = h5py.File(h5_MM_fn, 'r')
    except:
        cUTIL.ErrorExit('\nFATAL ERROR!\nInvalid MARMITES HDF5 file. Run MARMITES and/or MODFLOW again.')
    outPESTsmMM = open(os.path.join(MF_ws,'sm_MM4PEST.smp'), 'w')
    for o_ref in obs_list:
        for o in obs.keys():
            if o == o_ref:
                i = obs.get(o)['i']
                j = obs.get(o)['j']
                l = obs.get(o)['lay']
                obs_SM = obs.get(o)['obs_SM']
                if obs.get(o)['obs_sm_yn'] == 1:
                    cMF.cPROCESS.smMMname.append(o)
                MM_S = h5_MM['MM_S'][:,i,j,:,:]
                cMF.cPROCESS.ExportResultsMM4PEST(i, j, cMF.inputDate, _nslmax, MM_S, index_S, obs_SM, o)
    # write PEST smp file with MM output
    inputFile = cUTIL.readFile(MM_ws,inputObs_fn)
    ind = 0
    for i in range(len(inputFile)):
        line = inputFile[i].split()
        name = line[0]
        for j in cMF.cPROCESS.smMMname:
            if j == name:
                for l in cMF.cPROCESS.smMM[ind]:
                    outPESTsmMM.write(l)
        ind += 1
    outPESTsmMM.close()

# computing max and min values in MF fluxes for plotting
hmax = []
hmin = []
hmaxMM = []
hminMM = []
hdiff = []
cbcmax_d = []
cbcmin_d = []
axefact = 1.05
if os.path.exists(cMF.h5_MF_fn):
    # TODO missing STOuz and FLF (however this is not very relevant since these fluxes should not be the bigger in magnitude)
    try:
        h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
    except:
        cUTIL.ErrorExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.')
    # STO
    cbc_STO = h5_MF['STO_d']
    cbcmax_d.append(-1*np.ma.max(cbc_STO))
    cbcmin_d.append(-1*np.ma.min(cbc_STO))
    del cbc_STO
    # RCH
    cbc_RCH = h5_MF['RCH_d']
    RCHmax = np.ma.max(cbc_RCH)
    cbcmax_d.append(RCHmax)
    RCHmin = np.ma.min(cbc_RCH)
    print '\nMaximum GW recharge (%.2f mm/day) observed at:' % RCHmax
    tRCHmax = -1
    if RCHmax> 0.0:
        for l in range(cMF.nlay):
            for row in range(cMF.nrow):
                for t,col in enumerate(cbc_RCH[:,l,row,:]):
                    try:
                        if plt_out_obs == 1 and obs_list[-1] <> 'PzRCHmax':
                            j = list(col).index(RCHmax)
                            x = cMF.delc[row]*row + xllcorner
                            y = cMF.delr[j]*j + yllcorner
                            obs['Rm'] = {'x':x,'y':y, 'i': row, 'j': j, 'lay': l, 'hi':999, 'h0':999, 'RC':999, 'STO':999, 'outpathname':os.path.join(MM_ws_out,'_0Rm_ts.txt'), 'obs_h':[], 'obs_h_yn':0, 'obs_SM':[], 'obs_sm_yn':0, 'obs_Ro':[], 'obs_Ro_yn':0}
                            obs_list.append('Rm')
                            obs4map[0].append('Rm')
                            obs4map[1].append(row)
                            obs4map[2].append(j)
                            obs4map[3].append(l)
                            del x, y
                        print 'row %d, col %d, layer %d and day %d (%s)' % (row + 1, list(col).index(RCHmax) + 1, l+1, t+1, mpl.dates.num2date(cMF.inputDate[t] + 1.0).isoformat()[:10])
                        tRCHmax = t
                    except:
                        pass
    if tRCHmax < 0:
        print 'WARNING!\nNo recharge max found!'
    del cbc_RCH
    RCHmax = np.ma.max(RCHmax) #float(np.ceil(np.ma.max(RCHmax)))  #
    RCHmin = np.ma.min(RCHmin) #float(np.floor(np.ma.min(RCHmin))) #
    # WEL
    if cMF.wel_yn == 1:
        cbc_WEL = h5_MF['WEL_d']
        cbcmax_d.append(np.ma.max(cbc_WEL))
        cbcmin_d.append(np.ma.min(cbc_WEL))
    # DRN
    if cMF.drn_yn == 1:
        cbc_DRN = h5_MF['DRN_d']
        DRNmax = np.ma.max(cbc_DRN)
        cbcmax_d.append(DRNmax)
        DRNmin = np.ma.min(cbc_DRN)
        cbcmin_d.append(DRNmin)
        del cbc_DRN
        DRNmax = -1.0*DRNmin
        DRNmin = 0.0
    # GHB
    if cMF.ghb_yn == 1:
        cbc_GHB = h5_MF['GHB_d']
        GHBmax = np.ma.max(cbc_GHB)
        cbcmax_d.append(GHBmax)
        GHBmin = np.ma.min(cbc_GHB)
        cbcmin_d.append(GHBmin)
        del cbc_GHB
        GHBmax = -1.0*GHBmin
        GHBmin = 0.0
    cbcmax_d = np.ma.max(cbcmax_d) #float(np.ceil(np.ma.max(cbcmax_d)))
    cbcmin_d = np.ma.min(cbcmin_d) #float(np.floor(np.ma.min(cbcmin_d)))
    # h
    h_MF_m = np.ma.masked_values(np.ma.masked_values(h5_MF['heads_d'], cMF.hnoflo, atol = 0.09), cMF.hdry, atol = 1E+25)
    hmaxMF = np.ma.max(h_MF_m[:,:,:,:].flatten())
    hminMF = np.ma.min(h_MF_m[:,:,:,:].flatten())
    GWTD = cMF.elev[np.newaxis,np.newaxis,:,:] - h_MF_m
    GWTDmax = np.ma.max(GWTD.flatten())
    GWTDmin = np.ma.min(GWTD.flatten())
    h5_MF.close()
    del GWTD
else:
    DRNmax = GHBmax = cbcmax_d = 1.0
    DRNmin = GHBmin = cbcmin_d = -1.0
    hmaxMF = -9999.9
    hminMF = 9999.9

if os.path.exists(h5_MM_fn):
    try:
        h5_MM = h5py.File(h5_MM_fn, 'r')
    except:
        cUTIL.ErrorExit('\nFATAL ERROR!\nInvalid MARMITES HDF5 file. Run MARMITES and/or MODFLOW again.')
    # h
    headscorr_m = np.ma.masked_values(np.ma.masked_values(h5_MM['MM'][:,:,:,19], cMF.hnoflo, atol = 0.09), cMF.hdry, atol = 1E+25)
    hcorrmax = np.ma.max(headscorr_m.flatten())
    hcorrmin = np.ma.min(headscorr_m.flatten())
    GWTDcorr = cMF.elev[np.newaxis,:,:] - headscorr_m
    GWTDcorrmax = np.ma.max(GWTDcorr.flatten())
    GWTDcorrmin = np.ma.min(GWTDcorr.flatten())
    if GWTDcorrmin < 0.0:
        GWTDcorrmin = 0.0
    h5_MM.close()
    del headscorr_m, GWTDcorr
else:
    hcorrmax = -9999.9
    hcorrmin = 9999.9

if obs != None:
    x = 0
    for o in obs.keys():
        i = obs.get(o)['i']
        j = obs.get(o)['j']
        l = obs.get(o)['lay']
        obs_h = obs.get(o)['obs_h']
        if os.path.exists(cMF.h5_MF_fn):
            hmaxMF_tmp = np.ma.max(h_MF_m[:,l,i,j].flatten())
            hminMF_tmp = np.ma.min(h_MF_m[:,l,i,j].flatten())
        else:
            hmaxMF_tmp = -9999.9
            hminMF_tmp = -9999.9
        if obs_h != []:
            npa_m_tmp = np.ma.masked_values(obs_h, cMF.hnoflo, atol = 0.09)
            hmaxMM = np.ma.max(npa_m_tmp.flatten())
            hminMM = np.ma.min(npa_m_tmp.flatten())
            del npa_m_tmp
        else:
            hmaxMM = -9999.9
            hminMM = 9999.9
        hmax.append(np.ma.max((hmaxMF_tmp, hmaxMM)))
        hmin.append(np.ma.min((hminMF_tmp, hminMM)))
        hdiff.append(hmax[x]-hmin[x])
        x += 1
    hdiff = np.ma.max(hdiff)
    del hmaxMF_tmp, hminMF_tmp
else:
    hdiff = 2000

# #################################################
# plot SOIL/GW ts and balance at the catchment scale
# #################################################
# RMSE list to plot
rmseSM = []
rsrSM = []
rSM = []
obslstSM = []
rmseHEADS = []
rsrHEADS = []
nseHEADS = []
nseSM = []
rHEADS = []
obslstHEADS = []
tTgmin = -1
if os.path.exists(h5_MM_fn):
    try:
        h5_MM = h5py.File(h5_MM_fn, 'r')
    except:
        cUTIL.ErrorExit('\nFATAL ERROR!\nInvalid MARMITES HDF5 file. Run MARMITES and/or MODFLOW again.')
    # indexes of the HDF5 output arrays
    #  index = {'iRF':0, 'iPT':1, 'iPE':2, 'iRFe':3, 'iSsurf':4, 'iRo':5, 'iEXF':6, 'iEsurf':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSsurf':13, 'iETg':14, 'iETsoil':15, 'iSsoil_pc':16, 'idSsoil':17, 'iinf':18, 'iHEADScorr':19, 'idgwt':20, 'iuzthick':21}
    # index_S = {'iEsoil':0, 'iTsoil':1,'iSsoil_pc':2, 'iRp':3, 'iRexf':4, 'idSsoil':5, 'iSsoil':6, 'iSAT':7, 'iMB_l':8}
    # TODO Ro sould not be averaged by ncell_MM, as well as EXF
    print '\nExporting output time series plots and ASCII files at catchment scale...'
    flxlbl       = ['RF', 'I', 'RFe', 'dSsurf', 'Esurf', 'Ro', 'dSsoil', 'EXF']
    flxlbl1      = ['Esoil', 'Tsoil']
    flxlbl2      = ['ETsoil']
    flxlbl2a     = ['Eg']    
    flxlbl3      = ['Tg']
    flxlbl3a     = ['ETg']
    flxlbl4      = ['Rp']
    flxlbl_tex   = [r'$RF$', r'$I$', r'$RFe$', r'$\Delta S_{surf}$', r'$E_{surf}$', r'$Ro$', r'$\Delta S_{soil}$', r'$EXF_g$']
    flxlbl1_tex  = [r'$E_{soil}$', r'$T_{soil}$']
    flxlbl2_tex  = [r'$ET_{soil}$']
    flxlbl3a_tex = [r'$ET_g$']
    flxlbl4_tex  = [r'$Rp$']
    flx_Cat_TS   = []
    flxmax_d     = []
    flxmin_d     = []
    TopSoilAverage = np.ma.masked_array(cMF.elev*1000.0, maskAllL).sum()*.001/sum(ncell_MM)
    for z, (i, i_tex) in enumerate(zip(flxlbl, flxlbl_tex)):
        i = 'i'+i
        array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
        flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
        flxmax_d.append(np.ma.max(flx_tmp))
        flxmin_d.append(np.ma.min(flx_tmp))
        array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
        flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
        if i == 'iRo':
            Ro_TOT = np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)
        del flx_tmp, array_tmp, array_tmp1
    for z, (i, i_tex) in enumerate(zip(flxlbl1, flxlbl1_tex)):
        flxlbl_tex.append(i_tex)
        i = 'i'+i
        flx_tmp1 = 0.0
        array_tmp2 = np.zeros((sum(cMF.perlen)), dtype = np.float)
        for l in range(_nslmax):
            array_tmp = h5_MM['MM_S'][:,:,:,l,index_S.get(i)]
            flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
            flxmax_d.append(np.ma.max(flx_tmp))
            flxmin_d.append(np.ma.min(flx_tmp))
            flx_tmp1 += flx_tmp.sum()
            array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
            array_tmp2 += np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)
        del flx_tmp, array_tmp, array_tmp1
        flx_Cat_TS.append(array_tmp2/sum(ncell_MM))
        del flx_tmp1, array_tmp2
    for z, (i, i_tex) in enumerate(zip(flxlbl2, flxlbl2_tex)):
        flxlbl_tex.append(i_tex)
        i = 'i'+i
        array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
        flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
        flxmax_d.append(np.ma.max(flx_tmp))
        flxmin_d.append(np.ma.min(flx_tmp))
        array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
        flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
        del flx_tmp, array_tmp, array_tmp1
    for z, i in enumerate(flxlbl2a):
        i = 'i'+i
        array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
        flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
        flxmax_d.append(np.ma.max(flx_tmp))
        flxmin_d.append(np.ma.min(flx_tmp))
        Eg_tmp = np.zeros((sum(cMF.perlen)), dtype = np.float)
        for L in range(cMF.nlay):
            mask_temp = mask_Lsup == (cMF.nlay-L)
            ncell_temp = float(np.sum(mask_temp))
            array_tmp1 = np.sum(array_tmp*mask_temp, axis = 1)
            array_tmp2 = np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/ncell_temp
            flx_Cat_TS.append(array_tmp2)
            Eg_tmp += array_tmp2
            flxlbl_tex.append(r'$E_g^{L%d}$'%(L+1))
        flx_Cat_TS.append(Eg_tmp)
        flxlbl_tex.append(r'$E_g$')
        del flx_tmp, array_tmp, array_tmp1, array_tmp2, mask_temp
    for z, i in enumerate(flxlbl3):        
        i = 'i'+i
        if cMF.wel_yn == 1:
            array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
            flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
            flxmax_d.append(np.ma.max(flx_tmp))
            Tg_min = np.ma.min(flx_tmp)
            flxmin_d.append(Tg_min)
            if abs(Tg_min) > 1E-6:
                print '\nTg negative (%.2f) observed at:' % Tg_min
                for row in range(cMF.nrow):
                    for t,col in enumerate(flx_tmp[:,row,:]):
                        try:
                            print 'row %d, col %d and day %d' % (row + 1, list(col).index(Tg_min) + 1, t + 1)
                            tTgmin = t
                            if plt_out_obs == 1:
                                obs['PzTgmin'] = {'x':999,'y':999, 'i': row, 'j': list(col).index(Tg_min), 'lay': 0, 'hi':999, 'h0':999, 'RC':999, 'STO':999, 'outpathname':os.path.join(MM_ws_out,'_0PzTgmin_ts.txt'), 'obs_h':[], 'obs_h_yn':0, 'obs_SM':[], 'obs_sm_yn':0, 'obs_Ro':[], 'obs_Ro_yn':0}
                                obs_list.append('Tgm')
                                obs4map[0].append('Tgm')
                                obs4map[1].append(row)
                                obs4map[2].append(list(col).index(Tg_min))
                                obs4map[3].append(0)
                                try:
                                    hmin.append(hmin[0])
                                except:
                                    hmin.append(999.9)
                        except:
                            pass
            Tg_tmp = np.zeros((sum(cMF.perlen)), dtype = np.float)
            for L in range(cMF.nlay):
                mask_temp = mask_Lsup == (cMF.nlay-L)
                ncell_temp = float(np.sum(mask_temp))
                array_tmp1 = np.sum(array_tmp*mask_temp, axis = 1)
                array_tmp2 = np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/ncell_temp
                flx_Cat_TS.append(array_tmp2)
                Tg_tmp += array_tmp2
                flxlbl_tex.append(r'$T_g^{L%d}$'%(L+1))
            flx_Cat_TS.append(Tg_tmp)
            flxlbl_tex.append(r'$T_g$')
            del flx_tmp, array_tmp, array_tmp1, array_tmp2, mask_temp
        else:
            flx_Cat_TS.append(0.0)
    for z, (i, i_tex) in enumerate(zip(flxlbl3a, flxlbl3a_tex)):
        flxlbl_tex.append(i_tex)
        i = 'i'+i
        array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
        flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
        flxmax_d.append(np.ma.max(flx_tmp))
        flxmin_d.append(np.ma.min(flx_tmp))
        array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
        flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
        del flx_tmp, array_tmp, array_tmp1
    for z, (i, i_tex) in enumerate(zip(flxlbl4, flxlbl4_tex)):
        flxlbl_tex.append(i_tex)
        i = 'i'+i
        array_tmp = h5_MM['finf_d']
        flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
        flxmax_d.append(conv_fact*np.ma.max(flx_tmp))
        flxmin_d.append(conv_fact*np.ma.min(flx_tmp))
        array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
        Rp = conv_fact*np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM)
        flx_Cat_TS.append(Rp)
        del flx_tmp, array_tmp
    for z, (i, i_tex) in enumerate(zip(['Ssurf', 'PE', 'PT', 'inf'], [r'$S_{surf}$', r'$PE$', r'$PT$', r'$inf$'])):
        flxlbl_tex.append(i_tex)
        i = 'i'+i
        array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
        array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
        flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
        del array_tmp, array_tmp1
    # ADD SM averaged
    flxlbl_tex.append(r'$\theta$')
    array_tmp = h5_MM['MM'][:,:,:,index.get('iSsoil_pc')]
    array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
    flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
    del array_tmp, array_tmp1
    h5_MM.close()
    flxmax_d = float(np.ceil(np.ma.max(flxmax_d)))
    flxmin_d = float(np.floor(np.ma.min(flxmin_d)))
    del flxlbl1, flxlbl2, flxlbl3, flxlbl3a
    if os.path.exists(cMF.h5_MF_fn):
        # compute UZF_STO and store GW_RCH
        try:
            h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
        except:
            cUTIL.ErrorExit('\nFATAL ERROR!\nInvalid MF HDF5 file. Run MARMITES and/or MODFLOW again.')        
        cbc_RCH = h5_MF['RCH_d']
        array_tmp2 = np.zeros((sum(cMF.perlen)), dtype = np.float)
        rch_tot = 0
        # GW_RCH
        for l in range(cMF.nlay):
            flxlbl_tex.append(r'$R^{L%d}$' % (l+1))
            array_tmp = cbc_RCH[:,l,:,:]
            array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
            rch_tmp = np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/ncell_MF[l]
            rch_tot += rch_tmp
            flx_Cat_TS.append(rch_tmp)
        flxlbl_tex.append(r'$R$')
        flx_Cat_TS.append(rch_tot)
        flxlbl_tex.append(r'$\Delta S_{u}$')
        flx_Cat_TS.append(rch_tot - Rp)
        del array_tmp, array_tmp1, rch_tmp, rch_tot, cbc_RCH, Rp
        for l in range(cMF.nlay):
            # ADD heads averaged
            flxlbl_tex.append(r'$h^{L%d}$' % (l+1))
            array_tmp = h_MF_m[:,l,:,:]
            array_tmp1 = np.sum(array_tmp, axis = 1)
            flx_Cat_TS.append(np.sum(array_tmp1, axis = 1)/ncell_MF[l])
            del array_tmp
            # ADD depth GWT
            flxlbl_tex.append(r'$d^{L%d}$' % (l+1))
            flx_Cat_TS.append(flx_Cat_TS[-1] - TopSoilAverage)
            del array_tmp1
        for l in range(cMF.nlay):
            # GW_STO
            cbc_STO = h5_MF['STO_d'][:,l,:,:]
            array_tmp1 = np.sum(np.ma.masked_values(cbc_STO, cMF.hnoflo, atol = 0.09), axis = 1)
            flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
            del cbc_STO, array_tmp1
            flxlbl_tex.append(r'$\Delta S_{g}^{L%d}$'%(l+1))
            # GW_FLF
            cbc_FLF = h5_MF['FLF_d'][:,l,:,:]
            array_tmp1 = np.sum(np.ma.masked_values(cbc_FLF, cMF.hnoflo, atol = 0.09), axis = 1)
            flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
            del cbc_FLF, array_tmp1
            flxlbl_tex.append('$FLF^{L%d}$'%(l+1))
            # EXF
            cbc_EXF = h5_MF['EXF_d'][:,l,:,:]
            array_tmp1 = np.sum(np.ma.masked_values(cbc_EXF, cMF.hnoflo, atol = 0.09), axis = 1)
            flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
            del cbc_EXF, array_tmp1            
            flxlbl_tex.append('$EXF^{L%d}$'%(l+1))
            # WEL
            if cMF.wel_yn == 1:
                cbc_WEL = h5_MF['WEL_d'][:,l,:,:]
                array_tmp1 = np.sum(np.ma.masked_values(cbc_WEL, cMF.hnoflo, atol = 0.09), axis = 1)
                flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
                flxlbl_tex.append('$WEL^{L%d}$'%(l+1))
                del cbc_WEL               
           # DRN
            if cMF.drn_yn == 1:
                cbc_DRN = h5_MF['DRN_d'][:,l,:,:]
                if cMF.drncells[l]>0:
                    array_tmp1 = np.sum(np.ma.masked_values(cbc_DRN, cMF.hnoflo, atol = 0.09), axis = 1)
                    flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
                    del array_tmp1
                    flxlbl_tex.append('$DRN^{L%d}$'%(l+1))
                del cbc_DRN
            # GHB
            if cMF.ghb_yn == 1:
                cbc_GHB = h5_MF['GHB_d'][:,l,:,:]
                if cMF.ghbcells[l] > 0:
                    array_tmp1 = np.sum(np.ma.masked_values(cbc_GHB, cMF.hnoflo, atol = 0.09), axis = 1)
                    flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
                    flxlbl_tex.append('$GHB^{L%d}$'%(l+1))
                    del array_tmp1
                del cbc_GHB                
        h5_MF.close()
        plt_exportCATCH_fn = os.path.join(MM_ws_out, '_0CATCHMENT_ts.png')
        plt_titleCATCH = 'Time serie of fluxes averaged over the whole catchment'
        rmseHEADS_tmp, rmseSM_tmp, rsrHEADS_tmp, rsrSM_tmp, nseHEADS_tmp, nseSM_tmp, rHEADS_tmp, rSM_tmp = MMplot.plotTIMESERIES_CATCH(cMF, flx_Cat_TS, flxlbl_tex, Ro_TOT, plt_exportCATCH_fn, plt_titleCATCH, hmax = hmaxMF, hmin = hminMF, iniMonthHydroYear = iniMonthHydroYear, date_ini = StartDate, date_end = EndDate, obs_catch = obs_catch, obs_catch_list = obs_catch_list, TopSoilAverage = TopSoilAverage, MF = 1)
    if rmseHEADS_tmp <> None:
        rmseHEADS.append(rmseHEADS_tmp)
        rsrHEADS.append(rsrHEADS_tmp)
        nseHEADS.append(nseHEADS_tmp)
        rHEADS.append(rHEADS_tmp)
        obslstHEADS.append('catch.')
    if rmseSM_tmp <> None:
        rmseSM.append(rmseSM_tmp)
        rsrSM.append(rsrSM_tmp)
        nseSM.append(nseSM_tmp)
        rSM.append(rSM_tmp)
        obslstSM.append('catch.')
    del rmseHEADS_tmp, rmseSM_tmp, rsrHEADS_tmp, rsrSM_tmp, nseHEADS_tmp, nseSM_tmp, rHEADS_tmp, rSM_tmp
    # export average time serie of fluxes in txt
    plt_exportCATCH_txt_fn = os.path.join(MM_ws_out, '_0CATCHMENT_ts.txt')
    plt_exportCATCH_txt = open(plt_exportCATCH_txt_fn, 'w')
    flxlbl_CATCH_str = 'Date'
    for e in flxlbl_tex:
        flxlbl_CATCH_str += ',' + e
    flxlbl_CATCH_str += ',%s' % '$hobs$'
    flxlbl_CATCH_str += ',%s' % '$\\thetaobs$'
    plt_exportCATCH_txt.write(flxlbl_CATCH_str)
    plt_exportCATCH_txt.write('\n')
    for t in range(len(cMF.inputDate)):
        flx_Cat_TS_str = str(flx_Cat_TS[0][t])
        for e in (flx_Cat_TS[1:]):
            flx_Cat_TS_str += ',%s' % str(e[t])
        if obs_catch_list[0] == 1:
            flx_Cat_TS_str += ',%s' % (obs_catch.get('catch')['obs_h'][0][t])
        else:
            flx_Cat_TS_str += ',%s' % cMF.hnoflo
        if obs_catch_list[1] == 1:
            flx_Cat_TS_str += ',%s' % (obs_catch.get('catch')['obs_SM'][0][t])
        else:
            flx_Cat_TS_str += ',%s' % cMF.hnoflo
        out_line = '%s,%s' % (mpl.dates.num2date(cMF.inputDate[t]).isoformat()[:10], flx_Cat_TS_str)
        for l in out_line:
            plt_exportCATCH_txt.write(l)
        plt_exportCATCH_txt.write('\n')
    plt_exportCATCH_txt.close()
    print '-------\nExporting water balance at the catchment scale...'
    try:
        MMplot.plotWBsankey(path = MM_ws_out, fn = '_0CATCHMENT_ts.txt', index = HYindex, year_lst = year_lst, cMF = cMF, ncell_MM = ncell_MM)
    except:
        print "\nError in plotting the catchment water balance!"
    del flx_Cat_TS, flx_Cat_TS_str, out_line, plt_exportCATCH_fn, plt_exportCATCH_txt_fn, plt_titleCATCH, Ro_TOT

# #################################################
# EXPORT AT OBSERVATION POINTS
# exporting MM time series results to ASCII files and plots at observations cells
# #################################################
if plt_out_obs == 1 and os.path.exists(h5_MM_fn) and os.path.exists(cMF.h5_MF_fn):
    print '\nExporting output time series plots and ASCII files at observation points...'
    h5_MM = h5py.File(h5_MM_fn, 'r')
    h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
    clr_lst = ['darkgreen', 'firebrick', 'darkmagenta', 'goldenrod', 'green', 'tomato', 'magenta', 'yellow']
    x = 0
    for o_ref in obs_list:
        for o in obs.keys():
            if o == o_ref:
                i = obs.get(o)['i']
                j = obs.get(o)['j']
                l_obs = obs.get(o)['lay']
                obs_h = obs.get(o)['obs_h']
                obs_SM = obs.get(o)['obs_SM']
                obs_Ro = obs.get(o)['obs_Ro']
                outFileExport = open(obs.get(o)['outpathname'], 'w')
                outFileExport.write(header)
                SOILzone_tmp = gridSOIL[i,j]-1
                nsl = _nsl[SOILzone_tmp]
                soilnam = _nam_soil[SOILzone_tmp]
                slprop = _slprop[SOILzone_tmp]
                # thickness of soil layers
                Tl = list(gridSOILthick[i,j]*slprop)
                for ii, ll in enumerate(Tl):
                    Tl[ii] = float('%.3f' % ll)
                MM = h5_MM['MM'][:,i,j,:]
                MM_S = h5_MM['MM_S'][:,i,j,:,:]
                # SATFLOW
                cbc_RCH = h5_MF['RCH_d']
                h_satflow = MM_SATFLOW.run(cbc_RCH[:,l_obs,i,j], float(obs.get(o)['hi']),float(obs.get(o)['h0']),float(obs.get(o)['RC']),float(obs.get(o)['STO']))
                # export ASCII file at piezometers location
                #TODO extract heads at piezo location and not center of cell
                if obs_h != []:
                    obs_h_tmp = obs_h[0,:]
                else:
                    obs_h_tmp = []
                if obs_Ro != []:
                    obs_Ro_tmp = obs_Ro[0,:]
                else:
                    obs_Ro_tmp = []
                if cMF.wel_yn == 1:
                    cbc_WEL = np.sum(np.ma.masked_values(h5_MF['WEL_d'][:,:,i,j], cMF.hnoflo, atol = 0.09), axis = 1)
                else:
                    cbc_WEL = 0
                if irr_yn == 0:
                    index_veg = np.ones((NVEG, sum(cMF.perlen)), dtype = float)
                    for v in range(NVEG):
                        index_veg[v] = cMF.LAI_veg_d[gridMETEO[i,j]-1,v,:]*gridVEGarea[v,i,j]
                    index_veg = (np.sum(index_veg, axis = 0) > 0.0)*(-1)
                else:
                    IRRfield = gridIRR[i,j]
                    if IRRfield == 0:
                        index_veg = np.ones((NVEG, sum(cMF.perlen)), dtype = float)
                        for v in range(NVEG):
                            index_veg[v] = cMF.LAI_veg_d[gridMETEO[i,j]-1,v,:]*gridVEGarea[v,i,j]
                        index_veg = (np.sum(index_veg, axis = 0) > 0.0)*(-1)
                    else:
                        index_veg = cMF.crop_irr_d[gridMETEO[i,j]-1, IRRfield-1,:]
                # Export time series results at observations points as ASCII file
                cMF.cPROCESS.ExportResultsMM(i, j, cMF.inputDate, SP_d, _nslmax, MM, index, MM_S, index_S, cbc_RCH[:,l_obs,i,j], cbc_WEL, h_satflow, h_MF_m[:,l_obs,i,j], obs_h_tmp, obs_SM, obs_Ro_tmp, index_veg, outFileExport, o)
                del cbc_WEL
                outFileExport.close()
                # plot time series results as plot
                plt_suptitle = 'Time serie of fluxes at observation point %s' % o
                plt_title = 'i = %d, j = %d, l = %d, x = %.2f, y = %.2f, elev. = %.2f, %s\nSm = %s, Sfc = %s, Sr = %s, Ks = %s, thick. = %.3f %s' % (i+1, j+1, l_obs+1, obs.get(o)['x'], obs.get(o)['y'], cMF.elev[i,j], soilnam, _Sm[SOILzone_tmp], _Sfc[SOILzone_tmp], _Sr[SOILzone_tmp], _Ks[SOILzone_tmp], gridSOILthick[i,j], Tl)
                # index = {'iRF':0, 'iPT':1, 'iPE':2, 'iRFe':3, 'iSsurf':4, 'iRo':5, 'iEXF':6, 'iEsurf':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSsurf':13, 'iETg':14, 'iETsoil':15, 'iSsoil_pc':16, 'idSsoil':17, 'iinf':18, 'iHEADScorr':19, 'idgwt':20, 'iuzthick':21}
                # index_S = {'iEsoil':0, 'iTsoil':1,'iSsoil_pc':2, 'iRp':3, 'iRexf':4, 'idSsoil':5, 'iSsoil':6, 'iSAT':7, 'iMB_l':8}
                plt_export_fn = os.path.join(MM_ws_out, '_0'+ o + '_ts.png')
                # def plotTIMESERIES(DateInput, P, PT, PE, Pe, dPOND, POND, Ro, Eu, Tu, Eg, Tg, S, dS, Spc, Rp, EXF, ETg, Es, MB, MB_l, dgwt, SAT, R, h_MF, h_MF_corr, h_SF, hobs, Sobs, Sm, Sr, hnoflo, plt_export_fn, plt_title, colors_nsl, hmax, hmin):
                MMplot.plotTIMESERIES(
                cMF,
                MM[:,index.get('iRF')],
                MM[:,index.get('iPT')],
                MM[:,index.get('iPE')],
                MM[:,index.get('iRFe')],
                MM[:,index.get('idSsurf')],
                MM[:,index.get('iSsurf')],
                MM[:,index.get('iRo')],
                MM_S[:,0:nsl,index_S.get('iEsoil')],
                MM_S[:,0:nsl,index_S.get('iTsoil')],
                MM[:,index.get('iEg')],
                MM[:,index.get('iTg')],
                MM_S[:,0:nsl,index_S.get('iSsoil')],
                MM_S[:,0:nsl,index_S.get('idSsoil')],
                MM_S[:,0:nsl,index_S.get('iSsoil_pc')],
                MM_S[:,0:nsl,index_S.get('iRp')],
                MM[:,index.get('iEXF')],
                MM[:,index.get('iETg')],
                MM[:,index.get('iEsurf')],
                MM[:,index.get('iMB')],
                MM_S[:,0:nsl,index_S.get('iMB_l')],
                MM[:,index.get('idgwt')],
                MM[:,index.get('iuzthick')],
                MM_S[:,0:nsl,index_S.get('iSAT')],
                cbc_RCH[:,l_obs,i,j],
                h_MF_m[:,:,i,j], MM[:,index.get('iHEADScorr')], h_satflow, obs_h_tmp, obs_SM, obs_Ro_tmp,
                _Sm[gridSOIL[i,j]-1],
                _Sr[gridSOIL[i,j]-1],
                cMF.hnoflo,
                plt_export_fn, plt_suptitle, plt_title,
                clr_lst,
                max(hmax), #hmax[x] + hdiff/2
                min(hmin), #hmin[x] - hdiff/2
                o,
                cMF.elev[i,j],
                cMF.nlay,
                l_obs,
                nsl,
                iniMonthHydroYear, date_ini = StartDate, date_end = EndDate
                )
                x += 1
    del i, j, l_obs, SOILzone_tmp, outFileExport, nsl, soilnam, slprop, Tl, plt_export_fn
    h5_MM.close()
    h5_MF.close()

# #################################################
# PLOT SPATIAL MF and MM OUTPUT
# #################################################
if plt_out == 1:
    print '\nExporting output maps...'
    day_lst = []
    Date_lst = []
    JD_lst = []
    day = HYindex[1]
    while day < sum(cMF.perlen):  #len(h_MF_m):
        day_lst.append(day)
        Date_lst.append(cMF.inputDate[day])
        JD_lst.append(cMF.JD[day])
        day += plt_freq
    if animation < 1:
        if tTgmin < 0 and tRCHmax > 0:
            lst = [sum(cMF.perlen) - 1, tRCHmax] #[len(h_MF_m)-1, tRCHmax]
        else:
            lst = [sum(cMF.perlen) - 1, tRCHmax, tTgmin] # [len(h_MF_m)-1, tRCHmax, tTgmin]
        for e in lst:
            day_lst.append(e)
            Date_lst.append(cMF.inputDate[e])
            JD_lst.append(cMF.JD[e])
    del day

    # ##############################
    # #### MODFLOW OUTPUT ##########
    # ##############################
    if os.path.exists(cMF.h5_MF_fn):
        h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
        h5_MM = h5py.File(h5_MM_fn, 'r')
        cbc_RCH = h5_MF['RCH_d']
        cbc_EXF = h5_MF['EXF_d']
        cbc_WEL = h5_MF['WEL_d']
        if cMF.drn_yn == 1:
            cbc_DRN = h5_MF['DRN_d']
        if cMF.ghb_yn == 1:
            cbc_GHB = h5_MF['GHB_d']
        # ############################################
        # plot at specified day
        # ############################################

        # plot heads [m]
        V = np.zeros((len(day_lst), cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        mask_tmp = np.zeros((cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        maskAllL_tmp  = np.zeros((cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        Vmax = np.zeros((len(day_lst)), dtype = np.float)
        Vmin = np.zeros((len(day_lst)), dtype = np.float)
        Vmax1 = np.zeros((len(day_lst)), dtype = np.float)
        Vmin1 = np.zeros((len(day_lst)), dtype = np.float)
        for i, t in enumerate(day_lst):
            for L in range(cMF.nlay):
                V[i,L,:,:] = h_MF_m[t,L,:,:]
                mask_tmp[L,:,:] = mask[L]
                maskAllL_tmp[L,:,:] = maskAllL
            Vmax[i] = hmaxMF
            Vmin[i] = hminMF
        MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation - $h$ $(m)$', msg = 'DRY', plt_title = 'OUT_MF_HEADS', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsMF, Vmax = Vmax, Vmin = Vmin, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo, animation = animation)

        # plot GWTD [m]
        for i in range(len(day_lst)):
           for L in range(cMF.nlay):
            V[i,L,:,:] = cMF.elev - V[i,L,:,:]
            Vmin[i] = GWTDmin
            Vmax[i] = GWTDmax
        MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'depth to groundwater table - $d$ $(m)$', msg = 'DRY', plt_title = 'OUT_MF_GWTD', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsMF, Vmax = Vmax, Vmin = Vmin, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo, animation = animation)

        # plot heads corrigidas [m]
        headscorr_m = np.zeros((len(day_lst), cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        for i, t in enumerate(day_lst):
            headscorr_m[i,0,:,:] = np.ma.masked_values(np.ma.masked_values(h5_MM['MM'][t,:,:,19], cMF.hnoflo, atol = 0.09), cMF.hdry, atol = 1E+25)
            Vmin[i] = hcorrmin
            Vmax[i] = hcorrmax
        MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = headscorr_m,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation - $h$ $(m)$', msg = 'DRY', plt_title = 'OUT_MF_HEADScorr', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsMF, Vmax = Vmax, Vmin = Vmin, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo, animation = animation)

        # plot GWTD correct [m]
        GWTDcorr = np.zeros((len(day_lst), cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        for i, t in enumerate(day_lst):
            GWTDcorr[i,0,:,:] = cMF.elev-headscorr_m[i,0,:,:]
            Vmin[i] = GWTDcorrmin
            Vmax[i] = GWTDcorrmax
        MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = GWTDcorr,  cmap = plt.cm.Blues, CBlabel = 'depth to groundwater table - $d$ $(m)$', msg = 'DRY', plt_title = 'OUT_MF_GWTDcorr', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsMF, Vmax = Vmax, Vmin = Vmin, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo, animation = animation)

        # plot GW GROSS RCH [mm]
        R = np.zeros((len(day_lst), cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        for i, t in enumerate(day_lst):
            for L in range(cMF.nlay):
                R[i,L,:,:] = np.ma.masked_array(cbc_RCH[t,L,:,:], mask[L])
            Vmin[i] = RCHmin
            Vmax[i] = RCHmax
            Vmin1[i] = np.ma.min(R[i,:,:,:])
            Vmax1[i] = np.ma.max(R[i,:,:,:])
        MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = R,  cmap = plt.cm.Blues, CBlabel = 'groundwater gross recharge - $R$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_MF_R', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = Vmin, contours = ctrsMF, Vmax = Vmax, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = mask_tmp, animation = animation)
        MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = R,  cmap = plt.cm.Blues, CBlabel = 'groundwater gross recharge - $R$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_MF_R1', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = Vmin1, contours = ctrsMF, Vmax = Vmax1, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = mask_tmp, animation = animation)
        
        # plot GW EFFECTIVE RCH [mm]
        Re = np.zeros((len(day_lst), cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        for i, t in enumerate(day_lst):
            for L in range(cMF.nlay):
                Re[i,L,:,:] = R[i,L,:,:] + np.ma.masked_array(cbc_EXF[t,L,:,:], mask[L])
            Vmin[i] = np.ma.min(Re)
            Vmax[i] = np.ma.max(Re)
            Vmin1[i] = np.ma.min(Re[i,:,:,:])
            Vmax1[i] = np.ma.max(Re[i,:,:,:])
        MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = Re,  cmap = plt.cm.Blues, CBlabel = 'groundwater effective recharge - $Re$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_MF_Re', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = Vmin, contours = ctrsMF, Vmax = Vmax, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = mask_tmp, animation = animation)
        MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = Re,  cmap = plt.cm.Blues, CBlabel = 'groundwater effective recharge - $Re$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_MF_Re1', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = Vmin1, contours = ctrsMF, Vmax = Vmax1, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = mask_tmp, animation = animation)

        # plot GW EFFECTIVE RCH [mm]
        Rn = np.zeros((len(day_lst), cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        for i, t in enumerate(day_lst):
            for L in range(cMF.nlay):
                Rn[i,L,:,:] = Re[i,L,:,:] - np.ma.masked_array(cbc_WEL[t,L,:,:], mask[L])
            Vmin[i] = np.ma.min(Rn)
            Vmax[i] = np.ma.max(Rn)
            Vmin1[i] = np.ma.min(Rn[i,:,:,:])
            Vmax1[i] = np.ma.max(Rn[i,:,:,:])
        MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = Rn,  cmap = plt.cm.Blues, CBlabel = 'groundwater net recharge - $Rn$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_MF_Rn', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = Vmin, contours = ctrsMF, Vmax = Vmax, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = mask_tmp, animation = animation)
        MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = Rn,  cmap = plt.cm.Blues, CBlabel = 'groundwater net recharge - $Rn$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_MF_Rn1', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = Vmin1, contours = ctrsMF, Vmax = Vmax1, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = mask_tmp, animation = animation)  
    
        del R, Rn, Re        

        # plot GW drainage [mm]
        if cMF.drn_yn == 1:
            V = np.zeros((len(day_lst), cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
            for i, t in enumerate(day_lst):
                for L in range(cMF.nlay):
                    V[i,L,:,:] = np.ma.masked_array(cbc_DRN[t,L,:,:], mask[L])*(-1.0)
                Vmin[i] = DRNmin
                Vmax[i] = DRNmax
                Vmin1[i] = np.ma.min(V[i,:,:,:])
                Vmax1[i] = np.ma.max(V[i,:,:,:])
            MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage - $DRN$ $(mm/day)$', msg = '- no drainage', plt_title = 'OUT_MF_DRN', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = Vmin, contours = ctrsMF, Vmax = Vmax, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo, animation = animation)
            MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage - $DRN$ $(mm/day)$', msg = '- no drainage', plt_title = 'OUT_MF_DRN1', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = Vmin1, contours = ctrsMF, Vmax = Vmax1, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo, animation = animation)

        # plot GHB [mm]
        if cMF.ghb_yn == 1:
            V = np.zeros((len(day_lst), cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
            for i, t in enumerate(day_lst):
                for L in range(cMF.nlay):
                    V[i,L,:,:] = np.ma.masked_array(cbc_GHB[t,L,:,:], mask[L])*(-1.0)
                Vmin[i] = GHBmin
                Vmax[i] = GHBmax
                Vmin1[i] = np.ma.min(V[i,:,:,:])
                Vmax1[i] = np.ma.max(V[i,:,:,:])
            MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'general head bdry - $GHB$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_MF_GHB', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = Vmin, contours = ctrsMF, Vmax = Vmax, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = mask_tmp, animation = animation)
            MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'general head bdry - $GHB$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_MF_GHB1', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = Vmin1, contours = ctrsMF, Vmax = Vmax1, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo,mask = mask_tmp, animation = animation)

        del V
        del Vmin1, Vmax1, Vmin, Vmax

    # ############################################
    # plot average of all hydrologic years
    # ############################################

    # plot heads average [m]
    V = np.zeros((1, cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
    for L in range(cMF.nlay):
        V[0,L,:,:] = np.ma.masked_array(np.sum(h_MF_m[HYindex[1]:HYindex[-1],L,:,:], axis = 0)/(HYindex[-1]-HYindex[1]+1), mask[L])
    MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation - $h$ $(m)$', msg = 'DRY', plt_title = 'OUT_average_MF_HEADS', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsMF, Vmax = [hmaxMF], Vmin = [hminMF], ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)

    # plot GWTD average [m]
    for L in range(cMF.nlay):
        V[0,L,:,:] = cMF.elev-V[0,L,:,:]
    MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'depth to groundwater table - $d$ $(m)$', msg = 'DRY', plt_title = 'OUT_average_MF_GWTD', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsMF, Vmax = [GWTDmax], Vmin = [GWTDmin], ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)

    # plot heads corrigidas average [m]
    headscorr_m = np.zeros((1, cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
    headscorr_m[0,0,:,:] = np.ma.masked_values(np.ma.masked_values(np.sum(h5_MM['MM'][HYindex[1]:HYindex[-1],:,:,19], axis = 0)/(HYindex[-1]-HYindex[1]+1), cMF.hnoflo, atol = 0.09), cMF.hdry, atol = 1E+25)
    MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = headscorr_m,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation - $h$ $(m)$', msg = 'DRY', plt_title = 'OUT_average_MF_HEADScorr', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsMF, Vmax = [hcorrmax], Vmin = [hcorrmin], ntick = ntick, points = obs4map, mask = maskAllL_tmp, hnoflo = cMF.hnoflo)

    # plot GWTD correct average [m]
    GWTDcorr = np.zeros((1, cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
    GWTDcorr[0,0,:,:] = cMF.elev-headscorr_m[0,0,:,:]
    MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = GWTDcorr, cmap = plt.cm.Blues, CBlabel = 'depth to groundwater table - $d$ $(m)$', msg = 'DRY', plt_title = 'OUT_average_MF_GWTDcorr', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsMF, Vmax = [GWTDcorrmax], Vmin = [GWTDcorrmin], ntick = ntick, points = obs4map, mask = maskAllL_tmp, hnoflo = cMF.hnoflo)

    # plot GW GROSS RCH average [mm]
    R = np.zeros((1, cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
    for L in range(cMF.nlay):
        R[0,L,:,:] = np.ma.masked_array(np.sum(cbc_RCH[HYindex[1]:HYindex[-1],L,:,:], axis = 0)/(HYindex[-1]-HYindex[1]+1), mask[L])
    MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = R,  cmap = plt.cm.Blues, CBlabel = 'groundwater gross recharge - $R$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_average_MF_R', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = [RCHmax], Vmin = [RCHmin], contours = ctrsMF, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)
    Vmin_tmp1 = np.min(R)
    Vmax_tmp1 = np.max(R)
    MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = R,  cmap = plt.cm.Blues, CBlabel = 'groundwater gross recharge - $R$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_average_MF_R1', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = [Vmax_tmp1], Vmin = [Vmin_tmp1], contours = ctrsMF, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)
    
    # plot GW EFFECTIVE RCH average [mm]
    Re = np.zeros((1, cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
    for L in range(cMF.nlay):
        Re[0,L,:,:] = R[0,L,:,:] + np.ma.masked_array(np.sum(cbc_EXF[HYindex[1]:HYindex[-1],L,:,:], axis = 0)/(HYindex[-1]-HYindex[1]+1), mask[L])
    MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = Re,  cmap = plt.cm.Blues, CBlabel = 'groundwater effective recharge - $Re$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_average_MF_Re', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = [np.ma.max(Re)], Vmin = [np.ma.min(Re)], contours = ctrsMF, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)
    Vmin_tmp1 = np.min(Re)
    Vmax_tmp1 = np.max(Re)
    MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = Re,  cmap = plt.cm.Blues, CBlabel = 'groundwater effective recharge - $Re$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_average_MF_Re1', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = [Vmax_tmp1], Vmin = [Vmin_tmp1], contours = ctrsMF, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)
    
    # plot GW NET RCH average [mm]
    Rn = np.zeros((1, cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
    for L in range(cMF.nlay):
        Rn[0,L,:,:] = Re[0,L,:,:] - np.ma.masked_array(np.sum(cbc_WEL[HYindex[1]:HYindex[-1],L,:,:], axis = 0)/(HYindex[-1]-HYindex[1]+1), mask[L])
    MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = Rn,  cmap = plt.cm.Blues, CBlabel = 'groundwater net recharge - $Rn$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_average_MF_Rn', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = [np.ma.min(Rn)], Vmin = [np.ma.min(Rn)], contours = ctrsMF, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)
    Vmin_tmp1 = np.min(Rn)
    Vmax_tmp1 = np.max(Rn)
    MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = Rn,  cmap = plt.cm.Blues, CBlabel = 'groundwater net recharge - $Rn$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_average_MF_Rn1', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = [Vmax_tmp1], Vmin = [Vmin_tmp1], contours = ctrsMF, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)    
    
    del cbc_RCH, cbc_EXF, cbc_WEL, R, Rn, Re

    # plot GW drainage average [mm]
    if cMF.drn_yn == 1:
        V = np.zeros((1, cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        for L in range(cMF.nlay):
            V[0,L,:,:] = np.ma.masked_array(np.sum(cbc_DRN[HYindex[1]:HYindex[-1],L,:,:], axis = 0)/(HYindex[-1]-HYindex[1]+1)*(-1.0), mask[L])
        MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage - $DRN$ $(mm/day)$', msg = '- no drainage', plt_title = 'OUT_average_MF_DRN', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = [DRNmax], Vmin = [DRNmin], contours = ctrsMF, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)
        del cbc_DRN, DRNmax, DRNmin
        Vmin_tmp1 = np.min(V)
        Vmax_tmp1 = np.max(V)
        MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage - $DRN$ $(mm/day)$', msg = '- no drainage', plt_title = 'OUT_average_MF_DRN1', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = [Vmax_tmp1], Vmin = [Vmin_tmp1], contours = ctrsMF, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)

    # plot GHB average [mm]
    if cMF.ghb_yn == 1:
        V = np.zeros((1, cMF.nlay, cMF.nrow, cMF.ncol), dtype = np.float)
        for L in range(cMF.nlay):
            V[0,L,:,:] = np.ma.masked_array(np.sum(cbc_GHB[HYindex[1]:HYindex[-1],L,:,:], axis = 0)/(HYindex[-1]-HYindex[1]+1)*(-1.0), mask[L])
        MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'general head bdry - $GHB$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_average_MF_GHB', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = [GHBmax], Vmin = [GHBmin], contours = ctrsMF, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)
        del cbc_GHB, GHBmax, GHBmin
        Vmin_tmp1 = np.min(V)
        Vmax_tmp1 = np.max(V)
        MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'general head bdry - $GHB$ $(mm/day)$', msg = '- no flux', plt_title = 'OUT_average_MF_GHB1', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = [Vmax_tmp1], Vmin = [Vmin_tmp1], contours = ctrsMF, ntick = ntick, points = obs4map, mask = mask_tmp, hnoflo = cMF.hnoflo)

    h5_MF.close()
    h5_MM.close()
    del V, Vmax_tmp1, Vmin_tmp1
##
    # ##############################
    # #### MARMITES OUTPUT #########
    # ##############################
    if os.path.exists(h5_MM_fn):
        h5_MM = h5py.File(h5_MM_fn, 'r')
        flxlbl = ['RF', 'RFe', 'I', 'EXF', 'dSsurf', 'Ro', 'Esurf', 'Eg', 'Tg', 'ETg', 'ETsoil', 'dSsoil']
        flxlbl_tex = [r'$RF$', r'$RFe$', r'$I$', r'$EXF_g$', r'$\Delta S_{surf}$', r'$Ro$', r'$E_{surf}$', r'$E_g$', r'$T_g$', r'$ET_g$', r'$ET_{soil}$', r'$\Delta S_{soil}$']
        for z, (i, i_lbl) in enumerate(zip(flxlbl, flxlbl_tex)):
            # ############################################
            # plot average of all hydrologic years
            # ############################################
            i1 = 'i'+i
            MM = h5_MM['MM'][:,:,:,index.get(i1)]
            V = np.zeros((1, 1, cMF.nrow, cMF.ncol), dtype = np.float)
            V[0,0,:,:] = np.sum(np.ma.masked_values(MM[HYindex[1]:HYindex[-1],:,:], cMF.hnoflo, atol = 0.09), axis = 0)/(HYindex[-1]-HYindex[1]+1)
            V[0,0,:,:] = np.ma.masked_values(V[0,0,:,:], cMF.hnoflo, atol = 0.09)
            Vmax = [np.ma.max(V)]
            Vmin = [np.ma.min(V)]
            MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' $(mm/day)$'), msg = 'no flux', plt_title = ('OUT_average_MM_' + i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax, Vmin = Vmin, contours = ctrsMM, ntick = ntick, points = obs4map, mask = maskAllL_tmp, hnoflo = cMF.hnoflo)
            del V
            # ############################################
            # plot for selected time step
            # ############################################
            V = np.zeros((len(day_lst), 1, cMF.nrow, cMF.ncol), dtype = np.float)
            Vmax = np.zeros((len(day_lst)), dtype = np.float)
            Vmin = np.zeros((len(day_lst)), dtype = np.float)
            Vmax1 = np.zeros((len(day_lst)), dtype = np.float)
            Vmin1 = np.zeros((len(day_lst)), dtype = np.float)
            for ii, t in enumerate(day_lst):
                V[ii,0,:,:] = np.ma.masked_values(MM[t,:,:], cMF.hnoflo, atol = 0.09)
                V[ii,0,:,:] = np.ma.masked_values(V[ii,0,:,:], cMF.hnoflo, atol = 0.09)
                Vmin1[ii] = np.ma.min(np.ma.masked_values(V[ii,0,:,:], cMF.hnoflo, atol = 0.09))
                Vmax1[ii] = np.ma.max(np.ma.masked_values(V[ii,0,:,:], cMF.hnoflo, atol = 0.09))
            for ii, t in enumerate(day_lst):
                Vmax[ii] = np.ma.max(Vmax1)
                Vmin[ii] = np.ma.min(Vmin1)
            MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' $(mm/day)$'), msg = 'no flux', plt_title = ('OUT_MM_'+i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax, Vmin = Vmin, contours = ctrsMM, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = maskAllL_tmp, animation = animation)
            MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' $(mm/day)$'), msg = 'no flux', plt_title = ('OUT_MM_%s1'%i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax1, Vmin = Vmin1, contours = ctrsMM, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = maskAllL_tmp, animation = animation)
            del V, MM, Vmax, Vmin

        flxlbl = ['Esoil', 'Tsoil']
        flxlbl_tex = [r'$E_{soil}$', r'$T_{soil}$']
        for z, (i, i_lbl) in enumerate(zip(flxlbl, flxlbl_tex)):
            # ############################################
            # plot average of all hydrologic years
            # ############################################
            i1 = 'i'+i
            V = np.zeros((1, 1, cMF.nrow, cMF.ncol), dtype = np.float)
            for l in range(_nslmax):
                MM = h5_MM['MM_S'][:,:,:,l,index_S.get(i1)]
                V[0,0,:,:] += np.sum(np.ma.masked_values(MM[HYindex[1]:HYindex[-1],:,:], cMF.hnoflo, atol = 0.09), axis = 0)/(HYindex[-1]-HYindex[1]+1)
                V[0,0,:,:] = np.ma.masked_values(V[0,0,:,:], cMF.hnoflo, atol = 0.09)
            Vmax = [np.ma.max(V)]
            Vmin = [np.ma.min(V)]
            MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' $(mm/day)$'), msg = 'no flux', plt_title = ('OUT_average_MM_%s' % i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax, Vmin = Vmin, contours = ctrsMM, ntick = ntick, points = obs4map, mask = maskAllL_tmp, hnoflo = cMF.hnoflo)
            del V

            # ############################################
            # plot for selected time step
            # ############################################
            Vmax = np.zeros((len(day_lst)), dtype = np.float)
            Vmin = np.zeros((len(day_lst)), dtype = np.float)
            Vmax1 = np.zeros((len(day_lst)), dtype = np.float)
            Vmin1 = np.zeros((len(day_lst)), dtype = np.float)
            V = np.zeros((len(day_lst), 1, cMF.nrow, cMF.ncol), dtype = np.float)
            for ii, t in enumerate(day_lst):
                for l in range(_nslmax):
                    V[ii,0,:,:] += h5_MM['MM_S'][t,:,:,l,index_S.get(i1)]
                Vmin1[ii] = np.ma.min(np.ma.masked_values(np.ma.masked_array(V[ii,0,:,:], maskAllL), cMF.hnoflo, atol = 0.09))
                Vmax1[ii] = np.ma.max(np.ma.masked_values(np.ma.masked_array(V[ii,0,:,:], maskAllL), cMF.hnoflo, atol = 0.09))
            for ii, t in enumerate(day_lst):
                Vmax[ii] = np.ma.max(Vmax1)
                Vmin[ii] = np.ma.min(Vmin1)
            MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' $(mm/day)$'), msg = 'no flux', plt_title = ('OUT_MM_'+i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax, Vmin = Vmin, contours = ctrsMM, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = maskAllL_tmp, animation = animation)
            MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' $(mm/day)$'), msg = 'no flux', plt_title = ('OUT_MM_%s1'%i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax1, Vmin = Vmin1, contours = ctrsMM, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = maskAllL_tmp, animation = animation)

        # ##############################
        # ####   PERCOLATION   #########
        # ##############################
        i = 'Rp'
        i_lbl = r'$Rp_{soil}$'
        # ############################################
        # plot average of all hydrologic years
        # ############################################
        V = np.zeros((1, 1, cMF.nrow, cMF.ncol), dtype = np.float)
        MM = conv_fact*h5_MM['finf_d'][:,:,:]
        V[0,0,:,:] = np.sum(np.ma.masked_values(MM[HYindex[1]:HYindex[-1],:,:], cMF.hnoflo, atol = 0.09), axis = 0)/(HYindex[-1]-HYindex[1]+1)   
        Vmax = [np.ma.max(V)]
        Vmin = [np.ma.min(V)]
        MMplot.plotLAYER(timesteps = ['NA'], Date = ['NA'], JD = ['NA'], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' $(mm/day)$'), msg = 'no flux', plt_title = ('OUT_average_MM_%s'% i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax, Vmin = Vmin, contours = ctrsMM, ntick = ntick, points = obs4map, mask = maskAllL_tmp, hnoflo = cMF.hnoflo)

        # ############################################
        # plot for selected time step
        # ############################################
        Vmax = np.zeros((len(day_lst)), dtype = np.float)
        Vmin = np.zeros((len(day_lst)), dtype = np.float)
        Vmax1 = np.zeros((len(day_lst)), dtype = np.float)
        Vmin1 = np.zeros((len(day_lst)), dtype = np.float)
        V = np.zeros((len(day_lst), 1, cMF.nrow, cMF.ncol), dtype = np.float)
        for ii, t in enumerate(day_lst):
            MM = conv_fact*h5_MM['finf_d'][t,:,:]
            V[ii,0,:,:] = np.ma.masked_values(MM, cMF.hnoflo, atol = 0.09)
            Vmin1[ii] = np.ma.min(np.ma.masked_values(V[ii,0,:,:], cMF.hnoflo, atol = 0.09))
            Vmax1[ii] = np.ma.max(np.ma.masked_values(V[ii,0,:,:], cMF.hnoflo, atol = 0.09))
        for ii, t in enumerate(day_lst):
            Vmax[ii] = np.ma.max(Vmax1)
            Vmin[ii] = np.ma.min(Vmin1)
        MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' $(mm/day)$'), msg = 'no flux', plt_title = ('OUT_MM_'+i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax, Vmin = Vmin, contours = ctrsMM, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = maskAllL_tmp, animation = animation)
        MMplot.plotLAYER(timesteps = day_lst, Date = Date_lst, JD = JD_lst, ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' $(mm/day)$'), msg = 'no flux', plt_title = ('OUT_MM_%s1'%i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax1, Vmin = Vmin1, contours = ctrsMM, ntick = ntick, points = obs4map, hnoflo = cMF.hnoflo, mask = maskAllL_tmp, animation = animation)

        h5_MM.close()
        del V, MM, t, Vmax, Vmin
        del day_lst, flxlbl, i, i1, h_diff_surf
    del hmaxMF, hminMF, hmin, hdiff, cbcmax_d, cbcmin_d

# CALIBRATION CRITERIA
# RMSE, RSR, Nash-Sutcliffe efficiency NSE, Pearson's correlation coefficient r
# Moriasi et al, 2007, ASABE 50(3):885-900
print '\nComputing calibration criteria at observation points...'
h5_MM = h5py.File(h5_MM_fn, 'r')
for o_ref in obs_list:
    for o in obs.keys():
        if o == o_ref:
            i = obs.get(o)['i']
            j = obs.get(o)['j']
            SOILzone_tmp = gridSOIL[i,j]-1
            nsl = _nsl[SOILzone_tmp]
            l_obs = obs.get(o)['lay']
            obs_h = obs.get(o)['obs_h']
            obs_SM = obs.get(o)['obs_SM']
            # TODO insert Ro            
            if obs_SM != []:
                obs_SM_tmp = obs_SM[:,HYindex[1]:HYindex[-1]]
            else:
                obs_SM_tmp = []
            MM_S = h5_MM['MM_S'][HYindex[1]:HYindex[-1],i,j,0:nsl,index_S.get('iSsoil_pc')]
            if obs_h != []:
                obs_h_tmp = obs_h[l_obs,HYindex[1]:HYindex[-1]]
            else:
                obs_h_tmp = []
            h_MF = h_MF_m[HYindex[1]:HYindex[-1],:,i,j]
            rmseHEADS_tmp, rmseSM_tmp, rsrHEADS_tmp, rsrSM_tmp, nseHEADS_tmp, nseSM_tmp, rHEADS_tmp, rSM_tmp = cMF.cPROCESS.compCalibCrit(MM_S, h_MF, obs_SM_tmp, obs_h_tmp, cMF.hnoflo, o, nsl)
            if rmseHEADS_tmp <> None:
                rmseHEADS.append(rmseHEADS_tmp)
                rsrHEADS.append(rsrHEADS_tmp)
                nseHEADS.append(nseHEADS_tmp)
                rHEADS.append(rHEADS_tmp)
                obslstHEADS.append(o)
            if rmseSM_tmp <> None:
                rmseSM.append(rmseSM_tmp)
                rsrSM.append(rsrSM_tmp)
                nseSM.append(nseSM_tmp)
                rSM.append(rSM_tmp)
                obslstSM.append(o)
            del rmseHEADS_tmp, rmseSM_tmp, rsrHEADS_tmp, rsrSM_tmp, nseHEADS_tmp, nseSM_tmp, rHEADS_tmp, rSM_tmp, h_MF, MM_S
for cc, (calibcritSM, calibcritHEADS, calibcrit, title, calibcritSMmax, calibcritHEADSmax, ymin, units) in enumerate(zip([rmseSM, rsrSM, nseSM, rSM], [rmseHEADS, rsrHEADS, nseHEADS, rHEADS], ['RMSE', 'RSR', 'NSE', 'r'], ['Root mean square error', 'Root mean square error - observations standard deviation ratio', 'Nash-Sutcliffe efficiency', "Pearson's correlation coefficient"], [rmseSMmax, None, 1.0, 1.0], [rmseHEADSmax, None, 1.0, 1.0], [0, 0, None, -1.0], [['($m$)', '($\%%wc$)'], ['',''], ['',''], ['','']])):
    #try:
    MMplot.plotCALIBCRIT(calibcritSM = calibcritSM, calibcritSMobslst = obslstSM, calibcritHEADS = calibcritHEADS, calibcritHEADSobslst = obslstHEADS, plt_export_fn = os.path.join(MM_ws_out, '__plt_calibcrit%s.png'% calibcrit), plt_title = 'Calibration criteria between simulated and observed state variables\n%s'%title, calibcrit = calibcrit, calibcritSMmax = calibcritSMmax, calibcritHEADSmax = calibcritHEADSmax, ymin = ymin, units = units, hnoflo = cMF.hnoflo)
    #except:
    #    print 'Error in exporting %s at obs. pt. %s' % (calibcrit, obs_list[cc])
print '-------\nRMSE/RSR/NSE/r averages of the obs. pts. (except catch.)'
try:
    for cc, (rmse, rsr, nse, r, obslst, msg) in enumerate(zip([rmseSM,rmseHEADS],[rsrSM,rsrHEADS],[nseSM,nseHEADS],[rSM,rHEADS],[obslstSM,obslstHEADS],['SM (all layers): %.1f %% /','h: %.2f m /'])):
        if obslst[0] == 'catch.' and len(rmse)> 1:
            rmseaverage = list(itertools.chain.from_iterable(rmse[1:]))
            rsraverage = list(itertools.chain.from_iterable(rsr[1:]))
            nseaverage = list(itertools.chain.from_iterable(nse[1:]))
            raverage = list(itertools.chain.from_iterable(r[1:]))
            numobs = len(obslst[1:])
        else:
            rmseaverage = list(itertools.chain.from_iterable(rmse))
            rsraverage = list(itertools.chain.from_iterable(rsr))
            nseaverage = list(itertools.chain.from_iterable(nse))
            raverage = list(itertools.chain.from_iterable(r))
            numobs = len(obslst)
        if len(rmse)> 1:
            rmseaverage = sum(rmseaverage)/float(len(rmseaverage))
            rsraverage = sum(rsraverage)/float(len(rsraverage))
            nseaverage = sum(nseaverage)/float(len(nseaverage))
            raverage = sum(raverage)/float(len(raverage))
            msg = '%s %s' % (msg, '%.2f / %.2f / %.2f (%d obs. points)')
            print msg % (rmseaverage, rsraverage, nseaverage, raverage, numobs)
except:
    print 'Error! Check observations data.'
    
timeendExport = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
durationExport=(timeendExport-timestartExport)
durationTotal = (timeendExport-timestart)

# final report of successful run
print '\n##############\nMARMITES executed successfully!\n%s' % mpl.dates.DateFormatter.format_data(fmt_DH, mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat()))
if MMsoil_yn > 0:
    print '\nLOOP %d/%d' % (LOOP-1, ccnum)
    for txt in msg_end_loop:
        print txt
print '\n%d stress periods, %d days' % (cMF.nper,sum(cMF.perlen))
print '%d layers x %d rows x %d cols (%d cells)' % (cMF.nlay, cMF.nrow, cMF.ncol, cMF.nrow*cMF.ncol)
print '%d MM active cells in total' % (sum(ncell_MM))
l = 1
for n in ncell_MF:
    print  'LAYER %d' % l
    print '%d MF active cells' % (n)
    print '%d MM active cells' % (ncell_MM[l-1])
    l += 1
print ('\nApproximate run times:')
if MMsurf_yn > 0:
    print ('MARMITES surface: %s minute(s) and %.1f second(s)') % (str(int(durationMMsurf*24.0*60.0)), (durationMMsurf*24.0*60.0-int(durationMMsurf*24.0*60.0))*60)
if MMsoil_yn > 0:
    print ('MARMITES soil zone: %s minute(s) and %.1f second(s)') % (str(int(durationMMsoil*24.0*60.0)), (durationMMsoil*24.0*60.0-int(durationMMsoil*24.0*60.0))*60)
if MF_yn == 1:
    print ('MODFLOW: %s minute(s) and %.1f second(s)') % (str(int(durationMF*24.0*60.0)), (durationMF*24.0*60.0-int(durationMF*24.0*60.0))*60)
print ('Export: %s minute(s) and %.1f second(s)') % (str(int(durationExport*24.0*60.0)), (durationExport*24.0*60.0-int(durationExport*24.0*60.0))*60)
print ('Total: %s minute(s) and %.1f second(s)') % (str(int(durationTotal*24.0*60.0)), (durationTotal*24.0*60.0-int(durationTotal*24.0*60.0))*60)
print ('\nOutput written in folder: \n%s\n##############\n') % MM_ws_out

if verbose == 0:
    sys.stdout = s
    report.close()
    print '##############\nMARMITES terminated normally!\n%s\n##############' % mpl.dates.DateFormatter.format_data(fmt_DH, mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat()))
    del s