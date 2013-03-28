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

import sys, os, traceback, h5py, shutil
import matplotlib as mpl
if mpl.get_backend!='agg':
    mpl.use('agg')
import matplotlib.pyplot as plt
plt.ioff()
import numpy as np
import startMARMITESsurface as startMMsurf
import MARMITESsoil_v3 as MMsoil
import MARMITESprocess_v3 as MMproc
import ppMODFLOW_flopy_v3 as ppMF
import MARMITESplot_v3 as MMplot
import CreateColors

# TODO verify if thickness of MF layer 1 > thickness of soil

#####################################

def minmax(min_, max_, ctrs_):
    if max_ == min_:
        if max_ < 10E-9:
            max_ = 1.0
            min_ = -1.0
        else:
            max_ *= 1.15
            min_ *= 0.85
        ctrs_ = False
    else:
        ctrs_ = ctrs_
    return min_, max_, ctrs_

#####################################

timestart = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
fmt = mpl.dates.DateFormatter('%Y-%m-%d %H:%M')
print '\n##############\nMARMITESv0.3 started!\n%s\n##############' % mpl.dates.DateFormatter.format_data(fmt, timestart)

# workspace (ws) definition
# read input file (called _input.ini in the MARMITES workspace
# the first character on the first line has to be the character used to comment
# the file can contain any comments as the user wish, but the sequence of the input has to be respected
# 00_TESTS\MARMITESv3_r13c6l2  00_TESTS\r40c20  00_TESTS\r20c40  r130c60l2   r130c60l2new
# SARDON2013  CARRIZAL3 CARRIZAL3newera LAMATA LaMata_new
MM_ws = r'E:\00code_ws\00_TESTS\MARMITESv3_r13c6l2'
MM_ini_fn = '__inputMM_v3.ini'

inputFile = MMproc.readFile(MM_ws,MM_ini_fn)

fmt = mpl.dates.DateFormatter('%Y%m%d%H%M')
MM_ws_out = os.path.join(MM_ws,'out_%s'%mpl.dates.DateFormatter.format_data(fmt, timestart))
f = 0
while not os.path.exists(MM_ws_out):
    if f>0:
        MM_ws_out1 = os.path.join(MM_ws,'out_%s_%s' % (MM_ws_out,mpl.dates.DateFormatter.format_data(fmt, timestart),f))
    else:
        MM_ws_out1 = MM_ws_out
    os.makedirs(MM_ws_out1)
    f += 1
del MM_ws_out1

l=0
try:
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
    # read observations?
    l += 1
    plt_out_obs = int(inputFile[l].strip())
    l += 1
    plt_WB_unit = inputFile[l].strip()
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
    chunks = int(inputFile[l].strip())
    if MMsoil_yn == 1:
        MF_yn = 1
    if MMsurf_plot == 1:
        plt_out = 0
        plt_out_obs = 0
        print "\nYou required the MMsurf plots to appear on the screen. Due to backends limitation, MM and MF plots were disabled. Run again MM with MMsurf_plot = 0 to obtain the MM and MF plots."
except:
    raise SystemExit('\nFATAL ERROR!\nType error in the input file %s\nError description:\n%s' % (MM_ini_fn, traceback.print_exc(file=sys.stdout)))
del inputFile

shutil.copy2(os.path.join(MM_ws,MM_ini_fn), os.path.join(MM_ws_out,'_%s%s'% (mpl.dates.DateFormatter.format_data(fmt, timestart), MM_ini_fn)))
shutil.copy2(os.path.join(MM_ws, MF_ws,MF_ini_fn), os.path.join(MM_ws_out,'_%s%s'% (mpl.dates.DateFormatter.format_data(fmt, timestart), MF_ini_fn)))

if verbose == 0:
#capture interpreter output to be written in to a report file
    report_fn = os.path.join(MM_ws_out,'_MM_00report_%s.txt' % (mpl.dates.DateFormatter.format_data(fmt, timestart)))
    del fmt
    print '\nECHO OFF (no screen output).\nSee the report of the MM-MF run in file:\n%s\n' % report_fn
    s = sys.stdout
    report = open(report_fn, 'w')
    sys.stdout = report
    print '\n##############\nMARMITES started!\n%s\n##############' % mpl.dates.num2date(timestart).isoformat()[:19]
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
#try:
# #############################
# ###  MARMITES surface  ######
# #############################
print'\n##############'
print 'MARMITESsurf RUN'
if MMsurf_yn>0:
    durationMMsurf = 0.0
    timestartMMsurf = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    if irr_yn == 1 :
        outMMsurf_fn = startMMsurf.MMsurf(MMsurf_ws, inputFile_TS_fn, inputFile_PAR_fn, outputFILE_fn, MM_ws, outMMsurf_fn, MMsurf_plot, inputFile_TSirr_fn)
    else:
        outMMsurf_fn = startMMsurf.MMsurf(MMsurf_ws, inputFile_TS_fn, inputFile_PAR_fn, outputFILE_fn, MM_ws, outMMsurf_fn, MMsurf_plot)
    timeendMMsurf = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    durationMMsurf=(timeendMMsurf-timestartMMsurf)
inputFile = MMproc.readFile(MM_ws,outMMsurf_fn)
l=0
VegName = []
Zr = []
kTu_min = []
kTu_n = []
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
        kTu_min.append(float(line[v]))
    l += 1
    line = inputFile[l].split()
    for v in range(NVEG):
        kTu_n.append(float(line[v]))
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
        kTu_min_c = []
        kTu_n_c = []
        line = inputFile[l].split()
        for c in range(NCROP):
            Zr_c.append(float(line[c]))
        Zr_c = np.asarray(Zr_c)
        l += 1
        line = inputFile[l].split()
        for c in range(NCROP):
            kTu_min_c.append(float(line[c]))
        kTu_min_c = np.asarray(kTu_min_c)
        l += 1
        line = inputFile[l].split()
        for c in range(NCROP):
            kTu_n_c.append(float(line[c]))
        kTu_n_c = np.asarray(kTu_n_c)
        l += 1
        input_dSP_crop_irr_fn = str(inputFile[l].strip())
except:
    raise SystemExit('\nFATAL ERROR!\Error in reading file [%s].\nError description:\n%s' % (outMMsurf_fn, traceback.print_exc(file=sys.stdout)))
del inputFile
numDays = len(MMproc.readFile(MM_ws, inputDate_fn))

# #############################
# ###  READ MODFLOW CONFIG ####
# #############################
print'\n##############'
print 'MODFLOW initialization'
durationMF = 0.0
timestartMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
cMF = ppMF.MF(MM_ws, MF_ws, MF_ini_fn, xllcorner, yllcorner, numDays = numDays)
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
    raise SystemExit('FATAL ERROR!\nDefine the length unit in the MODFLOW ini file!\n (see USGS Open-File Report 00-92)')
    # TODO if lenuni!=2 apply conversion factor to delr, delc, etc...
if cMF.laytyp[0]==0:
    raise SystemExit('FATAL ERROR!\nThe first layer cannot be confined type!\nChange your parameter laytyp in the MODFLOW lpf package.\n(see USGS Open-File Report 00-92)')
if cMF.itmuni != 4:
    raise SystemExit('FATAL ERROR!\nTime unit is not in days!')
else:
    itmuni_str = 'day'
ncell_MF = []
ncell_MM = []
iboundBOL = np.ones(np.array(cMF.ibound).shape, dtype = bool)
mask_tmp = np.zeros((cMF.nrow, cMF.ncol), dtype = int)
mask = []
for l in range(cMF.nlay):
    ncell_MF.append((np.asarray(cMF.ibound)[:,:,l] != 0).sum())
    ncell_MM.append((np.asarray(cMF.iuzfbnd) == l+1).sum())
    iboundBOL[:,:,l] = (np.asarray(cMF.ibound)[:,:,l] != 0)
    mask.append(np.ma.make_mask(iboundBOL[:,:,l]-1))
    mask_tmp += (np.asarray(cMF.ibound)[:,:,l] != 0)
maskAllL = (mask_tmp == 0)
del iboundBOL, mask_tmp
timeendMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
durationMF +=  timeendMF-timestartMF

# #############################
# ### MF time processing
# #############################
# if required by user, compute nper, perlen,etc based on RF analysis in the METEO zones
if cMF.timedef >= 0:
    if isinstance(cMF.nper, str):
        try:
            perlenmax = int(cMF.nper.split()[1].strip())
        except:
            raise SystemExit('\nFATAL ERROR!\nError in nper format of the MODFLOW ini file!\nError description:\n%s' % traceback.print_exc(file=sys.stdout))
    if irr_yn == 0:
        cMF.ppMFtime(MM_ws, MF_ws, inputDate_fn, inputZON_dSP_RF_veg_fn, inputZON_dSP_RFe_veg_fn, inputZON_dSP_PT_fn,input_dSP_LAI_veg_fn, inputZON_dSP_PE_fn, inputZON_dSP_E0_fn, NMETEO, NVEG, NSOIL)
    else:
        cMF.ppMFtime(MM_ws, MF_ws, inputDate_fn, inputZON_dSP_RF_veg_fn, inputZON_dSP_RFe_veg_fn, inputZON_dSP_PT_fn, input_dSP_LAI_veg_fn, inputZON_dSP_PE_fn, inputZON_dSP_E0_fn, NMETEO, NVEG, NSOIL, inputZON_dSP_RF_irr_fn, inputZON_dSP_RFe_irr_fn, inputZON_dSP_PT_irr_fn, input_dSP_crop_irr_fn, NFIELD)

print'\n##############'
print 'MARMITESsoil initialization'
MM_SOIL = MMsoil.SOIL(hnoflo = cMF.hnoflo)
MM_SATFLOW = MMsoil.SATFLOW()

# READ input ESRI ASCII rasters
print "\nImporting ESRI ASCII files to initialize MARMITES..."
gridMETEO = cMF.MM_PROCESS.inputEsriAscii(grid_fn = gridMETEO_fn, datatype = int)
gridSOIL = cMF.MM_PROCESS.inputEsriAscii(grid_fn = gridSOIL_fn, datatype = int)
gridSOILthick = cMF.MM_PROCESS.inputEsriAscii(grid_fn = gridSOILthick_fn, datatype = float)
gridSsurfhmax = cMF.MM_PROCESS.inputEsriAscii(grid_fn = gridSsurfhmax_fn, datatype = float)
gridSsurfw = cMF.MM_PROCESS.inputEsriAscii(grid_fn = gridSsurfw_fn, datatype = float)
if irr_yn == 1:
    gridIRR = cMF.MM_PROCESS.inputEsriAscii(grid_fn = gridIRR_fn, datatype = int)
    if gridIRR.max() > NFIELD:
        raise SystemExit('\nFATAL ERROR!\nThere is more fields in the asc file than in the MMsurf file!')

# READ input time series and parameters
if irr_yn == 0:
    gridVEGarea, RF_veg_zoneSP, E0_zonesSP, PT_veg_zonesSP, RFe_veg_zonesSP, LAI_veg_zonesSP, PE_zonesSP = cMF.MM_PROCESS.inputSP(                           NMETEO                   = NMETEO,
                                NVEG                     = NVEG,
                                NSOIL                    = NSOIL,
                                nper                     = cMF.nper,
                                inputZON_SP_RF_veg_fn    = cMF.inputZON_SP_RF_veg_fn,
                                inputZON_SP_RFe_veg_fn   = cMF.inputZON_SP_RFe_veg_fn,
                                inputZON_SP_LAI_veg_fn   = cMF.inputZON_SP_LAI_veg_fn,
                                inputZON_SP_PT_fn        = cMF.inputZON_SP_PT_fn,
                                inputZON_SP_PE_fn        = cMF.inputZON_SP_PE_fn,
                                inputZON_SP_E0_fn        = cMF.inputZON_SP_E0_fn
                                )
else:
    gridVEGarea, RF_veg_zoneSP, E0_zonesSP, PT_veg_zonesSP, RFe_veg_zonesSP, LAI_veg_zonesSP, PE_zonesSP, RF_irr_zoneSP, RFe_irr_zoneSP, PT_irr_zonesSP, crop_irr_SP = cMF.MM_PROCESS.inputSP(
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
_nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _Su_ini, _Ks = cMF.MM_PROCESS.inputSoilParam(MM_ws = MM_ws, SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
_nslmax = max(_nsl)
for l in range(NSOIL):
    _slprop[l] = np.asarray(_slprop[l])

# compute thickness, top and bottom elevation of each soil layer
cMF.elev = np.asarray(cMF.elev)
cMF.top = cMF.elev - gridSOILthick
cMF.botm = np.asarray(cMF.botm)
for l in range(cMF.nlay):
    cMF.botm[:,:,l] -= gridSOILthick
botm_l0 = np.asarray(cMF.botm)[:,:,0]

# create MM array
h5_MM_fn = os.path.join(MM_ws,'_h5_MM.h5')
# indexes of the HDF5 output arrays
index = {'iRF':0, 'iPT':1, 'iPE':2, 'iRFe':3, 'iSsurf':4, 'iRo':5, 'iEXF':6, 'iEsurf':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSsurf':13, 'iETg':14, 'iETsoil':15, 'iSsoil_pc':16, 'idSsoil':17, 'iinf':18, 'iHEADScorr':19, 'idgwt':20, 'iuzthick':21}
index_S = {'iEsoil':0, 'iTsoil':1,'iSsoil_pc':2, 'iRp':3, 'iRexf':4, 'idSsoil':5, 'iSsoil':6, 'iSAT':7, 'iMB_l':8}

# READ observations time series (heads and soil moisture)
print "\nReading observations time series (hydraulic heads and soil moisture)..."
obs, obs_list = cMF.MM_PROCESS.inputObs(MM_ws  = MM_ws,
                              inputObs_fn      = inputObs_fn,
                              inputObsHEADS_fn = inputObsHEADS_fn,
                              inputObsSM_fn    = inputObsSM_fn,
                              inputDate        = cMF.inputDate,
                              _nslmax          = _nslmax,
                              nlay             = cMF.nlay
                              )
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
header='Date,MF_SP,veg_crop,RF,E0,PT,PE,RFe,I,' + Esoil_str + Tsoil_str + 'Eg,Tg,ETg,WEL_MF,Esurf,' + Ssoil_str + Ssoilpc_str + dSsoil_str + 'dSsurf,Ssurf,Ro,GW_EXF,' + Rp_str + Rexf_str + 'R_MF,hSATFLOW,hMF,hMFcorr,hmeas,dgwt,' + Smeasout + MB_str + 'MB\n'
if cMF.uzf_yn == 1:
    cMF.uzf_obs(obs = obs)

# EXPORT INPUT MAPS
if plt_input == 1:
    ibound = np.asarray(cMF.ibound)
    print'\n##############'
    print 'Exporting input maps...'
    i_lbl = 1
    top_tmp = np.zeros((cMF.nrow, cMF.ncol, cMF.nlay), dtype = np.float)
    top_tmp[:,:,0] = cMF.top
    for l in range(1,cMF.nlay):
        top_tmp[:,:,l] = cMF.botm[:,:,l-1]
    lst = [cMF.elev, top_tmp, cMF.botm, cMF.thick, np.asarray(cMF.strt), gridSOILthick, gridSsurfhmax, gridSsurfw, np.asarray(cMF.hk_actual), np.asarray(cMF.ss_actual), np.asarray(cMF.sy_actual), np.asarray(cMF.vka_actual)]
    lst_lbl = ['elev', 'top', 'botm', 'thick', 'strt', 'gridSOILthick', 'gridSsurfhmax', 'gridSsurfw', 'hk', 'Ss', 'Sy', 'vka']
    lst_lblCB = ['Elev.', 'Aq. top - top', 'Aq. bot. - botm', 'Aq. thick.', 'Init. heads - strt', 'Soil thick.', 'Stream max. heigth', 'Stream width', 'Horizontal hydraulic cond. - hk', 'Specific storage - Ss', 'Specific yield - Sy', 'Vertical hydraulic cond. - vka']
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
        lst_lblCB.append('Epsilon - eps')
        lst.append(np.asarray(cMF.thts_actual))
        lst_lbl.append('thts')
        lst_lblCB.append('Sat. water content - thts')
        if cMF.vks_actual != 0.0:
            lst.append(np.asarray(cMF.vks_actual))
            lst_lbl.append('vks')
            lst_lblCB.append('Sat. vert. hydraulic conductivity - vks')
        try:
            lst.append(np.asarray(cMF.thti_actual))
            lst_lbl.append('thti')
            lst_lblCB.append('Initial water content - thti')
        except:
            pass
        try:
            lst.append(np.asarray(cMF.thtr_actual))
            lst_lbl.append('thtr')
            lst_lblCB.append('Residual water content - thtr')
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
        #print lst_lbl[i]
        V = []
        if l.shape == (cMF.nrow, cMF.ncol, cMF.nlay) or l.shape == (cMF.nrow, cMF.ncol):
            for L in range(cMF.nlay):
                try:
                    V.append(np.ma.masked_values(np.ma.masked_array(l[:,:,L], mask[L]), cMF.hnoflo, atol = 0.09))
                    nplot = cMF.nlay
                except:
                    V.append(np.ma.masked_values(np.ma.masked_array(l, maskAllL), cMF.hnoflo, atol = 0.09))
                    nplot = 1
        else:
            if l.shape != ():
                for L in range(cMF.nlay):
                    V.append(np.ma.masked_values(np.ma.masked_array(l[L]*ibound[:,:,L], mask[L]), cMF.hnoflo, atol = 0.09))
                nplot = cMF.nlay
            else:
                V.append(np.ma.masked_values(np.ma.masked_array(l*np.invert(maskAllL), maskAllL), cMF.hnoflo, atol = 0.09))
                nplot = 1
        if lst_lbl[i] == 'Ss' or lst_lbl[i] == 'Sy' or lst_lbl[i] == 'hk' or lst_lbl[i] == 'thts' or lst_lbl[i] == 'thti' or lst_lbl[i] == 'thtr' or lst_lbl[i] == 'drn_cond' or lst_lbl[i] == 'ghb_cond':
            fmt = '%.3g'
        else:
            fmt = '%.2f'
        if lst_lbl[i] == 'Ss' or lst_lbl[i] == 'eps':
           CBlabel = lst_lblCB[i] + ' ([-])'
        elif lst_lbl[i] == 'hk' or lst_lbl[i] == 'drn_cond' or lst_lbl[i] == 'ghb_cond':
            CBlabel = lst_lblCB[i] + ' (%s/%s)' % (lenuni_str, itmuni_str)
        elif lst_lbl[i] == 'Sy' or lst_lbl[i] == 'thts' or lst_lbl[i] == 'thti' or lst_lbl[i] == 'thtr':
            CBlabel = lst_lblCB[i] + ' ($L^3$/$L^{-3}$)'
        elif lst_lbl[i] == 'vka':
            CBlabel = lst_lblCB[i] + ' ([-] or %s/%s, layvka = %s)' % (lenuni_str, itmuni_str, cMF.layvka)
        else:
            CBlabel = lst_lblCB[i] + ' (m)'
        if lst_lbl[i] == 'elev' or lst_lbl[i] == 'top' or lst_lbl[i] == 'botm':
            Vmax = elev_max
            Vmin = elev_min
        else:
            Vmax = np.ma.max(V) #float(np.ceil(np.ma.max(V)))  #
            Vmin = np.ma.min(V) #float(np.floor(np.ma.min(V))) #
        Vmin_tmp, Vmax_tmp, ctrsV_tmp = minmax(Vmin, Vmax, ctrsMM)
        MMplot.plotLAYER(SP = 0, Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = nplot, V = V,  cmap = plt.cm.gist_rainbow_r, CBlabel = CBlabel, msg = '', plt_title = '_INPUT_%03d_%s' % (i_lbl,lst_lbl[i]), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsV_tmp, Vmax = Vmax_tmp, Vmin = Vmin_tmp, ntick = ntick, axisbg = 'white', fmt = fmt)
        i_lbl += 1
    del V, lst, lst_lbl, nplot, Vmax, Vmin, Vmax_tmp, Vmin_tmp, ctrsV_tmp

    Vmax = 100.0
    Vmin = 0.0
    for v in range(NVEG):
        V = [np.ma.masked_values(np.ma.masked_array(gridVEGarea[v,:,:], maskAllL), cMF.hnoflo, atol = 0.09)]
        V_lbl = 'veg%02d_%s' %(v+1, VegName[v])
        V_lblCB = 'Frac. area of veg #%02d (%s)' %(v+1, VegName[v])
        #print V_lbl
        MMplot.plotLAYER(SP = 0, Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.gist_rainbow_r, CBlabel = '%s %s' % (V_lblCB, ' (%)'), msg = '', plt_title = '_INPUT_%03d_%s'% (i_lbl,V_lbl), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = False, Vmax = Vmax, Vmin = Vmin, ntick = ntick, axisbg = 'white', fmt = '%3.1f')
        i_lbl += 1
    del V, V_lbl, Vmax, Vmin

    lst = [ibound, gridSOIL, gridMETEO]
    lst_lbl = ['ibound', 'gridSOIL', 'gridMETEO']
    lst_lblCB = ['MF cell type - ibound', 'Soil type', 'Meteo. zone']
    if irr_yn == 1:
        lst.append(gridIRR)
        lst_lbl.append('gridIRR')
        lst_lblCB.append('Irrigation plots')
    if cMF.uzf_yn == 1:
        lst.append(np.asarray(cMF.iuzfbnd))
        lst_lbl.append('iuzfbnd')
        lst_lblCB.append('Recharge layer - iuzfbnd')
    for i, l in enumerate(lst):
        #print lst_lbl[i]
        V = []
        for L in range(cMF.nlay):
            try:
                V.append(np.ma.masked_array(l[:,:,L], mask[L]))
                nplot = cMF.nlay
            except:
                V.append(np.ma.masked_array(l, maskAllL))
                nplot = 1
        Vmax = np.ma.max(V) #float(np.ceil(np.ma.max(V)))  #
        Vmin = np.ma.min(V) #float(np.floor(np.ma.min(V))) #
        Vmin_tmp, Vmax_tmp, ctrsV_tmp = minmax(Vmin, Vmax, ctrsMM)
        MMplot.plotLAYER(SP = 0, Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = nplot, V = V,  cmap = plt.cm.gist_rainbow_r, CBlabel = lst_lblCB[i], msg = '', plt_title = '_INPUT_%03d_%s'% (i_lbl,lst_lbl[i]), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = False, Vmax = Vmax_tmp, Vmin = Vmin_tmp, ntick = ntick, axisbg = 'white')
        i_lbl += 1
    del V, lst, lst_lbl, lst_lblCB, nplot, Vmax, Vmin, Vmax_tmp, Vmin_tmp, ctrsV_tmp

# #############################
# ### 1st MODFLOW RUN with initial user-input recharge
# #############################
if MF_yn == 1 :
    timestartMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    print'\n##############'
    print 'MODFLOW RUN (initial user-input fluxes)'
    if verbose == 0:
        print '\n--------------'
        sys.stdout = s
        report.close()
        s = sys.stdout
        report = open(report_fn, 'a')
        sys.stdout = report
    cMF.ppMF(MM_ws, MF_ws, MF_ini_fn, finf_MM = (h5_MM_fn, 'finf'), wel_MM = (h5_MM_fn, 'ETg'), report = report, verbose = verbose, chunks = chunks, numDays = numDays)
    timeendMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    durationMFtmp =  timeendMF-timestartMF
    durationMF +=  durationMFtmp
    print 'MF run time: %02.fmn%02.fs' % (int(durationMFtmp*24.0*60.0), (durationMFtmp*24.0*60.0-int(durationMFtmp*24.0*60.0))*60)
    del durationMFtmp

if os.path.exists(cMF.h5_MF_fn):
    try:
        h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
    except:
        raise SystemExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.\nError description:\n%s' % traceback.print_exc(file=sys.stdout))
    # heads format is : timestep, nrow, ncol, nlay
    # cbc format is: (kstp), kper, textprocess, nrow, ncol, nlay
    cbc_nam = []
    cbc_uzf_nam = []
    ncells_package = []
    for c in h5_MF['cbc_nam']:
        cbc_nam.append(c.strip())
    if cMF.uzf_yn == 1:
        for c in h5_MF['cbc_uzf_nam']:
            cbc_uzf_nam.append(c.strip())
    elif cMF.rch_yn == 1:
        imfRCH = cbc_nam.index('RECHARGE')
        raise SystemExit('\nFATAL ERROR!\nMM has to be run together with the UZF1 package of MODFLOW 2005/ MODFLOW NWT, thus the RCH package should be desactivacted!\nExisting MM.')
    imfSTO = cbc_nam.index('STORAGE')
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

##except StandardError:  #Exception
##    raise SystemExit('\nFATAL ERROR!\nAbnormal MM run interruption in the initialization!\nError description:\n%s' % traceback.print_exc(file=sys.stdout))

# #############################
# 2nd phase : MM/MF loop #####
# #############################
h_diff_surf = None
#try:
if MMsoil_yn > 0:
    durationMMsoil = 0.0
    h_pSP = 0
    h_pSP_all = 0
    LOOP = 0
    LOOPlst = [LOOP]
    h_diff = [1000]
    h_diff_log = [1]
    h_diff_all = [1000]
    h_diff_all_log = [1]
    loopdry = 0
    plt_ConvLoop_fn = os.path.join(MM_ws_out, '_MM_00plt_MM_MF_ConvLoop.png')
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

    while abs(h_diff[LOOP]) > convcrit or abs(h_diff_all[LOOP]) > convcritmax:
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
        _nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _Su_ini, _Ks = cMF.MM_PROCESS.inputSoilParam(MM_ws = MM_ws, SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
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
            MM_SOIL.run(_nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _Su_ini, botm_l0, _Ks,
                          gridSOIL, gridSOILthick, cMF.elev*1000.0, gridMETEO,
                          index, index_S, gridSsurfhmax, gridSsurfw,
                          RF_veg_zoneSP, E0_zonesSP, PT_veg_zonesSP, RFe_veg_zonesSP, PE_zonesSP, gridVEGarea,
                          LAI_veg_zonesSP, Zr, kTu_min, kTu_n, NVEG,
                          cMF, conv_fact, h5_MF, h5_MM, irr_yn
                          )
        else:
            MM_SOIL.run(_nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _Su_ini, botm_l0, _Ks,
                          gridSOIL, gridSOILthick, cMF.elev*1000.0, gridMETEO,
                          index, index_S, gridSsurfhmax, gridSsurfw,
                          RF_veg_zoneSP, E0_zonesSP, PT_veg_zonesSP, RFe_veg_zonesSP, PE_zonesSP, gridVEGarea,
                          LAI_veg_zonesSP, Zr, kTu_min, kTu_n, NVEG,
                          cMF, conv_fact, h5_MF, h5_MM, irr_yn,
                          RF_irr_zoneSP, PT_irr_zonesSP, RFe_irr_zoneSP,
                          crop_irr_SP, gridIRR,
                          Zr_c, kTu_min_c, kTu_n_c, NCROP
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
        if LOOP <2:
            msg_end_loop.append('Initial average heads:\n%.3f m' % h_diff[LOOP])
        else:
            msg_end_loop.append('Heads diff. from previous conv. loop: %.3f m' % h_diff[LOOP])
            msg_end_loop.append('Maximum heads difference:             %.3f m' % h_diff_all[LOOP])
        if h_MF_average == 0.0:
            loopdry += 1
            if loopdry > 1:
                print '\nWARNING: first layer of the model DRY twice successively!\nLoop break, correct your MARMITES input value.'
                timeendMMloop = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
                durationMMloop = timeendMMloop-timestartMMloop
                print '\nMM run time: %02.fmn%02.fs' % (int(durationMMloop*24.0*60.0), (durationMMloop*24.0*60.0-int(durationMMloop*24.0*60.0))*60)
                durationMMsoil += durationMMloop
                break
            else:
                print '\nWARNING: first layer of the model DRY!'
        elif abs(h_diff[LOOP]) < convcrit and abs(h_diff_all[LOOP]) < convcritmax:
            msg_end_loop.append('Successfull convergence between MARMITES and MODFLOW!\n(Conv. criterion = %.4f and conv. crit. max. = %.4f)' % (convcrit, convcritmax))
            for txt in msg_end_loop:
                print txt
            timeendMMloop = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
            durationMMloop = timeendMMloop-timestartMMloop
            print '\nMM run time: %02.fmn%02.fs' % (int(durationMMloop*24.0*60.0), (durationMMloop*24.0*60.0-int(durationMMloop*24.0*60.0))*60)
            durationMMsoil += durationMMloop
            break
        elif LOOP>ccnum:
            msg_end_loop.append('No convergence between MARMITES and MODFLOW!\n(Conv. criterion = %.4f and conv. crit. max. = %.4f)' % (convcrit, convcritmax))
            for txt in msg_end_loop:
                print txt
            timeendMMloop = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
            durationMMloop = timeendMMloop-timestartMMloop
            print '\nMM run time: %02.fmn%02.fs' % (int(durationMMloop*24.0*60.0), (durationMMloop*24.0*60.0-int(durationMMloop*24.0*60.0))*60)
            durationMMsoil += durationMMloop
            break
        del h_MF_average
        for txt in msg_end_loop:
            print txt

        timeendMMloop = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        durationMMloop = timeendMMloop-timestartMMloop
        print '\nMM run time: %02.fmn%02.fs' % (int(durationMMloop*24.0*60.0), (durationMMloop*24.0*60.0-int(durationMMloop*24.0*60.0))*60)
        durationMMsoil += durationMMloop

        # MODFLOW RUN with MM-computed recharge
        timestartMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        print'\n##############'
        print 'MODFLOW RUN (MARMITES fluxes)'
        if verbose == 0:
            print '\n--------------'
            sys.stdout = s
            report.close()
            s = sys.stdout
            report = open(report_fn, 'a')
            sys.stdout = report
        cMF.ppMF(MM_ws, MF_ws, MF_ini_fn, finf_MM = (h5_MM_fn, 'finf'), wel_MM = (h5_MM_fn, 'ETg'), report = report, verbose = verbose, chunks = chunks, numDays = numDays)
        try:
            h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
        except:
            raise SystemExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.\nError description:\n%s' % traceback.print_exc(file=sys.stdout))

        timeendMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        durationMFtmp =  timeendMF-timestartMF
        durationMF +=  durationMFtmp
        print '\nMF run time: %02.fmn%02.fs' % (int(durationMFtmp*24.0*60.0), (durationMFtmp*24.0*60.0-int(durationMFtmp*24.0*60.0))*60)
        del durationMFtmp

    # #############################
    # ###  END CONVERGENCE LOOP ###
    # #############################

    # export loop plot
#        print'\n##############'
#        print 'Exporting plot of the convergence loop...'
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
    # ### MODFLOW RUN with MM-computed recharge
    # #############################
    if MF_yn == 1:
        try:
            h5_MF.close()
        except:
            pass
        timestartMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        print'\n##############'
        print 'MODFLOW RUN (MARMITES fluxes after conv. loop)'
        if verbose == 0:
            print '\n--------------'
            sys.stdout = s
            report.close()
            s = sys.stdout
            report = open(report_fn, 'a')
            sys.stdout = report
        cMF.ppMF(MM_ws, MF_ws, MF_ini_fn, finf_MM = (h5_MM_fn, 'finf'), wel_MM = (h5_MM_fn, 'ETg'), report = report, verbose = verbose, chunks = chunks, numDays = numDays)
        timeendMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        durationMFtmp =  timeendMF-timestartMF
        durationMF +=  durationMFtmp
        print '\nMF run time: %02.fmn%02.fs' % (int(durationMFtmp*24.0*60.0), (durationMFtmp*24.0*60.0-int(durationMFtmp*24.0*60.0))*60)
        del durationMFtmp

##except StandardError:  #Exception
##    raise SystemExit('\nFATAL ERROR!\nAbnormal MM run interruption in the MM/MF loop!\nError description:\n%s' % traceback.print_exc(file=sys.stdout))

# #############################
# 3rd phase : export results #####
# #############################

#try:

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
        raise SystemExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.\nError description:\n%s' % traceback.print_exc(file=sys.stdout))
    if cMF.rch_yn == 1:
        cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'RCH_d', conv_fact = conv_fact, index = imfRCH)
    cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'STO_d', conv_fact = conv_fact, index = imfSTO)
    if cMF.drn_yn == 1:
        cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'DRN_d', conv_fact = conv_fact, index = imfDRN)
    if cMF.wel_yn == 1:
        cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'WEL_d', conv_fact = conv_fact, index = imfWEL)
    if cMF.ghb_yn == 1:
        cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'GHB_d', conv_fact = conv_fact, index = imfGHB)
    cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'heads', ds_name_new = 'heads_d', conv_fact = conv_fact)
    if cMF.uzf_yn == 1:
        cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc_uzf', ds_name_new = 'RCH_d', conv_fact = conv_fact, index = imfRCH)
        cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc_uzf', ds_name_new = 'EXF_d', conv_fact = conv_fact, index = imfEXF)
    h5_MF.close()

SP_d = np.ones(sum(cMF.perlen), dtype = int)
t = 0
for n in range(cMF.nper):
    for x in range(cMF.perlen[n]):
        SP_d[t] = n + 1
        t += 1

if os.path.exists(cMF.h5_MF_fn):
    index_cbc = [imfSTO]
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
    cbc_DRN = cbc_STO = cbc_RCH = cbc_WEL = np.zeros((sum(cMF.perlen), cMF.nrow, cMF.ncol, cMF.nlay))
    imfDRN = imfSTO = imfRCH = imfWEL = 0

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
    for L in range(cMF.nlay):
        V.append(np.ma.masked_array(h_diff_surf[h_diff_n,:,:,L], mask[L]))
        Vmax.append(np.ma.max(V[L]))
        Vmin.append(np.ma.min(V[L]))
    Vmax = max(Vmax) #float(np.ceil(max(Vmax)))
    Vmin = min(Vmin) #float(np.floor(min(Vmin)))
    Vmin_tmp, Vmax_tmp, ctrs_tmp = minmax(Vmin, Vmax, ctrsMF)
    # TODO JD and Date are not correct since h_diff_n is # stress periods and not # of days (same in the plots of MF and MM)
    MMplot.plotLAYER(SP = h_diff_n, Date = cMF.inputDate[h_diff_n], JD = cMF.JD[h_diff_n], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = ('(m)'), msg = 'no value', plt_title = '_HEADSmaxdiff_ConvLoop', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax_tmp, Vmin = Vmin_tmp, contours = ctrs_tmp, ntick = ntick)
    del h_diff_n, Vmin, Vmax, Vmin_tmp, Vmax_tmp

# exporting sm computed by MM for PEST (smp format)
if os.path.exists(h5_MM_fn):
    try:
        h5_MM = h5py.File(h5_MM_fn, 'r')
    except:
        raise SystemExit('\nFATAL ERROR!\nInvalid MARMITES HDF5 file. Run MARMITES and/or MODFLOW again.\nError description:\n%s' % traceback.print_exc(file=sys.stdout))
    outPESTsmMM = open(os.path.join(MF_ws,'sm_MM4PEST.smp'), 'w')
    for o_ref in obs_list:
        for o in obs.keys():
            if o == o_ref:
                i = obs.get(o)['i']
                j = obs.get(o)['j']
                l = obs.get(o)['lay']
                obs_S = obs.get(o)['obs_S']
                if obs.get(o)['obs_sm_yn'] == 1:
                    cMF.MM_PROCESS.smMMname.append(o)
                MM_S = h5_MM['MM_S'][:,i,j,:,:]
                cMF.MM_PROCESS.ExportResultsMM4PEST(i, j, cMF.inputDate, _nslmax, MM_S, index_S, obs_S, o)
    # write PEST smp file with MM output
    inputFile = MMproc.readFile(MM_ws,inputObs_fn)
    ind = 0
    for i in range(len(inputFile)):
        line = inputFile[i].split()
        name = line[0]
        for j in cMF.MM_PROCESS.smMMname:
            if j == name:
                for l in cMF.MM_PROCESS.smMM[ind]:
                    outPESTsmMM.write(l)
        ind += 1
    outPESTsmMM.close()

if plt_out == 1 or plt_out_obs == 1:
    # computing max and min values in MF fluxes for plotting
    hmax = []
    hmin = []
    hmaxMM = []
    hminMM = []
    hdiff = []
    cbcmax_d = []
    cbcmin_d = []
    axefact = 1.05
    facTim = 365
    if os.path.exists(cMF.h5_MF_fn):
        # TODO missing STOuz (however this is not very relevant since these fluxes should not be the bigger in magnitude)
        try:
            h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
        except:
            raise SystemExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.\nError description:\n%s' % traceback.print_exc(file=sys.stdout))
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
                    for t,col in enumerate(cbc_RCH[:,row,:,l]):
                        try:
                            if plt_out_obs == 1:
                                j = list(col).index(RCHmax)
                                x = cMF.delc[row]*row + xllcorner
                                y = cMF.delr[j]*j + yllcorner
                                obs['PzRCHmax'] = {'x':x,'y':y, 'i': row, 'j': j, 'lay': l, 'hi':999, 'h0':999, 'RC':999, 'STO':999, 'outpathname':os.path.join(MM_ws,'_MM_0PzRCHmax.txt'), 'obs_h':[], 'obs_h_yn':0, 'obs_S':[], 'obs_sm_yn':0}
                                obs_list.append('PzRCHmax')
                                del x, y
                            print 'row %d, col %d, layer %d and day %d (%s)' % (row + 1, list(col).index(RCHmax) + 1, l+1, t, mpl.dates.num2date(cMF.inputDate[t] + 1.0).isoformat()[:10])
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
        GWTD = cMF.elev[np.newaxis,:,:,np.newaxis] - h_MF_m
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
            raise SystemExit('\nFATAL ERROR!\nInvalid MARMITES HDF5 file. Run MARMITES and/or MODFLOW again.\nError description:\n%s' % traceback.print_exc(file=sys.stdout))
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
                hmaxMF_tmp = np.ma.max(h_MF_m[:,i,j,l].flatten())
                hminMF_tmp = np.ma.min(h_MF_m[:,i,j,l].flatten())
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
    # plot SOIL/GW balance at the catchment scale
    # #################################################
    tTgmin = -1
    if os.path.exists(h5_MM_fn):
        try:
            h5_MM = h5py.File(h5_MM_fn, 'r')
        except:
            raise SystemExit('\nFATAL ERROR!\nInvalid MARMITES HDF5 file. Run MARMITES and/or MODFLOW again.\nError description:\n%s' % traceback.print_exc(file=sys.stdout))
        # indexes of the HDF5 output arrays
        #  index = {'iRF':0, 'iPT':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iEXF':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSs':13, 'iETg':14, 'iETu':15, 'iSu_pc':16, 'idSu':17, 'iinf':18, 'iHEADScorr':19, 'idgwt':20, 'iuzthick':21}
        # index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'iRexf':4, 'idSu':5, 'iSu':6, 'iSAT':7, 'iMB_l':8}
        flxlbl       = ['RF', 'I', 'RFe', 'dSsurf', 'Ro', 'Esurf', 'dSsoil', 'EXF']
        flxlbl1      = ['Esoil', 'Tsoil']
        flxlbl2      = ['ETsoil', 'Eg']
        flxlbl3      = ['Tg']
        flxlbl3a     = ['ETg']
        flxlbl4      = ['Rp']
        sign         = [   1,  -1,     1,     1,   -1,   -1,    1,     1, -1, -1, -1, -1, -1, -1, -1]
        flxlst       = []
        flx_Cat_TS   = []
        flxmax_d     = []
        flxmin_d     = []
        flxlbl_CATCH = []
        TopSoilAverage = np.ma.masked_array(cMF.elev*1000.0, maskAllL).sum()*.001/sum(ncell_MM)
        for i in flxlbl:
            flxlbl_CATCH.append(i)
            i = 'i'+i
            array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
            flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
            flxmax_d.append(np.ma.max(flx_tmp))
            flxmin_d.append(np.ma.min(flx_tmp))
            flxlst.append(facTim*(flx_tmp.sum())/sum(cMF.perlen)/sum(ncell_MM))
            array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
            flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
            del flx_tmp, array_tmp, array_tmp1
        for i in flxlbl1:
            flxlbl_CATCH.append(i)
            flxlbl.append(i)
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
            flxlst.append(facTim*flx_tmp1/sum(cMF.perlen)/sum(ncell_MM))
            flx_Cat_TS.append(array_tmp2/sum(ncell_MM))
            del flx_tmp1, array_tmp2
        for i in flxlbl2:
            flxlbl_CATCH.append(i)
            flxlbl.append(i)
            i = 'i'+i
            array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
            flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
            flxmax_d.append(np.ma.max(flx_tmp))
            flxmin_d.append(np.ma.min(flx_tmp))
            flxlst.append(facTim*(flx_tmp.sum())/sum(cMF.perlen)/sum(ncell_MM))
            array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
            flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
            del flx_tmp, array_tmp, array_tmp1
        for i in flxlbl3:
            flxlbl_CATCH.append(i)
            flxlbl.append(i)
            i = 'i'+i
            if cMF.wel_yn == 1:
                array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
                flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
                flxmax_d.append(np.ma.max(flx_tmp))
                Tg_min = np.ma.min(flx_tmp)
                flxmin_d.append(Tg_min)
                flxlst.append(facTim*(flx_tmp.sum())/sum(cMF.perlen)/sum(ncell_MM))
                if abs(Tg_min) > 1E-6:
                    print '\nTg negative (%.2f) observed at:' % Tg_min
                    for row in range(cMF.nrow):
                        for t,col in enumerate(flx_tmp[:,row,:]):
                            try:
                                print 'row %d, col %d and day %d' % (row + 1, list(col).index(Tg_min) + 1, t + 1)
                                tTgmin = t
                                if plt_out_obs == 1:
                                    obs['PzTgmin'] = {'x':999,'y':999, 'i': row, 'j': list(col).index(Tg_min), 'lay': 0, 'hi':999, 'h0':999, 'RC':999, 'STO':999, 'outpathname':os.path.join(MM_ws,'_MM_0PzTgmin.txt'), 'obs_h':[], 'obs_h_yn':0, 'obs_S':[], 'obs_sm_yn':0}
                                    obs_list.append('PzTgmin')
                                    try:
                                        hmin.append(hmin[0])
                                    except:
                                        hmin.append(999.9)
                            except:
                                pass
                array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
                flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
                del flx_tmp, array_tmp, array_tmp1
            else:
                flxlst.append(0.0)
                flx_Cat_TS.append(0.0)
        for i in flxlbl3a:
            flxlbl_CATCH.append(i)
            flxlbl.append(i)
            i = 'i'+i
            array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
            flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
            flxmax_d.append(np.ma.max(flx_tmp))
            flxmin_d.append(np.ma.min(flx_tmp))
            flxlst.append(facTim*flx_tmp.sum()/sum(cMF.perlen)/sum(ncell_MM))
            array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
            flx_Cat_TS.append(-1.0*np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
            del flx_tmp, array_tmp, array_tmp1
        for i in flxlbl4:
            flxlbl.append(i)
            i = 'i'+i
            array_tmp = h5_MM['finf'][:,:,:]
            flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
            flxmax_d.append(conv_fact*np.ma.max(flx_tmp))
            flxmin_d.append(conv_fact*np.ma.min(flx_tmp))
            inf = facTim*conv_fact*(flx_tmp.sum())/sum(cMF.perlen)/sum(ncell_MM)
            flxlst.append(inf)
            del flx_tmp, array_tmp
        for i in ['Ssurf', 'PE', 'PT', 'inf']:
            flxlbl_CATCH.append(i)
            i = 'i'+i
            array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
            array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
            flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
            del array_tmp, array_tmp1
        plt_exportCATCH_fn = os.path.join(MM_ws_out, '_plt_0CATCHMENT.png')
        plt_exportCATCH_txt_fn = os.path.join(MM_ws_out, '_plt_0CATCHMENT.txt')
        plt_titleCATCH = 'Time serie of fluxes averaged over the whole catchment'
        InMM  = flxlst[0] + flxlst[3] + flxlst[6] + flxlst[7]
        OutMM = flxlst[1] + flxlst[4] + flxlst[5] + flxlst[10] + flxlst[14]
        MB_MM = 100*(InMM - OutMM)/((InMM + OutMM)/2)
        # TODO delete next line after resolving MB_MM error and uncomment the previous one
        #MB_MM = InMM - OutMM
        # ADD SM averaged
        flxlbl_CATCH.append('Su')
        array_tmp = h5_MM['MM'][:,:,:,index.get('iSsoil_pc')]
        array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
        flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
        del array_tmp, array_tmp1
        h5_MM.close()
        flxmax_d = float(np.ceil(np.ma.max(flxmax_d)))
        flxmin_d = float(np.floor(np.ma.min(flxmin_d)))
        flxmax = axefact*float(np.ceil(max(flxlst)))
        flxmin = axefact*float(np.floor(min(flxlst)))
        for l,(x,y) in enumerate(zip(flxlst, sign)):
            flxlst[l] = x*y
        del flxlbl1, flxlbl2, flxlbl3, flxlbl3a, sign
        if os.path.exists(cMF.h5_MF_fn):
            plt_export_fn = os.path.join(MM_ws_out, '_plt_0CATCHMENT_SOILandGWbalances.png')
            # compute UZF_STO and store GW_RCH
            flxlbl.append('UZ_STO')
            rch_tmp = 0
            flxlst_tmp = []
            h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
            cbc_RCH = h5_MF['RCH_d']
            array_tmp2 = np.zeros((sum(cMF.perlen)), dtype = np.float)
            for l in range(cMF.nlay):
                array_tmp = cbc_RCH[:,:,:,l]
                array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
                array_tmp2 += np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)
                rch_tmp1 = facTim*(array_tmp.sum())/sum(cMF.perlen)/sum(ncell_MM)
                flxlst_tmp.append(rch_tmp1)
                rch_tmp += rch_tmp1
            flxlbl_CATCH.append('R')
            flx_Cat_TS.append(array_tmp2/sum(ncell_MM))
            del array_tmp, array_tmp1, array_tmp2, rch_tmp1, cbc_RCH
            flxlst.append(inf - rch_tmp)
            del rch_tmp, inf
            InUZF = -flxlst[14] - flxlst[15]
            OutUZF = 0
            InMF = 0
            OutMF = 0
            for l in range(cMF.nlay):
                # ADD heads averaged
                flxlbl_CATCH.append('h_MF_L%d' % (l+1))
                array_tmp = h_MF_m[:,:,:,l]
                array_tmp1 = np.sum(array_tmp, axis = 1)
                flx_Cat_TS.append(np.sum(array_tmp1, axis = 1)/ncell_MF[l])
                del array_tmp, array_tmp1
                # ADD depth GWT
                flxlbl_CATCH.append('DGWT_L%d' % (l+1))
                flx_Cat_TS.append(flx_Cat_TS[-1] - TopSoilAverage)
            for l in range(cMF.nlay):
                # GW_STO
                cbc_STO = h5_MF['STO_d'][:,:,:,l]
                flxlst.append(facTim*(cbc_STO.sum()/sum(cMF.perlen)/sum(ncell_MM)))  # -1*
                InMF += flxlst[-1]
                flxlbl_CATCH.append('$\Delta$Sg_L%d' % (l+1))
                array_tmp1 = np.sum(np.ma.masked_values(cbc_STO, cMF.hnoflo, atol = 0.09), axis = 1)
                flx_Cat_TS.append(-1.0*np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
                del cbc_STO, array_tmp1
                # GW_RCH
                flxlst.append(flxlst_tmp[l])
                OutUZF += flxlst_tmp[l]
                InMF += flxlst_tmp[l]
                # EXF
                cbc_EXF = h5_MF['EXF_d'][:,:,:,l]
                flxlst.append(facTim*(cbc_EXF.sum()/sum(cMF.perlen)/sum(ncell_MM)))
                OutMF += -flxlst[-1]
                del cbc_EXF
                if cMF.drn_yn == 1:
                    cbc_DRN = h5_MF['DRN_d'][:,:,:,l]
                    if cMF.drncells[l]>0:
                        flxlst.append(facTim*(cbc_DRN.sum()/sum(cMF.perlen)/sum(ncell_MM)))
                        OutMF += -flxlst[-1]
                        flxlbl_CATCH.append('DRN_L%d' % (l+1))
                        array_tmp1 = np.sum(np.ma.masked_values(cbc_DRN, cMF.hnoflo, atol = 0.09), axis = 1)
                        flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
                        del array_tmp1
                    else:
                        flxlst.append(0.0)
                    del cbc_DRN
                if cMF.wel_yn == 1:
                    cbc_WEL = h5_MF['WEL_d'][:,:,:,l]
                    if ncell_MM[l]>0:
                        flxlst.append(facTim*(cbc_WEL.sum()/sum(cMF.perlen)/sum(ncell_MM)))
                        OutMF += -1.0*flxlst[-1]
                    else:
                        flxlst.append(0.0)
                    del cbc_WEL
                if cMF.ghb_yn == 1:
                    cbc_GHB = h5_MF['GHB_d'][:,:,:,l]
                    if cMF.ghbcells[l] > 0:
                        flxlst.append(facTim*(cbc_GHB.sum()/sum(cMF.perlen)/sum(ncell_MM)))
                        OutMF += -flxlst[-1]
                        flxlbl_CATCH.append('GHB_L%d' % (l+1))
                        array_tmp1 = np.sum(np.ma.masked_values(cbc_GHB, cMF.hnoflo, atol = 0.09), axis = 1)
                        flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
                        del array_tmp1
                    else:
                        flxlst.append(0.0)
                    del cbc_GHB
                flxlbl.append('GW_' + cbc_nam[index_cbc[0]] + '_L' + str(l+1))
                for x in range(len(index_cbc_uzf)):
                    flxlbl.append('GW_' + cbc_uzf_nam[index_cbc_uzf[x]] + '_L' + str(l+1))
                for x in range(1,len(index_cbc)):
                    flxlbl.append('GW_' + cbc_nam[index_cbc[x]] + '_L' + str(l+1))
            flxmax1 = float(np.ceil(max(flxlst)))
            flxmin1 = float(np.floor(min(flxlst)))
            flxmax = axefact*max(flxmax, flxmax1)
            flxmin = axefact*min(flxmin, flxmin1)
            del flxlst_tmp, flxmax1, flxmin1
            h5_MF.close()
            MB_UZF = 100.0*(InUZF - OutUZF)/((InUZF+OutUZF)/2.0)
            MB_MF = 100.0*(InMF - OutMF)/((InMF+OutMF)/2.0)
            plt_title = 'MARMITES and MODFLOW water balance for the whole catchment\nMass balance error: MM = %1.2f%%, UZF = %1.2f%%, MF = %1.2f%%' % (MB_MM, MB_UZF, MB_MF)
            header_tmp = ['MM_MB','UZF_MB','MF_MB']
            MB_tmp = [MB_MM, MB_UZF,MB_MF]
            MMplot.plotTIMESERIES_CATCH(cMF.inputDate, flx_Cat_TS, flxlbl_CATCH, plt_exportCATCH_fn, plt_titleCATCH, hmax = hmaxMF, hmin = hminMF, cMF = cMF)
        else:
            MMplot.plotTIMESERIES_CATCH(cMF.inputDate, flx_Cat_TS, flxlbl_CATCH, plt_exportCATCH_fn, plt_titleCATCH, hmax = hmaxMF, hmin = hminMF)
            plt_export_fn = os.path.join(MM_ws_out, '_plt_0CATCHMENT_SOILbalance.png')
            plt_title = 'MARMITES water balance for the whole catchment\nMass balance error: MM = %1.2f%%' % (MB_MM)
            header_tmp = ['MM_MB']
            MB_tmp = [MB_MM]
        # export average time serie of fluxes in txt
        plt_exportCATCH_txt = open(os.path.join(plt_exportCATCH_txt_fn), 'w')
        flxlbl_CATCH_str = 'Date,'
        flxlbl_CATCH_str += flxlbl_CATCH[0]
        for e in (flxlbl_CATCH[1:]):
            flxlbl_CATCH_str += ',' + e
        plt_exportCATCH_txt.write(flxlbl_CATCH_str)
        plt_exportCATCH_txt.write('\n')
        for t in range(len(cMF.inputDate)):
            flx_Cat_TS_str = str(flx_Cat_TS[0][t])
            for e in (flx_Cat_TS[1:]):
                flx_Cat_TS_str += ',' + str(e[t])
            out_line = mpl.dates.num2date(cMF.inputDate[t]).isoformat()[:10] + ',' + flx_Cat_TS_str
            for l in out_line:
                plt_exportCATCH_txt.write(l)
            plt_exportCATCH_txt.write('\n')
        plt_exportCATCH_txt.close()
        del flx_Cat_TS, flx_Cat_TS_str, out_line
        colors_flx = CreateColors.main(hi=0, hf=180, numbcolors = len(flxlbl))
        MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title, fluxmax = flxmax, fluxmin = flxmin, unit = plt_WB_unit)
        plt_export_txt = open(os.path.join(MM_ws_out, '_plt_0CATCHMENT_SOILandGWbalances.txt'), 'w')
        flxlbl_str = flxlbl[0]
        for e in (flxlbl[1:] + header_tmp):
            flxlbl_str += ',' + e
        flxlst_str = str(flxlst[0])
        for e in (flxlst[1:] + MB_tmp):
            flxlst_str += ',' + str(e)
        plt_export_txt.write(flxlbl_str)
        plt_export_txt.write('\n')
        plt_export_txt.write(flxlst_str)
        plt_export_txt.close()
        del flxlst, header_tmp, MB_tmp, flxlst_str
        del flxmax_d, flxmin_d, flxmax, flxmin, flxlbl_CATCH, TopSoilAverage, i
        del plt_exportCATCH_fn, plt_exportCATCH_txt_fn, plt_titleCATCH, InMM, OutMM, InMF, OutMF, MB_MM, MB_UZF, plt_export_fn, plt_title
    # no MM, plot only MF balance if MF exists
    elif os.path.exists(cMF.h5_MF_fn):
        plt_export_fn = os.path.join(MM_ws_out, '_plt_0CATCHMENT_GWbalances.png')
        flxlbl = []
        flxlst = []
        InMF = 0
        OutMF = 0
        for l in range(cMF.nlay):
            h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
            # GW_STO
            cbc_STO = h5_MF['STO_d']
            flxlst.append(facTim*(cbc_STO[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))  # *-1
            InMF += flxlst[-1]
            del cbc_STO
            # GW_RCH
            cbc_RCH = h5_MF['RCH_d']
            flxlst.append(facTim*(cbc_RCH[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))
            InMF += flxlst[-1]
            del cbc_RCH
            # EXF
            cbc_EXF = h5_MF['EXF_d']
            flxlst.append(facTim*(cbc_EXF[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))
            OutMF += flxlst[-1]
            del cbc_EXF
            if cMF.drn_yn == 1:
                cbc_DRN = h5_MF['DRN_d']
                if cMF.drncells[l]>0:
                    flxlst.append(facTim*(cbc_DRN[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))
                    OutMF += flxlst[-1]
                else:
                    flxlst.append(0.0)
                del cbc_DRN
            if cMF.wel_yn == 1:
                cbc_WEL = h5_MF['WEL_d']
                if ncell_MM[l]>0:
                    flxlst.append(-1.0*facTim*(cbc_WEL[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))
                    OutMF += flxlst[-1]
                else:
                    flxlst.append(0.0)
                del cbc_WEL
            if cMF.ghb_yn == 1:
                cbc_GHB = h5_MF['GHB_d']
                if cMF.ghbcells[l] > 0:
                    flxlst.append(facTim*(cbc_GHB[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))
                    OutMF += flxlst[-1]
                else:
                    flxlst.append(0.0)
                del cbc_GHB
            h5_MF.close()
            flxlbl.append('GW_' + cbc_nam[index_cbc[0]] + '_L' + str(l+1))
            for x in range(len(index_cbc_uzf)):
                flxlbl.append('GW_' + cbc_uzf_nam[index_cbc_uzf[x]] + '_L' + str(l+1))
            for x in range(1,len(index_cbc)):
                flxlbl.append('GW_' + cbc_nam[index_cbc[x]] + '_L' + str(l+1))
        colors_flx = CreateColors.main(hi=0, hf=180, numbcolors = len(flxlbl))
        MB_MF = 100.0*(InMF - OutMF)/((InMF+OutMF)/2.0)
        # TODO cbcmax and cbcmin represent the max and min of the whole MF fluxes time series, not the max and min of the fluxes average
        plt_title = 'MODFLOW water balance for the whole catchment\nMass balance error: MF = %1.2f%%' % (MB_MF)
        MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title, fluxmax = cbcmax_d, fluxmin = cbcmin_d, unit = plt_WB_unit)
        del flxlst, MB_MF, InMF, OutMF

    # #################################################
    # EXPORT AT OBSERVATION POINTS
    # exporting MM time series results to ASCII files and plots at observations cells
    # #################################################
    if plt_out_obs == 1 and os.path.exists(h5_MM_fn) and os.path.exists(cMF.h5_MF_fn):
        print '\nExporting output time series plots and ASCII files at observation points...'
        h5_MM = h5py.File(h5_MM_fn, 'r')
        h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
        colors_nsl = CreateColors.main(hi=00, hf=180, numbcolors = (_nslmax+1))
        x = 0
        flxlst = []
        plt_exportBAL_fn  = []
        plt_export_txt_fn = []
        plt_titleBAL      = []
        InOut_tmp             = []
        for o_ref in obs_list:
            for o in obs.keys():
                if o == o_ref:
                    i = obs.get(o)['i']
                    j = obs.get(o)['j']
                    l = obs.get(o)['lay']
                    obs_h = obs.get(o)['obs_h']
                    obs_S = obs.get(o)['obs_S']
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
                    h_satflow = MM_SATFLOW.run(cbc_RCH[:,i,j,0], float(obs.get(o)['hi']),float(obs.get(o)['h0']),float(obs.get(o)['RC']),float(obs.get(o)['STO']))
                    # export ASCII file at piezometers location
                    #TODO extract heads at piezo location and not center of cell
                    if obs_h != []:
                        obs_h_tmp = obs_h[0,:]
                    else:
                        obs_h_tmp = []
                    if cMF.wel_yn == 1:
                        cbc_WEL = -1.0*np.sum(np.ma.masked_values(h5_MF['WEL_d'][:,i,j,:], cMF.hnoflo, atol = 0.09), axis = 1)
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
                    cMF.MM_PROCESS.ExportResultsMM(i, j, cMF.inputDate, SP_d, _nslmax, MM, index, MM_S, index_S, cbc_RCH[:,i,j,0], cbc_WEL, h_satflow, h_MF_m[:,i,j,l], obs_h_tmp, obs_S, index_veg, outFileExport, o)
                    del cbc_WEL
                    outFileExport.close()
                    # plot time series results as plot
                    plt_title = 'Time serie of fluxes at observation point %s\ni = %d, j = %d, l = %d, x = %.2f, y = %.2f, elev. = %.2f, %s\nSm = %s, Sfc = %s, Sr = %s, Ks = %s, thick. = %.3f %s' % (o, i+1, j+1, l+1, obs.get(o)['x'], obs.get(o)['y'], cMF.elev[i,j], soilnam, _Sm[gridSOIL[i,j]-1], _Sfc[gridSOIL[i,j]-1], _Sr[gridSOIL[i,j]-1], _Ks[gridSOIL[i,j]-1], gridSOILthick[i,j], Tl)
                    # index = {'iRF':0, 'iPT':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iEXF':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSs':13, 'iETg':14, 'iETu':15, 'iSu_pc':16, 'idSu':17, 'iinf':18, 'iHEADScorr':19, 'idgwt':20, 'iuzthick':21}
                    # index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'iRexf':4, 'idSu':5, 'iSu':6, 'iSAT':7, 'iMB_l':8}
                    plt_export_fn = os.path.join(MM_ws_out, '_plt_0'+ o + '.png')
                    # def plotTIMESERIES(DateInput, P, PT, PE, Pe, dPOND, POND, Ro, Eu, Tu, Eg, Tg, S, dS, Spc, Rp, EXF, ETg, Es, MB, MB_l, dgwt, SAT, R, h_MF, h_MF_corr, h_SF, hobs, Sobs, Sm, Sr, hnoflo, plt_export_fn, plt_title, colors_nsl, hmax, hmin):
                    MMplot.plotTIMESERIES(
                    cMF.inputDate,
                    MM[:,index.get('iRF')],
                    MM[:,index.get('iPT')],
                    MM[:,index.get('iPE')],
                    MM[:,index.get('iRFe')],
                    -MM[:,index.get('idSsurf')],
                    MM[:,index.get('iSsurf')],
                    MM[:,index.get('iRo')],
                    MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iEsoil')],
                    MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iTsoil')],
                    MM[:,index.get('iEg')],
                    MM[:,index.get('iTg')],
                    MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iSsoil')],
                    MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('idSsoil')],
                    MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iSsoil_pc')],
                    MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iRp')],
                    MM[:,index.get('iEXF')],
                    -MM[:,index.get('iETg')],
                    MM[:,index.get('iEsurf')],
                    MM[:,index.get('iMB')],
                    MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iMB_l')],
                    MM[:,index.get('idgwt')],
                    MM[:,index.get('iuzthick')],
                    MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iSAT')],
                    cbc_RCH[:,i,j,0],
                    h_MF_m[:,i,j,:], MM[:,index.get('iHEADScorr')], h_satflow, obs_h_tmp, obs_S,
                    _Sm[gridSOIL[i,j]-1],
                    _Sr[gridSOIL[i,j]-1],
                    cMF.hnoflo,
                    plt_export_fn,
                    plt_title,
                    colors_nsl,
                    max(hmax), #hmax[x] + hdiff/2
                    min(hmin), #hmin[x] - hdiff/2
                    o,
                    cMF.elev[i,j],
                    cMF.nlay
                    )
                    x += 1
                    # plot water balance at each obs. cell
                    #flxlbl   = ['RF', 'I', 'dSs', 'Ro', 'Es', 'dSu', 'EXF']
                    #flxlbl1  = ['Eu', 'Tu']
                    #flxlbl2  = ['ETu', 'Eg']
                    #flxlbl3  = ['Tg']
                    #flxlbl3a = ['ETg']
                    #flxlbl4  = ['Rp']
                    flxlst.append([
                         facTim*(MM[:,index.get('iRF')].sum()/sum(cMF.perlen)),
                        -1.0*facTim*(MM[:,index.get('iI')].sum()/sum(cMF.perlen)),
                         facTim*(MM[:,index.get('iRFe')].sum()/sum(cMF.perlen)),
                        -1.0*facTim*(MM[:,index.get('idSsurf')].sum()/sum(cMF.perlen)),
                         -1.0*facTim*(MM[:,index.get('iRo')].sum()/sum(cMF.perlen)),
                        -1.0*facTim*(MM[:,index.get('iEsurf')].sum()/sum(cMF.perlen)),
                         facTim*(MM[:,index.get('idSsoil')].sum()/sum(cMF.perlen)),
                         facTim*(MM[:,index.get('iEXF')].sum()/sum(cMF.perlen)),
                        -1.0*facTim*(MM_S[:,:,index_S.get('iEsoil')].sum()/sum(cMF.perlen)),
                        -1.0*facTim*(MM_S[:,:,index_S.get('iTsoil')].sum()/sum(cMF.perlen)),
                        -1.0*facTim*(MM[:,index.get('iETsoil')].sum()/sum(cMF.perlen)),
                        -1.0*facTim*(MM[:,index.get('iEg')].sum()/sum(cMF.perlen)),
                        -1.0*facTim*(MM[:,index.get('iTg')].sum()/sum(cMF.perlen)),
                        -1.0*facTim*(MM[:,index.get('iETg')].sum()/sum(cMF.perlen)),
                        -1.0*facTim*conv_fact*((cMF.perlen*h5_MM['finf'][:,i,j]).sum()/sum(cMF.perlen))
                        ])
                    InMM = flxlst[-1][0] + flxlst[-1][3] + flxlst[-1][6] + flxlst[-1][7]
                    OutMM = -(flxlst[-1][1] + flxlst[-1][4] + flxlst[-1][5] + flxlst[-1][10] + flxlst[-1][14])
                    InOut_MM = InMM - OutMM
                    InOut_tmp.append([InOut_MM])
                    if isinstance(cMF.h5_MF_fn, str):
                        plt_exportBAL_fn.append(os.path.join(MM_ws_out, '_plt_0'+ o + '_SOILandGWbalances.png'))
                        plt_export_txt_fn.append(os.path.join(MM_ws_out, '_plt_0'+ o + '_SOILandGWbalances.txt'))
                        # compute UZF_STO and store GW_RCH
                        rch_tmp = 0
                        flxlst_tmp = []
                        for l in range(cMF.nlay):
                            rch_tmp1 = facTim*(cbc_RCH[:,i,j,l].sum()/sum(cMF.perlen))
                            flxlst_tmp.append(rch_tmp1)
                            rch_tmp += rch_tmp1
                        flxlst[-1].append(-rch_tmp + facTim*conv_fact*((cMF.perlen*h5_MM['finf'][:,i,j]).sum()/sum(cMF.perlen)))
                        del rch_tmp, rch_tmp1, cbc_RCH
                        InUZF = -flxlst[-1][14] - flxlst[-1][15]
                        OutUZF = 0
                        InMF = 0
                        OutMF = 0
                        for l in range(cMF.nlay):
                            # GW STO
                            cbc_STO = h5_MF['STO_d']
                            flxlst[-1].append(facTim*(cbc_STO[:,i,j,l].sum()/sum(cMF.perlen)))    # -1*
                            InMF += flxlst[-1][-1]
                            del cbc_STO
                            # GW_RCH
                            flxlst[-1].append(flxlst_tmp[l])
                            OutUZF += flxlst_tmp[l]
                            InMF += flxlst_tmp[l]
                            # EXF
                            h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
                            cbc_EXF = h5_MF['EXF_d']
                            flxlst[-1].append(facTim*(cbc_EXF[:,i,j,l].sum()/sum(cMF.perlen)))
                            OutMF += -flxlst[-1][-1]
                            del cbc_EXF
                            if cMF.drn_yn == 1:
                                cbc_DRN = h5_MF['DRN_d']
                                flxlst[-1].append(facTim*(cbc_DRN[:,i,j,l].sum()/sum(cMF.perlen)))
                                OutMF += -flxlst[-1][-1]
                                del cbc_DRN
                            if cMF.wel_yn == 1:
                                cbc_WEL = h5_MF['WEL_d']
                                flxlst[-1].append(facTim*(cbc_WEL[:,i,j,l].sum()/sum(cMF.perlen)))
                                OutMF += -1.0*flxlst[-1][-1]
                                del cbc_WEL
                            if cMF.ghb_yn == 1:
                                cbc_GHB = h5_MF['GHB_d']
                                flxlst[-1].append(facTim*(cbc_GHB[:,i,j,l].sum()/sum(cMF.perlen)))
                                OutMF += -flxlst[-1][-1]
                                del cbc_GHB
                        InOut_UZF = InUZF - OutUZF
                        InOut_tmp[-1].append(InOut_UZF)
                        InOut_MF = InMF - OutMF
                        InOut_tmp[-1].append(InOut_MF)
                        plt_title = 'MARMITES and MODFLOW water flux balance at observation point %s\ni = %d, j = %d, l = %d, x = %d, y = %d, %s\n\nMass balance (In - Out,  [mm]): MM = %1.2f, UZF = %1.2f, MF = %1.2f' % (o, i+1, j+1, l+1, obs.get(o)['x'], obs.get(o)['y'], soilnam, InOut_MM, InOut_UZF, InOut_MF)
                        plt_titleBAL.append(plt_title)
                        del flxlst_tmp, plt_title
                    else:
                        plt_exportBAL_fn.append(os.path.join(MM_ws_out, '_plt_0'+ o + '_SOILbalance.png'))
                        plt_export_txt_fn.append(os.path.join(MM_ws_out, '_plt_0' + o + '_SOILbalances.txt'), 'w')
                        plt_title = 'MARMITES water flux balance at observation point %s\ni = %d, j = %d, l = %d, x = %d, y = %d, %s\n\nMass balance (In - Out,  [mm]): MM = %1.2f' % (o, i+1, j+1, l+1, obs.get(o)['x'], obs.get(o)['y'], soilnam, InOut_MM)
                        plt_titleBAL.append(plt_title)
                        del plt_title
                    del obs_h, obs_S
        flxmax = axefact*float(np.ceil(np.asarray(flxlst).max()))
        flxmin = axefact*float(np.floor(np.asarray(flxlst).min()))
        for l, (lst, fn, title, fn_txt, InOut) in enumerate(zip(flxlst, plt_exportBAL_fn, plt_titleBAL, plt_export_txt_fn, InOut_tmp)):
            MMplot.plotGWbudget(flxlst = lst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = fn, plt_title = title, fluxmax = flxmax, fluxmin = flxmin, unit = plt_WB_unit)
            plt_export_txt = open(fn_txt, 'w')
            flxlst_str = str(lst[0])
            for e in (lst[1:]  + InOut):
                flxlst_str += ',' + str(e)
            plt_export_txt.write(flxlbl_str)
            plt_export_txt.write('\n')
            plt_export_txt.write(flxlst_str)
            plt_export_txt.close()
        del flxlst, flxlbl, flxlbl_str, InOut_tmp
        del h_satflow, MM, MM_S
        del i, j, l, SOILzone_tmp, outFileExport, nsl, soilnam, slprop, Tl, plt_export_fn, InMM, OutMM, InUZF, OutUZF, plt_titleBAL, plt_exportBAL_fn, colors_flx
        h5_MM.close()
        h5_MF.close()
        del obs

    # #################################################
    # PLOT SPATIAL MF and MM OUTPUT
    # #################################################
    if plt_out == 1:
        print '\nExporting output maps...'
        SP_lst = []
        Date_lst = []
        JD_lst = []
        SP = 0
        while SP < sum(cMF.perlen):  #len(h_MF_m):
            SP_lst.append(SP)
            Date_lst.append(cMF.inputDate[SP])
            JD_lst.append(cMF.JD[SP])
            SP += plt_freq
        if os.path.exists(cMF.h5_MF_fn):
            if tTgmin < 0 and tRCHmax > 0:
                lst = [sum(cMF.perlen) - 1, tRCHmax] #[len(h_MF_m)-1, tRCHmax]
            else:
                lst = [sum(cMF.perlen) - 1, tRCHmax, tTgmin] # [len(h_MF_m)-1, tRCHmax, tTgmin]
            for e in lst:
                SP_lst.append(e)
                Date_lst.append(cMF.inputDate[e])
                JD_lst.append(cMF.JD[e])

        # ##############################
        # #### MODFLOW OUTPUT ##########
        # ##############################
        if os.path.exists(cMF.h5_MF_fn):

            # ############################################
            # plot at specified SP
            # ############################################

            h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
            h5_MM = h5py.File(h5_MM_fn, 'r')
            cbc_RCH = h5_MF['RCH_d']
            if cMF.drn_yn == 1:
                cbc_DRN = h5_MF['DRN_d']
            if cMF.ghb_yn == 1:
                cbc_GHB = h5_MF['GHB_d']

            # plot for selected time step
            t = 0
            for SP in SP_lst:

                # plot heads [m]
                V = []
                for L in range(cMF.nlay):
                    V.append(h_MF_m[SP,:,:,L])
                hminMF_tmp, hmaxMF_tmp, ctrsMF_tmp = minmax(hminMF, hmaxMF, ctrsMF)
                MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation (m)', msg = 'DRY', plt_title = 'MF_HEADS', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsMF_tmp, Vmax = hmaxMF_tmp, Vmin = hminMF_tmp, ntick = ntick)

                # plot GWTD [m]
                for L in range(cMF.nlay):
                    V[L] = cMF.elev-V[L]
                GWTDmin_tmp, GWTDmax_tmp, ctrsGWTD_tmp = minmax(GWTDmin, GWTDmax, ctrsMF)
                MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'depth to groundwater table (m)', msg = 'DRY', plt_title = 'MF_GWTD', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsGWTD_tmp, Vmax = GWTDmax_tmp, Vmin = GWTDmin_tmp, ntick = ntick)

                # plot heads corrigidas [m]
                headscorr_m = [np.ma.masked_values(np.ma.masked_values(h5_MM['MM'][SP,:,:,19], cMF.hnoflo, atol = 0.09), cMF.hdry, atol = 1E+25)]
                hcorrmin_tmp, hcorrmax_tmp, ctrs_tmp = minmax(hcorrmin, hcorrmax, ctrsMF)
                MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = headscorr_m,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation (m)', msg = 'DRY', plt_title = 'MF_HEADScorr', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrs_tmp, Vmax = hcorrmax_tmp, Vmin = hcorrmin_tmp, ntick = ntick)

                # plot GWTD correct [m]
                GWTDcorr = []
                GWTDcorr.append(cMF.elev-headscorr_m[0])
                GWTDcorrmin_tmp, GWTDcorrmax_tmp, ctrsGWTDcorr_tmp = minmax(GWTDcorrmin, GWTDcorrmax, ctrsMF)
                MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = GWTDcorr,  cmap = plt.cm.Blues, CBlabel = 'depth to groundwater table (m)', msg = 'DRY', plt_title = 'MF_GWTDcorr', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsGWTDcorr_tmp, Vmax = GWTDcorrmax_tmp, Vmin = GWTDcorrmin_tmp, ntick = ntick)

                # plot GW drainage [mm]
                if cMF.drn_yn == 1:
                    V = []
                    for L in range(cMF.nlay):
                        V.append(np.ma.masked_array(cbc_DRN[SP,:,:,L], mask[L])*(-1.0))
                    DRNmin_tmp, DRNmax_tmp, ctrs_tmp = minmax(DRNmin, DRNmax, ctrsMF)
                    MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage (mm/day)', msg = '- no drainage', plt_title = 'MF_DRN', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = DRNmin_tmp, contours = ctrs_tmp, Vmax = DRNmax_tmp, ntick = ntick)

                # plot GHB [mm]
                if cMF.ghb_yn == 1:
                    V = []
                    for L in range(cMF.nlay):
                        V.append(np.ma.masked_array(cbc_GHB[SP,:,:,L], mask[L])*(-1.0))
                    GHBmin_tmp, GHBmax_tmp, ctrs_tmp = minmax(GHBmin, GHBmax, ctrsMF)
                    MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'GHB flux (mm/day)', msg = '- no flux', plt_title = 'MF_GHB', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = GHBmin_tmp, contours = ctrs_tmp, Vmax = GHBmax_tmp, ntick = ntick)

                # plot GW RCH [mm]
                V = []
                for L in range(cMF.nlay):
                    V.append(np.ma.masked_array(cbc_RCH[SP,:,:,L], mask[L]))
                RCHmin_tmp, RCHmax_tmp, ctrs_tmp = minmax(RCHmin, RCHmax, ctrsMF)
                MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater recharge (mm/day)', msg = '- no flux', plt_title = 'MF_RCH', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = RCHmin_tmp, contours = ctrs_tmp, Vmax = RCHmax_tmp, ntick = ntick)
                RCHmin_tmp1 = np.min(np.asarray(V))
                RCHmax_tmp1 = np.max(np.asarray(V))
                RCHmin_tmp1, RCHmax_tmp1, ctrs_tmp = minmax(RCHmin_tmp1, RCHmax_tmp1, ctrsMF)
                MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater recharge (mm/day)', msg = '- no flux', plt_title = 'MF_RCH1', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmin = RCHmin_tmp1, contours = ctrs_tmp, Vmax = RCHmax_tmp1, ntick = ntick)
                t += 1
                del V
            del t

            # ############################################
            # plot average for the whole simulated period
            # ############################################

            # plot heads average [m]
            V = []
            for L in range(cMF.nlay):
                V.append(np.ma.masked_array(np.sum(h_MF_m[:,:,:,L], axis = 0)/sum(cMF.perlen), mask[L]))
            MMplot.plotLAYER(SP = SP, Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation (m)', msg = 'DRY', plt_title = 'MF_average_HEADS', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsMF_tmp, Vmax = hmaxMF_tmp, Vmin = hminMF_tmp, ntick = ntick)
            del h_MF_m, hminMF_tmp, hmaxMF_tmp, ctrsMF_tmp

            # plot GWTD average [m]
            for L in range(cMF.nlay):
                V[L] = cMF.elev-V[L]
            MMplot.plotLAYER(SP = SP, Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'depth to groundwater table (m)', msg = 'DRY', plt_title = 'MF_average_GWTD', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsGWTD_tmp, Vmax = GWTDmax_tmp, Vmin = GWTDmin_tmp, ntick = ntick)
            del V, GWTDmax_tmp, GWTDmin_tmp

            # plot heads corrigidas average [m]
            headscorr_m = [np.ma.masked_values(np.ma.masked_values(h5_MM['MM'][SP,:,:,19], cMF.hnoflo, atol = 0.09), cMF.hdry, atol = 1E+25)]
            MMplot.plotLAYER(SP = SP, Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = headscorr_m,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation (m)', msg = 'DRY', plt_title = 'MF_average_HEADScorr', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrs_tmp, Vmax = hcorrmax_tmp, Vmin = hcorrmin_tmp, ntick = ntick)
            del hcorrmin_tmp, hcorrmax_tmp

            # plot GWTD correct average [m]
            GWTDcorr = []
            GWTDcorr.append(cMF.elev-headscorr_m[0])
            MMplot.plotLAYER(SP = SP, Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = GWTDcorr, cmap = plt.cm.Blues, CBlabel = 'depth to groundwater table (m)', msg = 'DRY', plt_title = 'MF_average_GWTDcorr', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, contours = ctrsGWTDcorr_tmp, Vmax = GWTDcorrmax_tmp, Vmin = GWTDcorrmin_tmp, ntick = ntick)
            del GWTDcorr, headscorr_m, GWTDcorrmax_tmp, GWTDcorrmin_tmp

            # plot GW drainage average [mm]
            if cMF.drn_yn == 1:
                V = []
                for L in range(cMF.nlay):
                    V.append(np.ma.masked_array(np.sum(cbc_DRN[:,:,:,L], axis = 0)/sum(cMF.perlen)*(-1.0), mask[L]))
                MMplot.plotLAYER(SP = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage (mm/day)', msg = '- no drainage', plt_title = 'MF_average_DRN', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = DRNmax_tmp, Vmin = DRNmin_tmp, contours = ctrs_tmp, ntick = ntick)
                del cbc_DRN, DRNmax_tmp, DRNmin_tmp, DRNmax, DRNmin

            # plot GHB average [mm]
            if cMF.ghb_yn == 1:
                V = []
                for L in range(cMF.nlay):
                    V.append(np.ma.masked_array(np.sum(cbc_GHB[:,:,:,L], axis = 0)/sum(cMF.perlen)*(-1.0), mask[L]))
                MMplot.plotLAYER(SP = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'GHB flux (mm/day)', msg = '- no flux', plt_title = 'MF_average_GHB', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = GHBmax_tmp, Vmin = GHBmin_tmp, contours = ctrs_tmp, ntick = ntick)
                del cbc_GHB, GHBmax_tmp, GHBmin_tmp, GHBmax, GHBmin

            # plot GW RCH average [mm]
            V = []
            for L in range(cMF.nlay):
                V.append(np.ma.masked_array(np.sum(cbc_RCH[:,:,:,L], axis = 0)/sum(cMF.perlen), mask[L]))
            MMplot.plotLAYER(SP = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater recharge (mm/day)', msg = '- no flux', plt_title = 'MF_average_RCH', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = RCHmax_tmp, Vmin = RCHmin_tmp, contours = ctrs_tmp, ntick = ntick)
            RCHmin_tmp1 = np.min(np.asarray(V))
            RCHmax_tmp1 = np.max(np.asarray(V))
            RCHmin_tmp1, RCHmax_tmp1, ctrs_tmp = minmax(RCHmin_tmp1, RCHmax_tmp1, ctrsMF)
            MMplot.plotLAYER(SP = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater recharge (mm/day)', msg = '- no flux', plt_title = 'MF_average_RCH1', MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = RCHmax_tmp1, Vmin = RCHmin_tmp1, contours = ctrs_tmp, ntick = ntick)
            h5_MF.close()
            del V, cbc_RCH, RCHmax_tmp, RCHmin_tmp

        # ##############################
        # #### MARMITES OUTPUT #########
        # ##############################
        if os.path.exists(h5_MM_fn):

            flxlbl = ['RF', 'RFe', 'I', 'EXF', 'dSsurf', 'Ro', 'Esurf', 'Eg', 'Tg', 'ETg', 'ETsoil', 'dSsoil']
            for i in flxlbl:
                # ############################################
                # plot average for the whole simulated period
                # ############################################
                i1 = 'i'+i
                h5_MM = h5py.File(h5_MM_fn, 'r')
                MM = h5_MM['MM'][:,:,:,index.get(i1)]
                h5_MM.close()
                V = [np.sum(np.ma.masked_values(MM, cMF.hnoflo, atol = 0.09), axis = 0)/sum(cMF.perlen)]
                V[0] = np.ma.masked_values(V[0], cMF.hnoflo, atol = 0.09)
                Vmax = np.ma.max(V) #float(np.ceil(np.ma.max(V)))  #
                Vmin = np.ma.min(V) #float(np.floor(np.ma.min(V))) #
                Vmin, Vmax, ctrs_tmp = minmax(Vmin, Vmax, ctrsMM)
                if i == 'dSsoil':
                    i_lbl = '$\Delta$Ssoil' #'$\Delta$S$_{u}'
                elif i == 'dSsurf':
                    i_lbl = '$\Delta$Ssurf'
                else:
                    i_lbl = i
                MMplot.plotLAYER(SP = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_average_' + i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
                del V

                # ############################################
                # plot for selected time step
                # ############################################
                t = 0
                for SP in SP_lst:
                    V = [np.ma.masked_values(MM[SP,:,:], cMF.hnoflo, atol = 0.09)]
                    MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_'+i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
                    t += 1
                del V, MM, t, Vmax, Vmin, ctrs_tmp

            flxlbl = ['Esoil', 'Tsoil','Rp']
            for i in flxlbl:

                # ############################################
                # plot average for the whole simulated period
                # ############################################
                i1 = 'i'+i
                h5_MM = h5py.File(h5_MM_fn, 'r')
                if i != 'Rp':
                    V = [np.zeros([cMF.nrow,cMF.ncol])]
                    for l in range(_nslmax):
                        MM = h5_MM['MM_S'][:,:,:,l,index_S.get(i1)]
                        V1 = np.sum(np.ma.masked_values(MM, cMF.hnoflo, atol = 0.09), axis = 0)/sum(cMF.perlen)
                        V1 = np.ma.masked_values(V1, cMF.hnoflo, atol = 0.09)
                        V += V1
                    del V1
                else:
                    MM = h5_MM['MM_S'][:,:,:,-1,index_S.get(i1)]
                    V = [np.sum(np.ma.masked_values(MM, cMF.hnoflo, atol = 0.09), axis = 0)/sum(cMF.perlen)]
                    V[0] = np.ma.masked_values(V[0], cMF.hnoflo, atol = 0.09)
                    i = 'Rp_botlayer'
                h5_MM.close()
                Vmax = np.ma.max(V) #float(np.ceil(np.ma.max(V)))  #
                Vmin = np.ma.min(V) #float(np.floor(np.ma.min(V))) #
                Vmin, Vmax, ctrs_tmp = minmax(Vmin, Vmax, ctrsMM)
                MMplot.plotLAYER(SP = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_'+ 'average_' + i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
                del V

                # ############################################
                # plot for selected time step
                # ############################################
                t = 0
                for SP in SP_lst:
                    h5_MM = h5py.File(h5_MM_fn, 'r')
                    if i1 != 'iRp':
                        V = [np.zeros([cMF.nrow,cMF.ncol])]
                        for l in range(_nslmax):
                            MM = h5_MM['MM_S'][SP,:,:,l,index_S.get(i1)]
                            V += np.ma.masked_values(MM, cMF.hnoflo, atol = 0.09)
                    else:
                        MM = h5_MM['MM_S'][SP,:,:,:,index_S.get(i1)]
                        V = [np.ma.masked_values(MM[:,:,-1], cMF.hnoflo, atol = 0.09)]
                    h5_MM.close()
                    MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_'+i), MM_ws = MM_ws_out, interval_type = 'linspace', interval_num = 5, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
                    t += 1
                del V, MM, t, Vmax, Vmin, ctrs_tmp
            del SP_lst, flxlbl, i, i1, h_diff_surf
        del gridSOIL, cMF.inputDate
        del hmaxMF, hminMF, hmin, hdiff, cbcmax_d, cbcmin_d

timeendExport = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
durationExport=(timeendExport-timestartExport)
durationTotal = (timeendExport-timestart)

##except StandardError:  #Exception
##    try:
##        h5_MM.close()
##    except:
##        pass
##    try:
##        h5_MF.close()
##    except:
##        pass
##    raise SystemExit('\nFATAL ERROR!\nMM run interruption in the export phase!\nError description:\n%s' % traceback.print_exc(file=sys.stdout))
###    traceback.print_exc(limit=1, file=sys.stdout)

# final report of successful run
print '\n##############\nMARMITES executed successfully!\n%s' % mpl.dates.datetime.datetime.today().isoformat()[:19]
if MMsoil_yn > 0:
    print '\nLOOP %d/%d' % (LOOP-1, ccnum)
    for txt in msg_end_loop:
        print txt
print '\n%d stress periods, %d days' % (cMF.nper,sum(cMF.perlen))
print '%d rows x %d cols (%d cells) x %d layers' % (cMF.nrow, cMF.ncol, cMF.nrow*cMF.ncol, cMF.nlay)
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
    print '\nMARMITES terminated!\n%s\n' % mpl.dates.datetime.datetime.today().isoformat()[:19]
    del s