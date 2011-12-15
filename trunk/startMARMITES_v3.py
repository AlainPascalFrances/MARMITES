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

""" See info in MARMITESunsat_v3.py"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.2"
__date__ = "November 2010"

import sys, os, traceback, h5py
import matplotlib as mpl
if mpl.get_backend!='agg':
    mpl.use('agg')
import matplotlib.pyplot as plt
plt.ioff()
import numpy as np
import startMARMITESsurface as startMMsurf
import MARMITESunsat_v3 as MMunsat
import MARMITESprocess_v3 as MMproc
import ppMODFLOW_flopy_v3 as ppMF
import MARMITESplot_v3 as MMplot
import CreateColors

#####################################

duration = 0.0
timestart = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
print '\n##############\nMARMITES started!\n%s\n##############' % mpl.dates.num2date(timestart).isoformat()[:19]

# workspace (ws) definition
# read input file (called _input.ini in the MARMITES workspace
# the first character on the first line has to be the character used to comment
# the file can contain any comments as the user wish, but the sequence of the input has to be respected
MM_ws = r'E:\00code_ws\CARRIZAL'  # 00_TESTS\MARMITESv3_r13c6l2'  SARDON'  CARRIZAL'
MM_fn = '__inputMM_g.ini'

inputFile = MMproc.readFile(MM_ws,MM_fn)

l=0
try:
    # # ECHO ON/OFF (1 - interpreter verbose BUT NO report, 0 - NO interpreter verbose BUT report)
    verbose = int(inputFile[l].strip())
    l += 1
    # output plot (1 is YES, 0 is NO)
    plot_out  = int(inputFile[l].strip())
    l += 1
    plot_freq =  int(inputFile[l].strip())
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
    #run MARMITESsurface  (1 is YES, 0 is NO)
    MMsurf_yn = int(inputFile[l].strip())
    l += 1
    #plot MARMITESsurface  (1 is YES, 0 is NO)
    MMsurf_plot = int(inputFile[l].strip())
    l += 1
    #run MARMITESunsat  (1 is YES, 0 is NO)
    MMunsat_yn = int(inputFile[l].strip())
    l += 1
    #run MODFLOW  (1 is YES, 0 is NO)
    MF_yn = int(inputFile[l].strip())
    l += 1
    # Define MARMITESsurface folder
    MMsurf_ws = inputFile[l].strip()
    l += 1
    # METEO TIME SERIES file name
    inputFile_TS_fn = inputFile[l].strip()
    l += 1
    # METEO/VEGETATION/SOIL/WATER PARAMETERS file name
    inputFile_PAR_fn = inputFile[l].strip()
    l += 1
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
    gridIRR_fn = inputFile[l].strip()
    l += 1
    gridSshmax_fn =  inputFile[l].strip()
    l += 1
    gridSsw_fn =  inputFile[l].strip()
    l += 1
    SOILparam_fn = inputFile[l].strip()
    l += 1
    IRR_fn = inputFile[l].strip()
    l += 1
    inputObs_fn = inputFile[l].strip()
    l += 1
    inputObsHEADS_fn = inputFile[l].strip()
    l += 1
    inputObsSM_fn = inputFile[l].strip()
    l += 1
    convcrit = float(inputFile[l].strip())
    l += 1
    ccnum = int(inputFile[l].strip())
    l += 1
    chunks = int(inputFile[l].strip())
    if MMunsat_yn == 1 and MF_yn != 1:
        MF_yn == 1
    if MMsurf_plot == 1:
        plot_out = 0
        plt_out_obs = 0
        print "\nYou required the MMsurf plots to appear on the screen. Due to backends limitation, MM and MF plots were disabled. Run again MM with MMsurf_plot = 0 to obtain the MM and MF plots."
except:
    raise SystemExit('\nType error in the input file %s' % (MM_fn))
del inputFile

if verbose == 0:
#capture interpreter output to be written in to a report file
    report_fn = os.path.join(MM_ws,'_MM_00report.txt')
    print '\nECHO OFF (no screen output).\nSee MM run report in file:\n%s\n' % report_fn
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
    outMMsurf_fn = startMMsurf.MMsurf(MMsurf_ws, inputFile_TS_fn, inputFile_PAR_fn, outputFILE_fn, MM_ws, outMMsurf_fn, MMsurf_plot)

inputFile = MMproc.readFile(MM_ws,outMMsurf_fn)

l=0
TRANS_vdw = []
Zr = []
kTu_min = []
kTu_n = []
TRANS_sdw = []
try:
    NMETEO = int(inputFile[l].strip())
    l += 1
    NVEG = int(inputFile[l].strip())
    l += 1
    NSOIL = int(inputFile[l].strip())
    l += 1
    inputDate_fn = str(inputFile[l].strip())
    l += 1
    inputZON_dTS_RF_fn = str(inputFile[l].strip())
    l += 1
    inputZON_dTS_PET_fn = str(inputFile[l].strip())
    l += 1
    inputZON_dTS_RFe_fn = str(inputFile[l].strip())
    l += 1
    inputZON_dTS_PE_fn = str(inputFile[l].strip())
    l += 1
    inputZON_dTS_E0_fn = str(inputFile[l].strip())
    l += 1
    line = inputFile[l].split()
    for v in range(NVEG):
        TRANS_vdw.append(int(line[v]))
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
    l += 1
    line = inputFile[l].split()
    for v in range(NSOIL):
        TRANS_sdw.append(int(line[v]))
except:
    raise SystemExit('\nType error in file [' + outMMsurf_fn + ']')
del inputFile

numDays = len(MMproc.readFile(MM_ws, inputDate_fn))

# #############################
# ###  READ MODFLOW CONFIG ####
# #############################

print'\n##############'
print 'MODFLOW initialization'
cMF = ppMF.MF(MM_ws, MF_ws, MF_ini_fn, xllcorner, yllcorner, numDays = numDays)

# ####   SUMMARY OF MODFLOW READINGS   ####
# active cells in layer 1_____________________ibound[0]
# elevation___________________________________top_l0
# aquifer type of layer 1_____________________laytyp[0]
# numero total de time step___________________sum(nstp)
# code for dry cell___________________________hdry
# cell size___________________________________delr[0]

# compute cbc conversion factor from volume to mm
if cMF.lenuni == 1:
    conv_fact = 304.8
elif cMF.lenuni == 2:
    conv_fact = 1000.0
elif cMF.lenuni == 3:
    conv_fact = 10.0
else:
    raise SystemExit('FATAL ERROR!\nDefine the length unit in the MODFLOW ini file!\n (see USGS Open-File Report 00-92)')
    # TODO if lenuni!=2 apply conversion factor to delr, delc, etc...
if cMF.laytyp[0]==0:
    raise SystemExit('FATAL ERROR!\nThe first layer cannot be confined type!\nChange your parameter laytyp in the MODFLOW lpf package.\n(see USGS Open-File Report 00-92)')
if cMF.itmuni != 4:
    raise SystemExit('FATAL ERROR!\nTime unit is not in days!')
iboundBOL = np.ones(np.array(cMF.ibound).shape, dtype = bool)
ncell = []
mask = []
for l in range(cMF.nlay):
    ncell.append((np.asarray(cMF.ibound)[:,:,l] != 0).sum())
    iboundBOL[:,:,l] = (np.asarray(cMF.ibound)[:,:,l] != 0)
    mask.append(np.ma.make_mask(iboundBOL[:,:,l]-1))

# #############################
# ### MF time processing
# #############################
# if required by user, compute nper, perlen,etc based on RF analysis in the METEO zones
if cMF.timedef >= 0:
    if isinstance(cMF.nper, str):
        try:
            perlenmax = int(cMF.nper.split()[1].strip())
        except:
            raise SystemExit('\nError in nper format of the MODFLOW ini file!\n')
    cMF.ppMFtime(MM_ws, MF_ws, inputDate_fn, inputZON_dTS_RF_fn, inputZON_dTS_PET_fn, inputZON_dTS_RFe_fn, inputZON_dTS_PE_fn, inputZON_dTS_E0_fn, NMETEO, NVEG, NSOIL)

print'\n##############'
print 'MARMITESunsat initialization'

# MARMITES INITIALIZATION
MM_UNSAT = MMunsat.UNSAT(hnoflo = cMF.hnoflo)
MM_SATFLOW = MMunsat.SATFLOW()

# READ input ESRI ASCII rasters # missing gridIRR_fn
print "\nImporting ESRI ASCII files to initialize MARMITES..."
gridMETEO = cMF.MM_PROCESS.inputEsriAscii(grid_fn = gridMETEO_fn, datatype = int)

gridSOIL = cMF.MM_PROCESS.inputEsriAscii(grid_fn = gridSOIL_fn, datatype = int)

gridSOILthick = cMF.MM_PROCESS.inputEsriAscii(grid_fn = gridSOILthick_fn,
 datatype = float)

gridSshmax = cMF.MM_PROCESS.inputEsriAscii(grid_fn = gridSshmax_fn,
 datatype = float)

gridSsw = cMF.MM_PROCESS.inputEsriAscii(grid_fn = gridSsw_fn,
 datatype = float)

##gridIRR = cMF.MM_PROCESS.inputEsriAscii(grid_fn                  = gridIRR_fn)

# READ input time series and parameters   # missing IRR_fn
gridVEGarea, RFzonesTS, E0zonesTS, PETvegzonesTS, RFevegzonesTS, PEsoilzonesTS, inputDate, JD = cMF.MM_PROCESS.inputTS(
                                NMETEO                   = NMETEO,
                                NVEG                     = NVEG,
                                NSOIL                    = NSOIL,
                                nstp                     = cMF.nstp,
                                inputDate_fn             = inputDate_fn,
                                inputZON_TS_RF_fn        = cMF.inputZON_TS_RF_fn,
                                inputZON_TS_PET_fn       = cMF.inputZON_TS_PET_fn,
                                inputZON_TS_RFe_fn       = cMF.inputZON_TS_RFe_fn,
                                inputZON_TS_PE_fn        = cMF.inputZON_TS_PE_fn,
                                inputZON_TS_E0_fn        = cMF.inputZON_TS_E0_fn
                                ) # IRR_fn

# SOIL PARAMETERS
_nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _Su_ini, _Ks = cMF.MM_PROCESS.inputSoilParam(MM_ws = MM_ws, SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
_nslmax = max(_nsl)

# compute thickness, top and bottom elevation of each soil layer
TopAquif = np.asarray(cMF.top) * 1000.0   # conversion from m to mm
# topography elevation
TopSoil = TopAquif + gridSOILthick*1000.0
del TopAquif

# create MM array
h5_MM_fn = os.path.join(MM_ws,'_h5_MM.h5')
# indexes of the HDF5 output arrays
index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iEXF':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSs':13, 'iETg':14, 'iETu':15, 'idSu':16, 'iHEADScorr':17, 'idtwt':18, 'iuzthick':19}
index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'iRexf':4, 'idSu':5, 'iSu':6, 'iSAT':7, 'iMB_l':8}

# READ observations time series (heads and soil moisture)
if plt_out_obs == 1:
    print "\nReading observations time series (hydraulic heads and soil moisture)..."
    obs, outpathname, obs_h, obs_S = cMF.MM_PROCESS.inputObs(
                                     MM_ws            = MM_ws,
                                     inputObs_fn      = inputObs_fn,
                                     inputObsHEADS_fn = inputObsHEADS_fn,
                                     inputObsSM_fn    = inputObsSM_fn,
                                     inputDate        = inputDate,
                                     _nslmax          = _nslmax,
                                     nlay             = cMF.nlay
                                     )
    obs['PzRCHmax'] = {'x':999,'y':999,'i': 999, 'j': 999, 'lay': 999, 'hi':999, 'h0':999, 'RC':999, 'STO':999}
    outpathname.append(os.path.join(MM_ws,'_MM_0PzRCHmax.txt'))
    obs_h.append([])
    obs_S.append([])
    # To write MM output in a txt file
    outFileExport = []
    for o in range(len(obs.keys())):
        outFileExport.append(open(outpathname[o], 'w'))
        Su_str   = ''
        Supc_str = ''
        dSu_str  = ''
        Rp_str   = ''
        Rexf_str = ''
        Eu_str   = ''
        Tu_str   = ''
        Smeasout = ''
        MB_str   = ''
        for l in range(_nslmax):
            Su_str = Su_str + 'Su_l' + str(l+1) + ','
            Supc_str = Supc_str + 'Supc_l' + str(l+1) + ','
            dSu_str = dSu_str + 'dSu_l' + str(l+1) + ','
            Eu_str = Eu_str + 'Eu_l' + str(l+1) + ','
            Tu_str = Tu_str + 'Tu_l' + str(l+1) + ','
            Rp_str = Rp_str + 'Rp_l' + str(l+1) + ','
            Rexf_str = Rexf_str + 'Rexf_l' + str(l+1) + ','
            MB_str = MB_str + 'MB_l' + str(l+1) + ','
            Smeasout = Smeasout + 'Smeas_' + str(l+1) + ','
        header='Date,RF,E0,PET,PE,RFe,I,' + Eu_str + Tu_str + 'Eg,Tg,ETg,WEL_MF,Es,' + Su_str + Supc_str + dSu_str + 'dSs,Ss,Ro,GW_EXF,' + Rp_str + Rexf_str + 'R_MF,hSATFLOW,hMF,hMFcorr,hmeas,dtwt,' + Smeasout + MB_str + 'MB\n'
        outFileExport[o].write(header)
    outPESTheads_fn      = 'PESTheads.dat'
    outPESTsm_fn         = 'PESTsm.dat'
    outPESTheads=open(os.path.join(MM_ws,outPESTheads_fn), 'w')
    outPESTsm=open(os.path.join(MM_ws,outPESTsm_fn), 'w')
    if cMF.uzf_yn == 1:
        cMF.uzf_obs(obs = obs)
else:
    print "\nNo reading of observations time series (hydraulic heads and soil moisture) required."
    obs = None

# #############################
# ### 1st MODFLOW RUN with initial user-input recharge
# #############################
if MF_yn == 1 :
    durationMF = 0.0
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
    durationMF +=  timeendMF-timestartMF

if isinstance(cMF.h5_MF_fn, str):
    try:
        h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
    except:
        raise SystemExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.')
    # heads format is : timestep, nrow, ncol, nlay
    # cbc format is: (kstp), kper, textprocess, nrow, ncol, nlay
    cbc_nam = []
    cbc_uzf_nam = []
    for c in h5_MF['cbc_nam']:
        cbc_nam.append(c.strip())
    if cMF.uzf_yn == 1:
        for c in h5_MF['cbc_uzf_nam']:
            cbc_uzf_nam.append(c.strip())
    elif cMF.rch_yn == 1:
        imfRCH = cbc_nam.index('RECHARGE')
        raise SystemExit('\nFATAL ERROR!\nMM has to be run together with the UZF1 package of MODFLOW 2005/ MODFLOW NWT, thus the    package should be desactivacted!\nExisting MM.')
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
    if MMunsat_yn == 0:
        h5_MF.close()

##except StandardError, e:  #Exception
##    raise SystemExit('\nFATAL ERROR!\nAbnormal MM run interruption in the initialization!\nError description:\n%s' % traceback.print_exc(file=sys.stdout))

# #############################
# 2nd phase : MM/MF loop #####
# #############################
h_diff_surf = None
#try:
if MMunsat_yn > 0:

    h_pSP = 0
    h_pSP_all = 0
    LOOP = 0
    LOOPlst = [LOOP]
    h_diff = [1000]
    h_diff_log = [1]
    h_diff_all = [1000]
    h_diff_all_log = [1]
    loopdry = 0
    plt_ConvLoop_fn = os.path.join(MM_ws, '_MM_00plt_MM_MF_ConvLoop.png')

    # #############################
    # ###  CONVERGENCE LOOP   #####
    # #############################

    while abs(h_diff[LOOP]) > convcrit:
        print '\n##############\nCONVERGENCE LOOP #%s\n##############' % str(LOOP)
        h_MF_average = 0.0
        timestartMMloop = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        # ###########################
        # ###  MARMITES INPUT #######
        # ###########################

        print'\n##############'
        print 'MARMITESunsat RUN'

        # SOIL PARAMETERS
        _nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _Su_ini, _Ks = cMF.MM_PROCESS.inputSoilParam(MM_ws = MM_ws, SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
        _nslmax = max(_nsl)

        # ###############
        # Create HDF5 arrays to store MARMITES output
        h5_MM = h5py.File(h5_MM_fn, 'w')
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
        # ###############
        # # main loop: calculation of soil water balance in each cell-grid for each time step inside each stress period
#        t0=0
        print '\nComputing...'

        # initial values of SP
        Su_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol,_nslmax])
        Rp_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol,_nslmax])
        Ss_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol])
        Ss_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol])
        for n in range(cMF.nper):
            tstart_MM = 0
            for t in range(n):
                tstart_MM += cMF.perlen[t]
            tend_MM = tstart_MM + cMF.perlen[n]
            tstart_MF = 0
            for t in range(n):
                tstart_MF += cMF.nstp[t]
            tend_MF = tstart_MF + cMF.nstp[n]
            MM = np.zeros([cMF.perlen[n],cMF.nrow,cMF.ncol,len(index)], dtype=float)
            MM_S = np.zeros([cMF.perlen[n],cMF.nrow,cMF.ncol,_nslmax,len(index_S)], dtype=float)
            MM_finf_MF = np.zeros([cMF.nrow,cMF.ncol], dtype=float)
            MM_wel_MF = np.zeros([cMF.nrow,cMF.ncol], dtype=float)
            h_MF_per = h5_MF['heads'][tstart_MF:tend_MF,:,:,0]
            exf_MF_per = h5_MF['cbc_uzf'][tstart_MF:tend_MF,imfEXF,:,:,0]
            # loop into the grid
            for i in range(cMF.nrow):
                for j in range(cMF.ncol):
                    SOILzone_tmp = gridSOIL[i,j]-1
                    METEOzone_tmp = gridMETEO[i,j]-1
                    if np.abs(iboundBOL[i,j,:]).sum() != 0.0:
                        slprop = _slprop[SOILzone_tmp]
                        nsl = _nsl[SOILzone_tmp]
                        # thickness of soil layers
                        Tl = gridSOILthick[i,j]*np.asarray(slprop)*1000.0
                        # elevation of top and bottom of soil layers
                        TopSoilLay = np.zeros([nsl], dtype=float)
                        BotSoilLay = np.zeros([nsl], dtype=float)
                        for l in range(nsl):
                            if l==0:
                                TopSoilLay[l] = TopSoil[i,j]
                                BotSoilLay[l] = TopSoil[i,j]-Tl[l]
                            else:
                                TopSoilLay[l] = BotSoilLay[l-1]
                                BotSoilLay[l] = TopSoilLay[l]-Tl[l]
                        Su_ini_tmp = []
                        if n == 0:
                            for l in range(nsl):
                                Su_ini_tmp.append(_Su_ini[SOILzone_tmp][l])
                            Ss_ini_tmp  = 0.0
                            dti = 1.0
                        else:
                            Su_ini_tmp    = Su_ini_tmp_array[i,j,:]
                            Ss_ini_tmp    = Ss_ini_tmp_array[i,j]
                        PEsoilzonesTS_tmp = PEsoilzonesTS[METEOzone_tmp,SOILzone_tmp,tstart_MF:tend_MF]
                        PEsoilzonesTS_tmp = np.asarray(PEsoilzonesTS_tmp)
                        PETvegzonesTS_tmp = []
                        RFevegzonesTS_tmp = []
                        for z in range(NVEG):
                            PETvegzonesTS_tmp.append(PETvegzonesTS[METEOzone_tmp,z,tstart_MF:tend_MF])
                            RFevegzonesTS_tmp.append(RFevegzonesTS[METEOzone_tmp,z,tstart_MF:tend_MF])
                        PETvegzonesTS_tmp = np.asarray(PETvegzonesTS_tmp)
                        RFevegzonesTS_tmp = np.asarray(RFevegzonesTS_tmp)
                        VEGarea_tmp=np.zeros([NVEG], dtype=np.float)
                        for v in range(NVEG):
                            VEGarea_tmp[v]=gridVEGarea[v,i,j]
                        h_MF_tmp = h_MF_per[:,i,j]
                        exf_MF_tmp = -exf_MF_per[:,i,j]
                        #conv_fact_tmp = -conv_fact/(delr[j]*delc[i])
                        # for the conv loop
                        # cal functions for reservoirs calculations
                        MM_tmp, MM_S_tmp = MM_UNSAT.run(
                                                     i, j, n,
                                                     nsl     = nsl,
                                                     st      = _st[SOILzone_tmp],
                                                     Sm      = _Sm[SOILzone_tmp],
                                                     Sfc     = _Sfc[SOILzone_tmp],
                                                     Sr      = _Sr[SOILzone_tmp],
                                                     Su_ini  = Su_ini_tmp,
                                                     Ss_ini  = Ss_ini_tmp,
                                                     Rp_ini  = Rp_ini_tmp_array[i,j,:],
                                                     botm_l0 = np.asarray(cMF.botm)[i,j,0],
                                                     TopSoilLay   = TopSoilLay,
                                                     BotSoilLay   = BotSoilLay,
                                                     Tl      = Tl,
                                                     Ks      = _Ks[SOILzone_tmp],
                                                     Ss_max  = 1000*1.12*gridSshmax[i,j]*gridSsw[i,j]/cMF.delr[j],
                                                     Ss_ratio= 1.12*gridSsw[i,j]/cMF.delr[j],
                                                     HEADS   = h_MF_tmp,
                                                     EXF     = exf_MF_tmp,
                                                     RF      = RFzonesTS[METEOzone_tmp][tstart_MF:tend_MF],
                                                     E0      = E0zonesTS[METEOzone_tmp][tstart_MF:tend_MF],
                                                     PETveg  = PETvegzonesTS_tmp,
                                                     RFeveg  = RFevegzonesTS_tmp,
                                                     PEsoil  = PEsoilzonesTS_tmp,
                                                     VEGarea = VEGarea_tmp,
                                                     Zr      = Zr,
                                                     nstp    = cMF.nstp[n],
                                                     perlen  = cMF.perlen[n],
                                                     dti     = dti,
                                                     hdry    = cMF.hdry,
                                                     kTu_min = kTu_min,
                                                     kTu_n   = kTu_n)
                        if (float(cMF.perlen[n])/float(cMF.nstp[n])) != 1.0:
                            for stp in range(cMF.nstp[n]):
                                ts = float(cMF.perlen[n])/float(cMF.nstp[n])
                                tstart =    stp*ts
                                tend   =   (1+stp)*ts
                                for k in range(len(index)):
                                    MM[tstart:tend,i,j,k] = MM_tmp[stp,k]
                                for k in range(len(index_S)):
                                    for l in range(nsl):
                                        MM_S[tstart:tend,i,j,l,k] = MM_S_tmp[stp,l,k]
                        else:
                            for k in range(len(index)):
                                MM[:,i,j,k] = MM_tmp[:,k]
                            for k in range(len(index_S)):
                                for l in range(nsl):
                                    MM_S[:,i,j,l,k] = MM_S_tmp[:,l,k]
                        MM_finf_MF[i,j] = MM_S_tmp[:,nsl-1,index_S.get('iRp')].sum()/conv_fact
                        MM_wel_MF[i,j] = MM_tmp[:,index.get('iETg')].sum()/conv_fact
                        del MM_tmp, MM_S_tmp
                        # setting initial conditions for the next SP
                        Su_ini_tmp_array[i,j,:]  = MM_S[cMF.nstp[n]-1,i,j,:,index_S.get('iSu')]
                        Rp_ini_tmp_array[i,j,:]  = MM_S[cMF.nstp[n]-1,i,j,:,index_S.get('iRp')]
                        Ss_ini_tmp_array[i,j]    = MM[cMF.nstp[n]-1,i,j,index.get('iSs')]
                    else:
                        if cMF.perlen[n]>1:
                            MM[:,i,j,:] = cMF.hnoflo
                            MM_S[:,i,j,:,:] = cMF.hnoflo
                        else:
                            MM[:,i,j,:] = cMF.hnoflo
                            MM_S[:,i,j,:,:] = cMF.hnoflo
                        MM_finf_MF[i,j] = cMF.hnoflo
                        MM_wel_MF[i,j] = cMF.hnoflo
            dti = float(cMF.perlen[n])/float(cMF.nstp[n])
            h5_MM['MM'][tstart_MM:tend_MM,:,:,:]     = MM[:,:,:,:]
            h5_MM['MM_S'][tstart_MM:tend_MM,:,:,:,:] = MM_S[:,:,:,:,:]
            h5_MM['finf'][n,:,:] = MM_finf_MF
            h5_MM['ETg'][n,:,:] = MM_wel_MF
        del MM, MM_S, MM_finf_MF, MM_wel_MF, exf_MF_per, Su_ini_tmp_array, Rp_ini_tmp_array, Ss_ini_tmp_array, dti
        h5_MM.close()

        # CHECK MM amd MF CONVERG.
        h_MF_per_m = np.ma.masked_values(np.ma.masked_values(h_MF_per, cMF.hdry, atol = 1E+25), cMF.hnoflo, atol = 0.09)
        del h_MF_per
        h_MF_average = np.ma.average(h_MF_per_m)
        h_diff.append(h_MF_average - h_pSP)
        h_diff_surf = h_MF_per_m - h_pSP_all
        h_diff_all_max = np.ma.max(h_diff_surf)
        h_diff_all_min = np.ma.min(h_diff_surf)
        if abs(h_diff_all_max)>abs(h_diff_all_min):
            h_diff_all.append(h_diff_all_max)
        else:
            h_diff_all.append(h_diff_all_min)
        LOOPlst.append(LOOP)
        LOOP += 1
        h_pSP = h_MF_average
        h_pSP_all = h_MF_per_m
        del h_MF_per_m
        if np.absolute(h_diff[LOOP])>0.0:
            h_diff_log.append(np.log10(np.absolute(h_diff[LOOP])))
            h_diff_all_log.append(np.log10(np.absolute(h_diff_all[LOOP])))
        else:
            h_diff_log.append(np.log10(convcrit))
            h_diff_all_log.append(np.log10(convcrit))

        timeendMMloop = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        durationMMloop = timeendMMloop-timestartMMloop
        print '\nMM loop run time: %02.fmn%02.fs' % (int(durationMMloop*24.0*60.0), (durationMMloop*24.0*60.0-int(durationMMloop*24.0*60.0))*60)

        if LOOP <2:
            print "\nInitial average heads:\n%.3f m" % h_diff[LOOP]
        else:
            print "\nHeads diff. from previous conv. loop: %.3f m" % h_diff[LOOP]
            print 'Maximum heads difference:             %.3f m' % h_diff_all[LOOP]
        if h_MF_average == 0.0:
            loopdry += 1
            if loopdry > 1:
                print '\nWARNING: first layer of the model DRY twice successively!\nLoop break, correct your MARMITES input value.'
                break
            else:
                print '\nWARNING: first layer of the model DRY!'
        elif abs(h_diff[LOOP]) < convcrit:
            print '\nSuccessfull convergence between MARMITES and MODFLOW!\n(Conv. criterion = %.3G)' % convcrit
            break
        elif LOOP>ccnum:
            print'\nNo convergence between MARMITES and MODFLOW!\n(Conv. criterion = %.3G)' % convcrit
            break
        del h_MF_average

        # MODFLOW RUN with MM-computed recharge
        h5_MF.close()
        durationMF = 0.0
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
            raise SystemExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.')

        timeendMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        durationMF += (timeendMF-timestartMF)

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
    del fig, LOOP, LOOPlst, h_diff, h_diff_log, h_pSP_all

    timeend = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    duration += (timeend-timestart) -durationMF

    # #############################
    # ### MODFLOW RUN with MM-computed recharge
    # #############################
    if MF_yn == 1 :
        h5_MF.close()
        durationMF = 0.0
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
        durationMF += (timeendMF-timestartMF)

##except StandardError, e:  #Exception
##    raise SystemExit('\nFATAL ERROR!\nAbnormal MM run interruption in the MM/MF loop!\nError description:\n%s' % traceback.print_exc(file=sys.stdout))

# #############################
# 3rd phase : export results #####
# #############################

#try:
del gridVEGarea
del RFzonesTS
del E0zonesTS
del PETvegzonesTS
del RFevegzonesTS
del PEsoilzonesTS
del gridMETEO
del gridSOILthick
del gridSshmax
del gridSsw
del TopSoil

# #############################
# ###  OUTPUT EXPORT   ########
# #############################

timestartExport = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
if plot_out == 1 or plt_out_obs == 1:
    print '\n##############\nMARMITES exporting...'
    # reorganizing MF output in daily data
    if MF_yn == 1 and isinstance(cMF.h5_MF_fn, str):
        print '\nConverting MODFLOW output into daily time step...'
        try:
            h5_MF = h5py.File(cMF.h5_MF_fn)
        except:
            raise SystemExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.')
        cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'STO_d', conv_fact = conv_fact, index = imfSTO)
        if cMF.drn_yn == 1:
            cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'DRN_d', conv_fact = conv_fact, index = imfDRN)
        if cMF.uzf_yn == 1:
            cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc_uzf', ds_name_new = 'RCH_d', conv_fact = conv_fact, index = imfRCH)
        if cMF.wel_yn == 1:
            cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'cbc', ds_name_new = 'WEL_d', conv_fact = conv_fact, index = imfWEL)
        cMF.MM_PROCESS.procMF(cMF = cMF, h5_MF = h5_MF, ds_name = 'heads', ds_name_new = 'heads_d', conv_fact = conv_fact)
        h5_MF.close()
    if isinstance(cMF.h5_MF_fn, str):
        top_m = np.ma.masked_values(cMF.top, cMF.hnoflo, atol = 0.09)
        index_cbc_uzf = [imfRCH]
        # TODO confirm the code below if wel_yn = 1 and drn_yn = 0
        if cMF.wel_yn == 1:
            if cMF.drn_yn == 1:
                index_cbc = [imfSTO, imfDRN, imfWEL]
            else:
                index_cbc = [imfSTO, imfWEL]
        else:
            if cMF.drn_yn == 1:
                index_cbc = [imfSTO, imfDRN]
            else:
                index_cbc = [imfSTO]
    else:
        cbc_DRN = cbc_STO = cbc_RCH = cbc_WEL = np.zeros((sum(cMF.perlen), cMF.nrow, cMF.ncol, cMF.nlay))
        imfDRN = imfSTO = imfRCH = imfWEL = 0
        top_m = np.zeros((cMF.nrow, cMF.ncol))

    print '\nExporting ASCII files and plots...'
    # computing max and min values in MF fluxes for plotting
    hmax = []
    hmin = []
    hmaxMM = []
    hminMM = []
    hdiff = []
    cbcmax = []
    cbcmin = []
    if isinstance(cMF.h5_MF_fn, str):
        # TODO missing STOuz (however this is not very relevant since these fluxes should not be the bigger in magnitude)
        try:
            h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
            # DRN
            if cMF.drn_yn == 1:
                cbc_DRN = h5_MF['DRN_d']
                DRNmax = np.ma.max(cbc_DRN)
                cbcmax.append(DRNmax)
                DRNmin = np.ma.min(cbc_DRN)
                del cbc_DRN
                cbcmin.append(DRNmin)
                DRNmax = -float(np.ceil(np.ma.max(DRNmin)))
                DRNmin = 0.0
            # STO
            cbc_STO = h5_MF['STO_d']
            cbcmax.append(np.ma.max(cbc_STO))
            cbcmin.append(np.ma.min(cbc_STO))
            del cbc_STO
            # RCH
            cbc_RCH = h5_MF['RCH_d']
            RCHmax = np.ma.max(cbc_RCH)
            cbcmax.append(RCHmax)
            RCHmin = np.ma.min(cbc_RCH)
            print '\nMaximum GW recharge (%.2f mm) observed at:' % RCHmax
            for l in range(cMF.nlay):
                for row in range(cMF.nrow):
                    for t,col in enumerate(cbc_RCH[:,row,:,l]):
                        try:
                            if plt_out_obs == 1:
                                obs['PzRCHmax'] = {'x':999,'y':999, 'i': row, 'j': list(col).index(RCHmax), 'lay': l, 'hi':999, 'h0':999, 'RC':999, 'STO':999}
                            print 'row %d, col %d and day %d (%s)' % (row, list(col).index(RCHmax), t, mpl.dates.num2date(inputDate[t]).isoformat()[:10])
                            tRCHmax = t
                        except:
                            pass
            del cbc_RCH
            RCHmax = float(np.ceil(np.ma.max(RCHmax)))
            RCHmin = float(np.floor(np.ma.min(RCHmin)))
            # WEL
            if cMF.wel_yn == 1:
                cbc_WEL = h5_MF['WEL_d']
                cbcmax.append(np.ma.max(cbc_WEL))
                cbcmin.append(np.ma.min(cbc_WEL))
            # GHB
            if cMF.ghb_yn == 1:
                cbc_GHB = h5_MF['GHB_d']
                cbcmax.append(np.ma.max(cbc_GHB))
                cbcmin.append(np.ma.min(cbc_GHB))
                del cbc_GHB
            cbcmax = float(np.ceil(np.ma.max(cbcmax)))
            cbcmin = float(np.floor(np.ma.min(cbcmin)))
            # h
            h_MF_m = np.ma.masked_values(np.ma.masked_values(h5_MF['heads_d'], cMF.hnoflo, atol = 0.09), cMF.hdry, atol = 1E+25)
            hmaxMF = float(np.ceil(np.ma.max(h_MF_m[:,:,:,:].flatten())))
            hminMF = float(np.floor(np.ma.min(h_MF_m[:,:,:,:].flatten())))
            h5_MF.close()
        except:
            raise SystemExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.')
    else:
        DRNmax = cbcmax = 1
        DRNmin = cbcmin = -1
        hmaxMF = -9999.9
        hminMF = 9999.9

    if obs != None:
        for o in range(len(obs.keys())):
            i = obs.get(obs.keys()[o])['i']
            j = obs.get(obs.keys()[o])['j']
            l = obs.get(obs.keys()[o])['lay']
            hmaxMF_tmp = float(np.ceil(np.ma.max(h_MF_m[:,i,j,l].flatten())))
            hminMF_tmp = float(np.floor(np.ma.min(h_MF_m[:,i,j,l].flatten())))
            if obs_h[o] != []:
                npa_m_tmp = np.ma.masked_values(obs_h[o], cMF.hnoflo, atol = 0.09)
                hmaxMM = float(np.ceil(np.ma.max(npa_m_tmp.flatten())))
                hminMM = float(np.floor(np.ma.min(npa_m_tmp.flatten())))
                del npa_m_tmp
            else:
                hmaxMM = -9999.9
                hminMM = 9999.9
            hmax.append(np.ma.max((hmaxMF_tmp, hmaxMM)))
            hmin.append(np.ma.min((hminMF_tmp, hminMM)))
            hdiff.append(hmax[o]-hmin[o])
        hdiff = np.ma.max(hdiff)
        del hmaxMF_tmp, hminMF_tmp, hmax
    else:
        hdiff = 2000

    # plot UNSAT/GW balance at the catchment scale
    if plt_out_obs == 1:
        if os.path.exists(h5_MM_fn):
            try:
                h5_MM = h5py.File(h5_MM_fn, 'r')
            except:
                raise SystemExit('\nFATAL ERROR!\nInvalid MARMITES HDF5 file. Run MARMITES and/or MODFLOW again.')
            # indexes of the HDF5 output arrays
            #    index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iEXF':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSs':13, 'iETg':14, 'iETu':15, 'idSu':16, 'iHEADScorr':17, 'idtwt':18}
            #   index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'iRexf':4, 'idSu':5, 'iSu':6, 'iSAT':7, 'iMB_l':8}
            flxlbl   = ['RF', 'I', 'EXF', 'dSs', 'Ro', 'Es', 'dSu']
            flxlbl1  = ['Eu', 'Tu']
            flxlbl2  = ['ETu', 'Eg']
            flxlbl3  = ['Tg']
            flxlbl3a = ['ETg']
            flxlbl4  = ['Rp']
            sign = [1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1]
            flxlst = []
            flxmax = []
            flxmin = []
            for i in flxlbl:
                i = 'i'+i
                flx_tmp = np.ma.masked_values(h5_MM['MM'][:,:,:,index.get(i)], cMF.hnoflo, atol = 0.09)
                flxmax.append(np.ma.max(flx_tmp))
                flxmin.append(np.ma.min(flx_tmp))
                flxlst.append(flx_tmp.sum()/sum(cMF.perlen)/ncell[0])
            for i in flxlbl1:
                tmp = 0.0
                flxlbl.append(i)
                i = 'i'+i
                for l in range(_nslmax):
                    flx_tmp = np.ma.masked_values(h5_MM['MM_S'][:,:,:,l,index_S.get(i)], cMF.hnoflo, atol = 0.09).sum()
                    flxmax.append(np.ma.max(flx_tmp))
                    flxmin.append(np.ma.min(flx_tmp))
                    tmp += flx_tmp
                flxlst.append(tmp/sum(cMF.perlen)/ncell[0])
                del tmp
            for i in flxlbl2:
                flxlbl.append(i)
                i = 'i'+i
                flx_tmp = np.ma.masked_values(h5_MM['MM'][:,:,:,index.get(i)], cMF.hnoflo, atol = 0.09)
                flxmax.append(np.ma.max(flx_tmp))
                flxmin.append(np.ma.min(flx_tmp))
                flxlst.append(flx_tmp.sum()/sum(cMF.perlen)/ncell[0])
            for i in flxlbl3:
                flxlbl.append(i)
                i = 'i'+i
                if cMF.wel_yn == 1:
                    flx_tmp = np.ma.masked_values(h5_MM['MM'][:,:,:,index.get(i)], cMF.hnoflo, atol = 0.09)
                    flxmax.append(np.ma.max(flx_tmp))
                    Tg_min = np.ma.min(flx_tmp)
                    flxmin.append(Tg_min)
                    flxlst.append(flx_tmp.sum()/sum(cMF.perlen)/ncell[0])
                    tTgmin = -1
                    if Tg_min < 0.0:
                        print '\nTg negative (%.2f) observed at:' % Tg_min
                        for row in range(cMF.nrow):
                            for t,col in enumerate(flx_tmp[:,row,:]):
                                try:
                                    print 'row %d, col %d and day %d' % (row, list(col).index(Tg_min), t)
                                    tTgmin = t
                                    if plt_out_obs == 1:
                                        obs['PzTgmin'] = {'x':999,'y':999, 'i': row, 'j': list(col).index(Tg_min), 'lay': 0, 'hi':999, 'h0':999, 'RC':999, 'STO':999}
                                        outpathname.append(os.path.join(MM_ws,'_MM_0PzTgmin.txt'))
                                        outFileExport.append(open(outpathname[len(outpathname)-1], 'w'))
                                        obs_h.append([])
                                        obs_S.append([])
                                        try:
                                            hmin.append(hmin[0])
                                        except:
                                            hmin.append(999.9)
                                        try:
                                            outFileExport[len(outFileExport)-1].write(header)
                                        except:
                                            pass
                                except:
                                    pass
                else:
                    flxlst.append(0.0)
            for i in flxlbl3a:
                flxlbl.append(i)
                i = 'i'+i
                flx_tmp = np.ma.masked_values(h5_MM['MM'][:,:,:,index.get(i)], cMF.hnoflo, atol = 0.09)
                flxmax.append(np.ma.max(flx_tmp))
                flxmin.append(np.ma.min(flx_tmp))
                flxlst.append(flx_tmp.sum()/sum(cMF.perlen)/ncell[0])
            for i in flxlbl4:
                flxlbl.append(i)
                i = 'i'+i
                flx_tmp = np.ma.masked_values(h5_MM['MM_S'][:,:,:,-1,index_S.get(i)], cMF.hnoflo, atol = 0.09)
                flxmax.append(np.ma.max(flx_tmp))
                flxmin.append(np.ma.min(flx_tmp))
                inf = flx_tmp.sum()/sum(cMF.perlen)/ncell[0]
                flxlst.append(inf)
            del flx_tmp
            h5_MM.close()
            for l,(x,y) in enumerate(zip(flxlst, sign)):
                flxlst[l] = x*y
            del flxlbl1, flxlbl2, flxlbl3, flxlbl3a, sign
            flxmax = float(np.ceil(np.ma.max(flxmax)))
            flxmin = float(np.floor(np.ma.min(flxmin)))
            #flxmax_avg = float(np.ceil(flxlst))
            #flxmin_avg = float(np.floor(flxlst))
            # TODO check where it is pertient to show flxmax, cbcmax, flxmax_avg... and same for min
            minfact = 0.95
            maxfact = 1.05
            if np.min((cbcmax, flxmax)) < 0:
                minfact = 1.05
            if np.max((cbcmax, flxmax)) < 0:
                maxfact = 0.95
            flxmax = maxfact*max(cbcmax, flxmax)
            flxmin = minfact*min(cbcmin, flxmin)
            if plt_out_obs == 1 and isinstance(cMF.h5_MF_fn, str):
                plt_export_fn = os.path.join(MM_ws, '_plt_0UNSATandGWbalances.png')
                # compute UZF_STO and store GW_RCH
                flxlbl.append('UZ_STO')
                rch_tmp = 0
                flxlst_tmp = []
                h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
                cbc_RCH = h5_MF['RCH_d']
                for l in range(cMF.nlay):
                    rch_tmp1 = cbc_RCH[:,:,:,l].sum()/sum(cMF.perlen)/ncell[l]
                    flxlst_tmp.append(rch_tmp1)
                    rch_tmp += rch_tmp1
                flxlst.append(-rch_tmp + inf)
                del rch_tmp, rch_tmp1, cbc_RCH
                for l in range(cMF.nlay):
                    # GW_RCH
                    flxlst.append(flxlst_tmp[l])
                    cbc_STO = h5_MF['STO_d']
                    flxlst.append(cbc_STO[:,:,:,l].sum()/sum(cMF.perlen)/ncell[l])
                    del cbc_STO
                    if cMF.drn_yn == 1:
                        cbc_DRN = h5_MF['DRN_d']
                        flxlst.append(cbc_DRN[:,:,:,l].sum()/sum(cMF.perlen)/ncell[l])
                        del cbc_DRN
                    if cMF.wel_yn == 1:
                        cbc_WEL = h5_MF['WEL_d']
                        flxlst.append(cbc_WEL[:,:,:,l].sum()/sum(cMF.perlen)/ncell[l])
                        del cbc_WEL
                    for x in range(len(index_cbc_uzf)):
                        flxlbl.append('GW_' + cbc_uzf_nam[index_cbc_uzf[x]] + '_L' + str(l+1))
                    for x in range(len(index_cbc)):
                        flxlbl.append('GW_' + cbc_nam[index_cbc[x]] + '_L' + str(l+1))
                    flxmax = float(np.ceil(np.ma.max(flxlst)))
                    flxmin = float(np.floor(np.ma.min(flxlst)))
                    flxmax = 1.15*max(cbcmax, flxmax)
                    flxmin = 1.15*min(cbcmin, flxmin)
                del flxlst_tmp
                h5_MF.close()
#                MB_tmp = sum(flxlst) - (flxlst[11] + flxlst[13] + flxlst[16] + flxlst[17] + flxlst[19] + flxlst[20])
                MB_tmp = cMF.hnoflo
                plt_title = 'MARMITES and MODFLOW water flux balance for the whole catchment\nsum of fluxes: %.4f' % MB_tmp
            else:
                plt_export_fn = os.path.join(MM_ws, '_plt_0UNSATbalance.png')
                #MB_tmp = sum(flxlst) - (flxlst[11] + flxlst[13])
                MB_tmp = cMF.hnoflo
                plt_title = 'MARMITES water flux balance for the whole catchment\nsum of fluxes: %.4f' % MB_tmp
            colors_flx = CreateColors.main(hi=0, hf=180, numbcolors = len(flxlbl))
            MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title, fluxmax = flxmax, fluxmin = flxmin)
            del flxlst
        # no MM, plot only MF balance if MF exists
        elif plt_out_obs == 1 and isinstance(cMF.h5_MF_fn, str):
            plt_export_fn = os.path.join(MM_ws, '_plt_GWbalances.png')
            flxlbl = []
            flxlst = []
            for l in range(cMF.nlay):
                h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
                cbc_RCH = h5_MF['RCH_d']
                flxlst.append(cbc_RCH[:,:,:,l].sum()/sum(cMF.perlen)/ncell[l])
                del cbc_RCH
                cbc_STO = h5_MF['STO_d']
                flxlst.append(cbc_STO[:,:,:,l].sum()/sum(cMF.perlen)/ncell[l])
                del cbc_STO
                if cMF.drn_yn == 1:
                    cbc_DRN = h5_MF['DRN_d']
                    flxlst.append(cbc_DRN[:,:,:,l].sum()/sum(cMF.perlen)/ncell[l])
                    del cbc_DRN
                if cMF.wel_yn == 1:
                    cbc_WEL = h5_MF['WEL_d']
                    flxlst.append(cbc_WEL[:,:,:,l].sum()/sum(cMF.perlen)/ncell[l])
                    del cbc_WEL
                h5_MF.close()
                for x in range(len(index_cbc_uzf)):
                    flxlbl.append('GW_' + cbc_uzf_nam[index_cbc_uzf[x]] + '_L' + str(l+1))
                for x in range(len(index_cbc)):
                    flxlbl.append('GW_' + cbc_nam[index_cbc[x]] + '_L' + str(l+1))
            colors_flx = CreateColors.main(hi=0, hf=180, numbcolors = len(flxlbl))
            #MB_tmp = sum(flxlst) - (flxlst[11] + flxlst[13] + flxlst[16] + flxlst[17] + flxlst[19] + flxlst[20])
            MB_tmp = cMF.hnoflo
            plt_title = 'MODFLOW water flux balance for the whole catchment\nsum of fluxes: %.4f' % MB_tmp
            MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title, fluxmax = cbcmax, fluxmin = cbcmin)
            del flxlst

    # exporting MM time series results to ASCII files and plots at observations cells
    if plt_out_obs == 1 and os.path.exists(h5_MM_fn):
        h5_MM = h5py.File(h5_MM_fn, 'r')
        h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
        colors_nsl = CreateColors.main(hi=00, hf=180, numbcolors = (_nslmax+1))
        for o in range(len(obs.keys())):
            i = obs.get(obs.keys()[o])['i']
            j = obs.get(obs.keys()[o])['j']
            l = obs.get(obs.keys()[o])['lay']
            SOILzone_tmp = gridSOIL[i,j]-1
            nsl = _nsl[SOILzone_tmp]
            MM = h5_MM['MM'][:,i,j,:]
            MM_S = h5_MM['MM_S'][:,i,j,:,:]
            # SATFLOW
            cbc_RCH = h5_MF['RCH_d']
            h_satflow = MM_SATFLOW.run(cbc_RCH[:,i,j,0], float(obs.get(obs.keys()[o])['hi']),float(obs.get(obs.keys()[o])['h0']),float(obs.get(obs.keys()[o])['RC']),float(obs.get(obs.keys()[o])['STO']))
            # export ASCII file at piezometers location
            #TODO extract heads at piezo location and not center of cell
            if obs_h[o] != []:
                obs_h_tmp = obs_h[o][0,:]
            else:
                obs_h_tmp = []
            if obs_S[o] != []:
                obs_S_tmp = obs_S[o]
            else:
                obs_S_tmp = []
            if cMF.wel_yn == 1:
                cbc_WEL = -h5_MF['WEL_d'][:,i,j,0]
            else:
                cbc_WEL = 0
            # Export time series results at observations points as ASCII file
            cMF.MM_PROCESS.ExportResultsMM(i, j, inputDate, _nslmax, MM, index, MM_S, index_S, cbc_RCH[:,i,j,0], cbc_WEL, h_satflow, h_MF_m[:,i,j,l], obs_h_tmp, obs_S_tmp, outFileExport[o], obs.keys()[o])
            del cbc_WEL
            outFileExport[o].close()
            # Export time series results at observations points as ASCII file for PEST
            # TODO reformulate the export format, it should be [date, SM_l1, SM_l2,...], i.e. the same format as the obs_SM and obs_heads files
            cMF.MM_PROCESS.ExportResultsPEST(i, j, inputDate, _nslmax, MM[:,index.get('iHEADScorr')], obs_h_tmp, obs_S_tmp, outPESTheads, outPESTsm, obs.keys()[o], MM_S[:,:,index_S.get('iSu_pc')])
            # plot time series results as plot
            if plt_out_obs == 1:
                plt_title = obs.keys()[o]
                # index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iEXF':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSs':13, 'iETg':14, 'iETu':15, 'idSu':16, 'iHEADScorr':17, 'idtwt':18}
                # index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'iRexf':4, 'idSu':5, 'iSu':6, 'iSAT':7, 'iMB_l':8}
                plt_export_fn = os.path.join(MM_ws, '_plt_0'+ obs.keys()[o] + '.png')
                # def plotTIMESERIES(DateInput, P, PET, PE, Pe, dPOND, POND, Ro, Eu, Tu, Eg, Tg, S, dS, Spc, Rp, EXF, ETg, Es, MB, MB_l, dtwt, SAT, R, h_MF, h_MF_corr, h_SF, hobs, Sobs, Sm, Sr, hnoflo, plt_export_fn, plt_title, colors_nsl, hmax, hmin):
                MMplot.plotTIMESERIES(
                inputDate,
                MM[:,index.get('iRF')],
                MM[:,index.get('iPET')],
                MM[:,index.get('iPE')],
                MM[:,index.get('iRFe')],
                MM[:,index.get('idSs')],
                MM[:,index.get('iSs')],
                MM[:,index.get('iRo')],
                MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iEu')],
                MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iTu')],
                MM[:,index.get('iEg')],
                MM[:,index.get('iTg')],
                MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iSu')],
                MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('idSu')],
                MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iSu_pc')],
                MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iRp')],
                MM[:,index.get('iEXF')],
                -MM[:,index.get('iETg')],
                MM[:,index.get('iEs')],
                MM[:,index.get('iMB')],
                MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iMB_l')],
                MM[:,index.get('idtwt')],
                MM[:,index.get('iuzthick')],
                MM_S[:,0:_nsl[gridSOIL[i,j]-1],index_S.get('iSAT')],
                cbc_RCH[:,i,j,0],
                h_MF_m[:,i,j,l], MM[:,index.get('iHEADScorr')], h_satflow, obs_h_tmp, obs_S_tmp,
                _Sm[gridSOIL[i,j]-1],
                _Sr[gridSOIL[i,j]-1],
                cMF.hnoflo,
                plt_export_fn,
                plt_title,
                colors_nsl,
                hmin[o] + hdiff,
                hmin[o],
                obs.keys()[o]
                )
                # plot water balance at each obs. cell
                plt_export_fn = os.path.join(MM_ws, '_plt_0'+ obs.keys()[o] + '_UNSATbalance.png')
                # flxlbl = ['RF', 'I', 'EXF', 'dSs', 'Ro', 'Es', 'dSu']
                # flxlbl1 = ['Eu', 'Tu']
                # flxlbl1a = ['Rp','Rexf']
                # flxlbl2 = ['ETu', 'Eg', 'Tg', 'ETg']
                flxlst =[MM[:,index.get('iRF')].sum()/sum(cMF.perlen),
                    -MM[:,index.get('iI')].sum()/sum(cMF.perlen),
                     MM[:,index.get('iEXF')].sum()/sum(cMF.perlen),
                     MM[:,index.get('idSs')].sum()/sum(cMF.perlen),
                    -MM[:,index.get('iRo')].sum()/sum(cMF.perlen),
                    -MM[:,index.get('iEs')].sum()/sum(cMF.perlen),
                     MM[:,index.get('idSu')].sum()/sum(cMF.perlen),
                    -MM_S[:,:,index_S.get('iEu')].sum()/sum(cMF.perlen),
                    -MM_S[:,:,index_S.get('iTu')].sum()/sum(cMF.perlen),
                    -MM[:,index.get('iETu')].sum()/sum(cMF.perlen),
                    -MM[:,index.get('iEg')].sum()/sum(cMF.perlen),
                    -MM[:,index.get('iTg')].sum()/sum(cMF.perlen),
                    -MM[:,index.get('iETg')].sum()/sum(cMF.perlen),
                    -MM_S[:,-1,index_S.get('iRp')].sum()/sum(cMF.perlen)]
                #MB_tmp = sum(flxlst) - (flxlst[11] + flxlst[13])
                MB_tmp = cMF.hnoflo
                if plt_out_obs == 1 and isinstance(cMF.h5_MF_fn, str):
                    plt_export_fn = os.path.join(MM_ws, '_plt_0'+ obs.keys()[o] + '_UNSATandGWbalances.png')
                    # compute UZF_STO and store GW_RCH
                    rch_tmp = 0
                    flxlst_tmp = []
                    for l in range(cMF.nlay):
                        rch_tmp1 = cbc_RCH[:,i,j,l].sum()/sum(cMF.perlen)
                        flxlst_tmp.append(rch_tmp1)
                        rch_tmp += rch_tmp1
                    flxlst.append(-rch_tmp + MM_S[:,-1,index_S.get('iRp')].sum()/sum(cMF.perlen))
                    del rch_tmp, rch_tmp1, cbc_RCH
                    for l in range(cMF.nlay):
                        # GW_RCH
                        flxlst.append(flxlst_tmp[l])
                        h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
                        cbc_STO = h5_MF['STO_d']
                        flxlst.append(cbc_STO[:,i,j,l].sum()/sum(cMF.perlen))
                        del cbc_STO
                        if cMF.drn_yn == 1:
                            cbc_DRN = h5_MF['DRN_d']
                            flxlst.append(cbc_DRN[:,i,j,l].sum()/sum(cMF.perlen))
                            del cbc_DRN
                        if cMF.wel_yn == 1:
                            cbc_WEL = h5_MF['WEL_d']
                            flxlst.append(cbc_WEL[:,i,j,l].sum()/sum(cMF.perlen))
                            del cbc_WEL
                    #MB_tmp -= flxlst[16] + flxlst[17] + flxlst[19] + flxlst[20]
                    MB_tmp = cMF.hnoflo
                MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title + '\nsum of fluxes: %.4f' % MB_tmp, fluxmax = flxmax, fluxmin = flxmin)
                del flxlst, flxlst_tmp
        del h_satflow, MM, MM_S
        h5_MM.close()
        h5_MF.close()
        # output for PEST
        outPESTheads.close()
        outPESTsm.close()
        del obs, obs_h, obs_S

    if plot_out == 1:
        TS_lst = []
        Date_lst = []
        JD_lst = []
        TS = 0
        while TS < len(h_MF_m):
            TS_lst.append(TS)
            Date_lst.append(inputDate[TS])
            JD_lst.append(JD[TS])
            TS += plot_freq
        if tTgmin < 0:
            lst = [len(h_MF_m)-1, tRCHmax]
        else:
            lst = [len(h_MF_m)-1, tRCHmax, tTgmin]
        for e in lst:
            TS_lst.append(e)
            Date_lst.append(inputDate[e])
            JD_lst.append(JD[e])

    # plot MF output
    if plot_out == 1 and isinstance(cMF.h5_MF_fn, str):
        # plot heads (grid + contours), DRN, etc... at specified TS
        h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
        # plot for selected time step
        t = 0
        for TS in TS_lst:
            # plot heads [m]
            V = []
            for L in range(cMF.nlay):
                V.append(h_MF_m[TS,:,:,L])
            if hmaxMF == hminMF:
                if hmaxMF == 0.0:
                    hmaxMF = 1.0
                    hminMF = -1.0
                else:
                    hmaxMF *= 1.15
                    hminMF *= 0.85
                ctrs_tmp = False
            else:
                ctrs_tmp = ctrsMF
            MMplot.plotLAYER(TS = TS, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation (m)', msg = 'DRY', plt_title = 'MF_HEADS', MM_ws = MM_ws, interval_type = 'arange', interval_diff = (hmaxMF - hminMF)/nrangeMF, contours = ctrs_tmp, Vmax = hmaxMF, Vmin = hminMF, ntick = ntick)
            # plot GW drainage [mm]
            if cMF.drn_yn == 1:
                V = []
                cbc_DRN = h5_MF['DRN_d']
                for L in range(cMF.nlay):
                    V.append(np.ma.masked_array(cbc_DRN[TS,:,:,L], mask[L])*(-1.0))
                DRNmax_tmp = np.ma.max(V)
                DRNmin_tmp = np.ma.min(V)
                if DRNmax_tmp == DRNmin_tmp:
                    if DRNmax_tmp == 0.0:
                        DRNmax_tmp = 1.0
                        DRNmin_tmp = -1.0
                    else:
                        DRNmax_tmp *= 1.15
                        DRNmin_tmp *= 0.85
                    ctrs_tmp = False
                else:
                    ctrs_tmp = ctrsMF
                MMplot.plotLAYER(TS = TS, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage (mm/day)', msg = '- no drainage', plt_title = 'MF_DRN', MM_ws = MM_ws, interval_type = 'linspace', interval_num = 5, Vmin = DRNmin_tmp, contours = ctrs_tmp, Vmax = DRNmax_tmp, fmt='%.3G', ntick = ntick)
            # plot GW RCH [mm]
            V = []
            cbc_RCH = h5_MF['RCH_d']
            for L in range(cMF.nlay):
                V.append(np.ma.masked_array(cbc_RCH[TS,:,:,L], mask[L]))
                RCHmax_tmp = np.ma.max(V)
                RCHmin_tmp = np.ma.min(V)
            if RCHmax_tmp == RCHmin_tmp:
                if RCHmax_tmp == 0.0:
                    RCHmax_tmp = 1.0
                    RCHmin_tmp = -1.0
                else:
                    RCHmax_tmp *= 1.15
                    RCHmin_tmp *= 0.85
                ctrs_tmp = False
            else:
                ctrs_tmp = ctrsMF
            MMplot.plotLAYER(TS = TS, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater recharge (mm/day)', msg = '- no flux', plt_title = 'MF_RCH', MM_ws = MM_ws, interval_type = 'linspace', interval_num = 5, Vmin = RCHmin_tmp, contours = ctrs_tmp, Vmax = RCHmax_tmp, fmt='%.3G', ntick = ntick)
            t += 1
            del V
        del t
        # plot average for the whole simulated period
        # plot GW drainage [mm]
        if cMF.drn_yn == 1:
            V = []
            for L in range(cMF.nlay):
                V.append(np.ma.masked_array(np.add.accumulate(cbc_DRN[:,:,:,L], axis = 0)[sum(cMF.perlen)-1,:,:]/sum(cMF.perlen)*(-1.0), mask[L]))
            DRNmax_tmp = np.ma.max(V)
            DRNmin_tmp = np.ma.min(V)
            if DRNmax_tmp == DRNmin_tmp:
                if DRNmax_tmp == 0.0:
                    DRNmax_tmp = 1.0
                    DRNmin_tmp = -1.0
                else:
                    DRNmax_tmp *= 1.15
                    DRNmin_tmp *= 0.85
                ctrs_tmp = False
            else:
                ctrs_tmp = ctrsMF
            MMplot.plotLAYER(TS = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage (mm/day)', msg = '- no drainage', plt_title = 'MF_average_DRN', MM_ws = MM_ws, interval_type = 'arange', interval_diff = (DRNmax_tmp - DRNmin_tmp)/nrangeMM, Vmax = DRNmax_tmp, Vmin = DRNmin_tmp, contours = ctrs_tmp, ntick = ntick)
        # plot GW RCH [mm]
        V = []
        for L in range(cMF.nlay):
            V.append(np.ma.masked_array(np.add.accumulate(cbc_RCH[:,:,:,L], axis = 0)[sum(cMF.perlen)-1,:,:]/sum(cMF.perlen), mask[L]))
        RCHmax_tmp = np.ma.max(V)
        RCHmin_tmp = np.ma.min(V)
        if RCHmax_tmp == RCHmin_tmp:
            if RCHmax_tmp == 0.0:
                RCHmax_tmp = 1.0
                RCHmin_tmp = -1.0
            else:
                RCHmax_tmp *= 1.15
                RCHmin_tmp *= 0.85
            ctrs_tmp = False
        else:
            ctrs_tmp = ctrsMF
        MMplot.plotLAYER(TS = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater recharge (mm/day)', msg = '- no flux', plt_title = 'MF_average_RCH', MM_ws = MM_ws, interval_type = 'arange', interval_diff = (RCHmax_tmp - RCHmin_tmp)/nrangeMM, Vmax = RCHmax_tmp, Vmin = RCHmin_tmp, contours = ctrs_tmp, ntick = ntick)
        h5_MF.close()
        del V, cbc_RCH
        if cMF.drn_yn == 1:
            del cbc_DRN, DRNmax_tmp, DRNmin_tmp
        del RCHmax_tmp, RCHmin_tmp

    # plot MM output
    if plot_out == 1 and os.path.exists(h5_MM_fn):
        flxlbl = ['RF', 'RFe', 'I', 'EXF', 'dSs', 'Ro', 'Es', 'Eg', 'Tg', 'ETg', 'ETu', 'dSu']
        for i in flxlbl:
            # plot average for the whole simulated period
            i1 = 'i'+i
            h5_MM = h5py.File(h5_MM_fn, 'r')
            MM = h5_MM['MM'][:,:,:,index.get(i1)]
            h5_MM.close()
            V = [np.add.accumulate(np.ma.masked_values(MM, cMF.hnoflo, atol = 0.09), axis = 0)[sum(cMF.perlen)-1,:,:]/sum(cMF.perlen)]
            V[0] = np.ma.masked_values(V[0], cMF.hnoflo, atol = 0.09)
            Vmax = np.ma.max(V[0]) #float(np.ceil(np.ma.max(V)))
            Vmin = np.ma.min(V[0]) #float(np.floor(np.ma.min(V)))
            if Vmax == Vmin:
                if Vmax == 0.0:
                    Vmax = 1.0
                    Vmin = -1.0
                else:
                    Vmax *= 1.15
                    Vmin *= 0.85
                ctrs_tmp = False
            else:
                ctrs_tmp = ctrsMM
            if i == 'dSu':
                i_lbl = '$\Delta$Su' #'$\Delta$S$_{u}'
            elif i == 'dSs':
                i_lbl = '$\Delta$Ss'
            else:
                i_lbl = i
            MMplot.plotLAYER(TS = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_average_' + i), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
            del V
            # plot for selected time step
            t = 0
            for TS in TS_lst:
                V = [np.ma.masked_values(MM[TS,:,:], cMF.hnoflo, atol = 0.09)]
                Vmax = np.ma.max(V[0]) #float(np.ceil(np.ma.max(V)))
                Vmin = np.ma.min(V[0]) #float(np.floor(np.ma.min(V)))
                if Vmax == Vmin:
                    if Vmax == 0.0:
                        Vmax = 1.0
                        Vmin = -1.0
                    else:
                        Vmax *= 1.15
                        Vmin *= 0.85
                    ctrs_tmp = False
                else:
                    ctrs_tmp = ctrsMM
                MMplot.plotLAYER(TS = TS, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_'+i), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
                t += 1
            del V, MM, t
        flxlbl = ['Eu', 'Tu','Rp']
        for i in flxlbl:
            # plot average for the whole simulated period
            i1 = 'i'+i
            h5_MM = h5py.File(h5_MM_fn, 'r')
            if i != 'Rp':
                V = [np.zeros([cMF.nrow,cMF.ncol])]
                for l in range(_nslmax):
                    MM = h5_MM['MM_S'][:,:,:,l,index_S.get(i1)]
                    V1 = np.add.accumulate(np.ma.masked_values(MM, cMF.hnoflo, atol = 0.09), axis = 0)[sum(cMF.perlen)-1,:,:]/sum(cMF.perlen)
                    V1 = np.ma.masked_values(V1, cMF.hnoflo, atol = 0.09)
                    V += V1
                del V1
            else:
                MM = h5_MM['MM_S'][:,:,:,-1,index_S.get(i1)]
                V = [np.add.accumulate(np.ma.masked_values(MM, cMF.hnoflo, atol = 0.09), axis = 0)[sum(cMF.perlen)-1,:,:]/sum(cMF.perlen)]
                V[0] = np.ma.masked_values(V[0], cMF.hnoflo, atol = 0.09)
                i = 'Rp_botlayer'
            h5_MM.close()
            Vmax = np.ma.max(V) #float(np.ceil(np.ma.max(V)))
            Vmin = np.ma.min(V) #float(np.floor(np.ma.min(V)))
            if Vmax == Vmin:
                if Vmax == 0.0:
                    Vmax = 1.0
                    Vmin = -1.0
                else:
                    Vmax *= 1.15
                    Vmin *= 0.85
                ctrs_tmp = False
            else:
                ctrs_tmp = ctrsMM
            MMplot.plotLAYER(TS = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_'+ 'average_' + i), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
            del V
            # plot for selected time step
            t = 0
            for TS in TS_lst:
                h5_MM = h5py.File(h5_MM_fn, 'r')
                if i1 != 'iRp':
                    V = [np.zeros([cMF.nrow,cMF.ncol])]
                    for l in range(_nslmax):
                        MM = h5_MM['MM_S'][TS,:,:,l,index_S.get(i1)]
                        V += np.ma.masked_values(MM, cMF.hnoflo, atol = 0.09)
                else:
                    MM = h5_MM['MM_S'][TS,:,:,:,index_S.get(i1)]
                    V = [np.ma.masked_values(MM[:,:,-1], cMF.hnoflo, atol = 0.09)]
                h5_MM.close()
                Vmax = np.ma.max(V) #float(np.ceil(np.ma.max(V)))
                Vmin = np.ma.min(V) #float(np.floor(np.ma.min(V)))
                if Vmax == Vmin:
                    if Vmax == 0.0:
                        Vmax = 1.0
                        Vmin = -1.0
                    else:
                        Vmax *= 1.15
                        Vmin *= 0.85
                    ctrs_tmp = False
                else:
                    ctrs_tmp = ctrsMM
                MMplot.plotLAYER(TS = TS, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_'+i), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
                t += 1
            del V, MM, t
        if h_diff_surf != None:
            flxlbl = ['h_diff_surf']
            Vmax = np.ma.max(h_diff_surf) #float(np.ceil(np.ma.max(V)))
            Vmin = np.ma.min(h_diff_surf) #float(np.floor(np.ma.min(V)))
            if Vmax == Vmin:
                if Vmax == 0.0:
                    Vmax = 1.0
                    Vmin = -1.0
                else:
                    Vmax *= 1.15
                    Vmin *= 0.85
                ctrs_tmp = False
            else:
                ctrs_tmp = ctrsMM
            MMplot.plotLAYER(TS = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = 1, nplot = 1, V = list(h_diff_surf),  cmap = plt.cm.Blues, CBlabel = ('(m)'), msg = 'no value', plt_title = ('_HEADSmaxdiff_ConvLoop'), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
        del TS_lst, flxlbl, i, i1, h_diff_surf

    del gridSOIL, inputDate
    del hmaxMF, hminMF, hmin, hdiff, cbcmax, cbcmin
    if cMF.drn_yn == 1:
        del DRNmax, DRNmin

timeendExport = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
durationExport=(timeendExport-timestartExport)

# final report of successful run
print '\n##############\nMARMITES executed successfully!\n%s\n' % mpl.dates.num2date(timestart).isoformat()[:19]
print '%d stress periods, %d days' % (cMF.nper,sum(cMF.perlen))
print '%d rows x %d cols (%d cells)' % (cMF.nrow,cMF.ncol,cMF.nrow*cMF.ncol)
l = 1
for n in ncell:
    print '%d active cells in layer %d' % (n, l)
    l += 1
print ('\nMARMITES run time: %s minute(s) and %.1f second(s)') % (str(int(duration*24.0*60.0)), (duration*24.0*60.0-int(duration*24.0*60.0))*60)
if MF_yn == 1:
    print ('MODFLOW run time: %s minute(s) and %.1f second(s)') % (str(int(durationMF*24.0*60.0)), (durationMF*24.0*60.0-int(durationMF*24.0*60.0))*60)
print ('Export run time: %s minute(s) and %.1f second(s)') % (str(int(durationExport*24.0*60.0)), (durationExport*24.0*60.0-int(durationExport*24.0*60.0))*60)
print ('\nOutput written in folder: \n%s\n##############\n') % MM_ws

##except StandardError, e:  #Exception
##    try:
##        h5_MM.close()
##        h5_MF.close()
##    except:
##        pass
##    raise SystemExit('\nFATAL ERROR!\nMM run interruption in the export phase!\nError description:\n%s' % traceback.print_exc(file=sys.stdout))
###    traceback.print_exc(limit=1, file=sys.stdout)

if verbose == 0:
    sys.stdout = s
    report.close()
    print '\nMARMITES terminated!\n%s\n' % mpl.dates.datetime.datetime.today().isoformat()[:19]
    del s
    try:
        h5_MM.close()
        h5_MF.close()
    except:
        pass