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

timestart = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
print '\n##############\nMARMITES started!\n%s\n##############' % mpl.dates.num2date(timestart).isoformat()[:19]

# workspace (ws) definition
# read input file (called _input.ini in the MARMITES workspace
# the first character on the first line has to be the character used to comment
# the file can contain any comments as the user wish, but the sequence of the input has to be respected
# 00_TESTS\MARMITESv3_r13c6l2'  00_TESTS\r40c20'  00_TESTS\r20c40'
# SARDON'  CARRIZAL' LAMATA'
MM_ws = r'E:\00code_ws\00_TESTS\MARMITESv3_r13c6l2'
MM_fn = '__inputMM.ini'

inputFile = MMproc.readFile(MM_ws,MM_fn)

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
    convcrit = float(inputFile[l].strip())
    l += 1
    ccnum = int(inputFile[l].strip())
    l += 1
    # Define MARMITESsurface folder
    MMsurf_ws = inputFile[l].strip()
    l += 1
    # METEO TIME SERIES file name
    inputFile_SP_fn = inputFile[l].strip()
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
    chunks = int(inputFile[l].strip())
    if MMunsat_yn == 1 and MF_yn != 1:
        MF_yn == 1
    if MMsurf_plot == 1:
        plt_out = 0
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
try:
# #############################
# ###  MARMITES surface  ######
# #############################
    print'\n##############'
    print 'MARMITESsurf RUN'
    if MMsurf_yn>0:
        durationMMsurf = 0.0
        timestartMMsurf = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        outMMsurf_fn = startMMsurf.MMsurf(MMsurf_ws, inputFile_SP_fn, inputFile_PAR_fn, outputFILE_fn, MM_ws, outMMsurf_fn, MMsurf_plot)
        timeendMMsurf = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        durationMMsurf=(timeendMMsurf-timestartMMsurf)
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
        inputZON_dSP_RF_fn = str(inputFile[l].strip())
        l += 1
        inputZON_dSP_PT_fn = str(inputFile[l].strip())
        l += 1
        inputZON_dSP_RFe_fn = str(inputFile[l].strip())
        l += 1
        inputZON_dSP_PE_fn = str(inputFile[l].strip())
        l += 1
        inputZON_dSP_E0_fn = str(inputFile[l].strip())
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
    durationMF = 0.0
    timestartMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    cMF = ppMF.MF(MM_ws, MF_ws, MF_ini_fn, xllcorner, yllcorner, numDays = numDays)
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
    ncell = []
    ncell_MM = []
    iboundBOL = np.ones(np.array(cMF.ibound).shape, dtype = bool)
    mask = []
    for l in range(cMF.nlay):
        ncell.append((np.asarray(cMF.ibound)[:,:,l] != 0).sum())
        ncell_MM.append((np.asarray(cMF.iuzfbnd) == l+1).sum())
        iboundBOL[:,:,l] = (np.asarray(cMF.ibound)[:,:,l] != 0)
        mask.append(np.ma.make_mask(iboundBOL[:,:,l]-1))
    del iboundBOL
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
                raise SystemExit('\nError in nper format of the MODFLOW ini file!\n')
        cMF.ppMFtime(MM_ws, MF_ws, inputDate_fn, inputZON_dSP_RF_fn, inputZON_dSP_PT_fn, inputZON_dSP_RFe_fn, inputZON_dSP_PE_fn, inputZON_dSP_E0_fn, NMETEO, NVEG, NSOIL)

    print'\n##############'
    print 'MARMITESunsat initialization'
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
    gridVEGarea, RFzonesSP, E0zonesSP, PTvegzonesSP, RFevegzonesSP, PEsoilzonesSP, inputDate, JD = cMF.MM_PROCESS.inputSP(
                                    NMETEO                   = NMETEO,
                                    NVEG                     = NVEG,
                                    NSOIL                    = NSOIL,
                                    nstp                     = cMF.nstp,
                                    inputDate_fn             = inputDate_fn,
                                    inputZON_SP_RF_fn        = cMF.inputZON_SP_RF_fn,
                                    inputZON_SP_PT_fn       = cMF.inputZON_SP_PT_fn,
                                    inputZON_SP_RFe_fn       = cMF.inputZON_SP_RFe_fn,
                                    inputZON_SP_PE_fn        = cMF.inputZON_SP_PE_fn,
                                    inputZON_SP_E0_fn        = cMF.inputZON_SP_E0_fn
                                    ) # IRR_fn

    # SOIL PARAMETERS
    _nsl, _nam_soil, _st, _facEg, _slprop, _Sm, _Sfc, _Sr, _Su_ini, _Ks = cMF.MM_PROCESS.inputSoilParam(MM_ws = MM_ws, SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
    _nslmax = max(_nsl)

    # compute thickness, top and bottom elevation of each soil layer
    TopAquif = np.asarray(cMF.top) * 1000.0   # conversion from m to mm
    botm_l0 = np.asarray(cMF.botm)[:,:,0]
    # topography elevation
    TopSoil = TopAquif + gridSOILthick*1000.0
    del TopAquif

    # create MM array
    h5_MM_fn = os.path.join(MM_ws,'_h5_MM.h5')
    # indexes of the HDF5 output arrays
    index = {'iRF':0, 'iPT':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iEXF':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSs':13, 'iETg':14, 'iETu':15, 'idSu':16, 'iinf':17, 'iHEADScorr':18, 'idtwt':19, 'iuzthick':20}
    index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'iRexf':4, 'idSu':5, 'iSu':6, 'iSAT':7, 'iMB_l':8}

    # READ observations time series (heads and soil moisture)
    if plt_out_obs == 1:
        print "\nReading observations time series (hydraulic heads and soil moisture)..."
        obs = cMF.MM_PROCESS.inputObs(MM_ws            = MM_ws,
                                      inputObs_fn      = inputObs_fn,
                                      inputObsHEADS_fn = inputObsHEADS_fn,
                                      inputObsSM_fn    = inputObsSM_fn,
                                      inputDate        = inputDate,
                                      _nslmax          = _nslmax,
                                      nlay             = cMF.nlay
                                      )
        # To write MM output in a txt file
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
        header='Date,RF,E0,PT,PE,RFe,I,' + Eu_str + Tu_str + 'Eg,Tg,ETg,WEL_MF,Es,' + Su_str + Supc_str + dSu_str + 'dSs,Ss,Ro,GW_EXF,' + Rp_str + Rexf_str + 'R_MF,hSATFLOW,hMF,hMFcorr,hmeas,dtwt,' + Smeasout + MB_str + 'MB\n'
        outPESTheads_fn      = 'PESTheads.dat'
        outPESTsm_fn         = 'PESTsm.dat'
        outPESTheads = open(os.path.join(MM_ws,outPESTheads_fn), 'w')
        outPESTsm = open(os.path.join(MM_ws,outPESTsm_fn), 'w')
        if cMF.uzf_yn == 1:
            cMF.uzf_obs(obs = obs)
    else:
        print "\nNo reading of observations time series (hydraulic heads and soil moisture) required."
        obs = None

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
        if MMunsat_yn == 0:
            h5_MF.close()

except StandardError, e:  #Exception
    raise SystemExit('\nFATAL ERROR!\nAbnormal MM run interruption in the initialization!\nError description:\n%s' % traceback.print_exc(file=sys.stdout))

# #############################
# 2nd phase : MM/MF loop #####
# #############################
h_diff_surf = None
try:
    if MMunsat_yn > 0:
        durationMMunsat = 0.0
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

        # #############################
        # ###  CONVERGENCE LOOP   #####
        # #############################

        while abs(h_diff[LOOP]) > convcrit:
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
            print 'MARMITESunsat RUN'
            # SOIL PARAMETERS
            _nsl, _nam_soil, _st, _facEg, _slprop, _Sm, _Sfc, _Sr, _Su_ini, _Ks = cMF.MM_PROCESS.inputSoilParam(MM_ws = MM_ws, SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
            _nslmax = max(_nsl)
            for l in range(NSOIL):
                _slprop[l] = np.asarray(_slprop[l])
            if LOOP > 0:
                h5_MM = h5py.File(h5_MM_fn)
            # ###############
            # # main loop: calculation of soil water balance in each cell-grid for each time step inside each stress period
    #        t0=0
            print '\nComputing...'
            # initial values of SP
            Su_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol,_nslmax])
            Rp_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol,_nslmax])
            Ss_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol])
            Ss_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol])
            MM_finf_MF = np.zeros([cMF.nrow,cMF.ncol], dtype=float)
            MM_wel_MF  = np.zeros([cMF.nrow,cMF.ncol], dtype=float)
            h_MF = None
            h_MF_mem = 'slow'
            try:
                h_MF = h5_MF['heads4MM'][:,:,:]
                h_MF_mem = 'fast'
            except:
                print '\nRAM memory too small compared to the size of the heads array -> slow computing.'
            if cMF.uzf_yn == 1:
                exf_MF = None
                exf_MF_mem = 'slow'
                try:
                    # TODO this below assume that there is no grid refinement
                    exf_MF = h5_MF['exf4MM'][:,:,:]*conv_fact/(cMF.delr[0]*cMF.delc[0])
                    exf_MF_mem = 'fast'
                except:
                    print '\nRAM memory too small compared to the size of the exfiltration array -> slow computing.'
            for n in range(cMF.nper):
                tstart_MM = 0
                for t in range(n):
                    tstart_MM += cMF.perlen[t]
                tend_MM = tstart_MM + cMF.perlen[n]
                tstart_MF = 0
                for t in range(n):
                    tstart_MF += cMF.nstp[t]
                tend_MF = tstart_MF + cMF.nstp[n]
                MM         = np.zeros([cMF.perlen[n],cMF.nrow,cMF.ncol,len(index)], dtype=float)
                MM_S       = np.zeros([cMF.perlen[n],cMF.nrow,cMF.ncol,_nslmax,len(index_S)], dtype=float)
                if h_MF == None:
                    h_MF       = h5_MF['heads4MM'][tstart_MF:tend_MF,:,:]
                if cMF.uzf_yn == 1:
                    if exf_MF == None:
                        exf_MF = h5_MF['exf4MM'][tstart_MF:tend_MF,:,:]*conv_fact/(cMF.delr[j]*cMF.delc[i])
                # loop into the grid
                for i in range(cMF.nrow):
                    for j in range(cMF.ncol):
                        SOILzone_tmp = gridSOIL[i,j]-1
                        METEOzone_tmp = gridMETEO[i,j]-1
                        if cMF.iuzfbnd[i][j] != 0.0:
                            _layer = cMF.iuzfbnd[i][j] - 1
                            facEg = _facEg[SOILzone_tmp]
                            slprop = _slprop[SOILzone_tmp]
                            nsl = _nsl[SOILzone_tmp]
                            # thickness of soil layers
                            Tl = gridSOILthick[i,j]*slprop*1000.0
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
                            PEsoilzonesSP_tmp = PEsoilzonesSP[METEOzone_tmp,SOILzone_tmp,tstart_MF:tend_MF]
                            PTvegzonesSP_tmp = np.zeros((NVEG,tend_MF-tstart_MF), dtype = np.float)
                            RFevegzonesSP_tmp = np.zeros((NVEG,tend_MF-tstart_MF), dtype = np.float)
                            for z in range(NVEG):
                                PTvegzonesSP_tmp[z,:] = PTvegzonesSP[METEOzone_tmp,z,tstart_MF:tend_MF]
                                RFevegzonesSP_tmp[z,:] = RFevegzonesSP[METEOzone_tmp,z,tstart_MF:tend_MF]
                            VEGarea_tmp=np.zeros([NVEG], dtype=np.float)
                            for v in range(NVEG):
                                VEGarea_tmp[v] = gridVEGarea[v,i,j]
                            if h_MF_mem == 'slow':
                                h_MF_tmp   = h_MF[:,i,j]
                            elif h_MF_mem == 'fast':
                                h_MF_tmp   = h_MF[tstart_MF:tend_MF,i,j]
                            if cMF.uzf_yn == 1:
                                if exf_MF_mem == 'slow':
                                    exf_MF_tmp = exf_MF[:,i,j]
                                elif exf_MF_mem == 'fast':
                                    exf_MF_tmp = exf_MF[tstart_MF:tend_MF,i,j]
                            else:
                                exf_MF_tmp = 0.0
                            MM_tmp, MM_S_tmp = MM_UNSAT.run(i          = i,
                                                            j          = j,
                                                            n          = n,
                                                            nsl        = nsl,
                                                            st         = _st[SOILzone_tmp],
                                                            Sm         = _Sm[SOILzone_tmp],
                                                            Sfc        = _Sfc[SOILzone_tmp],
                                                            Sr         = _Sr[SOILzone_tmp],
                                                            Su_ini     = Su_ini_tmp,
                                                            Ss_ini     = Ss_ini_tmp,
                                                            Rp_ini     = Rp_ini_tmp_array[i,j,:],
                                                            botm_l0    = botm_l0[i,j],
                                                            TopSoilLay = TopSoilLay,
                                                            BotSoilLay = BotSoilLay,
                                                            Tl         = Tl,
                                                            Ks         = _Ks[SOILzone_tmp],
                                                            Ss_max     = 1000*1.12*gridSshmax[i,j]*gridSsw[i,j]/cMF.delr[j],
                                                            Ss_ratio   = 1.12*gridSsw[i,j]/cMF.delr[j],
                                                            HEADS      = h_MF_tmp,
                                                            EXF        = -exf_MF_tmp,
                                                            RF         = RFzonesSP[METEOzone_tmp][tstart_MF:tend_MF],
                                                            E0         = E0zonesSP[METEOzone_tmp][tstart_MF:tend_MF],
                                                            PTveg      = PTvegzonesSP_tmp,
                                                            RFeveg     = RFevegzonesSP_tmp,
                                                            PEsoil     = PEsoilzonesSP_tmp,
                                                            VEGarea    = VEGarea_tmp,
                                                            Zr         = Zr,
                                                            nstp       = cMF.nstp[n],
                                                            perlen     = cMF.perlen[n],
                                                            dti        = dti,
                                                            hdry       = cMF.hdry,
                                                            kTu_min    = kTu_min,
                                                            kTu_n      = kTu_n,
                                                            facEg      = facEg)
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
                            MM_wel_MF[i,j] = MM_tmp[:,index.get('iETg')].sum()/(conv_fact)
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
                h5_MM['finf'][n,:,:]                     = MM_finf_MF
                h5_MM['ETg'][n,:,:]                      = MM_wel_MF
            del MM, MM_S, MM_finf_MF, MM_wel_MF, Su_ini_tmp_array, Rp_ini_tmp_array, Ss_ini_tmp_array, dti
            del h_MF, exf_MF, h_MF_tmp, exf_MF_tmp
            h5_MM.close()

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
                h_diff_all_log.append(np.log10(convcrit))

            timeendMMloop = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
            durationMMloop = timeendMMloop-timestartMMloop
            print '\nMM loop run time: %02.fmn%02.fs\n' % (int(durationMMloop*24.0*60.0), (durationMMloop*24.0*60.0-int(durationMMloop*24.0*60.0))*60)
            durationMMunsat += durationMMloop

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
                    break
                else:
                    print '\nWARNING: first layer of the model DRY!'
            elif abs(h_diff[LOOP]) < convcrit:
                msg_end_loop.append('Successfull convergence between MARMITES and MODFLOW!\n(Conv. criterion = %.3G)' % convcrit)
                for txt in msg_end_loop:
                    print txt
                break
            elif LOOP>ccnum:
                msg_end_loop.append('No convergence between MARMITES and MODFLOW!\n(Conv. criterion = %.3G)' % convcrit)
                for txt in msg_end_loop:
                    print txt
                break
            del h_MF_average
            for txt in msg_end_loop:
                print txt

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
        del fig, LOOPlst, h_diff, h_diff_log, h_pSP_all

        # #############################
        # ### MODFLOW RUN with MM-computed recharge
        # #############################
        if MF_yn == 1 :
            try:
                h5_MF.close()
            except:
                pass
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

except StandardError, e:  #Exception
    raise SystemExit('\nFATAL ERROR!\nAbnormal MM run interruption in the MM/MF loop!\nError description:\n%s' % traceback.print_exc(file=sys.stdout))

# #############################
# 3rd phase : export results #####
# #############################

##try:
del gridVEGarea
del RFzonesSP
del E0zonesSP
del PTvegzonesSP
del RFevegzonesSP
del PEsoilzonesSP
del gridMETEO
del gridSOILthick
del gridSshmax
del gridSsw
del TopSoil

# #############################
# ###  OUTPUT EXPORT   ########
# #############################

print '\n##############\nMARMITES exporting...'

timestartExport = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
# reorganizing MF output in daily data
if MF_yn == 1 and isinstance(cMF.h5_MF_fn, str):
    print '\nConverting MODFLOW output into daily time step...'
    try:
        h5_MF = h5py.File(cMF.h5_MF_fn)
    except:
        raise SystemExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.')
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
if isinstance(cMF.h5_MF_fn, str):
    top_m = np.ma.masked_values(cMF.top, cMF.hnoflo, atol = 0.09)
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
    top_m = np.zeros((cMF.nrow, cMF.ncol))

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
    Vmax = max(Vmax) #float(np.ceil(np.ma.max(V)))
    Vmin = min(Vmin) #float(np.floor(np.ma.min(V)))
    if Vmax == Vmin:
        if Vmax < 10E-9:
            Vmax =  1.0
            Vmin = -1.0
        else:
            Vmax *= 1.15
            Vmin *= 0.85
        ctrs_tmp = False
    else:
        ctrs_tmp = ctrsMF
    # TODO JD and Date are not correct since h_diff_n is # stress periods and not # of days (same in the plots of MF and MM)
    MMplot.plotLAYER(SP = h_diff_n, Date = inputDate[h_diff_n], JD = JD[h_diff_n], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = ('(m)'), msg = 'no value', plt_title = ('_HEADSmaxdiff_ConvLoop'), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMF, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
    del h_diff_n

if plt_out == 1 or plt_out_obs == 1:
    print '\nExporting ASCII files and plots...'
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
    if isinstance(cMF.h5_MF_fn, str):
        # TODO missing STOuz (however this is not very relevant since these fluxes should not be the bigger in magnitude)
        try:
            h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
            # DRN
            if cMF.drn_yn == 1:
                cbc_DRN = h5_MF['DRN_d']
                DRNmax = np.ma.max(cbc_DRN)
                cbcmax_d.append(DRNmax)
                DRNmin = np.ma.min(cbc_DRN)
                del cbc_DRN
                cbcmin_d.append(DRNmin)
                DRNmax = -float(np.ceil(np.ma.max(DRNmin)))
                DRNmin = 0.0
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
            if RCHmax> 0.0:
                for l in range(cMF.nlay):
                    for row in range(cMF.nrow):
                        for t,col in enumerate(cbc_RCH[:,row,:,l]):
                            try:
                                if plt_out_obs == 1:
                                    obs['PzRCHmax'] = {'x':999,'y':999, 'i': row, 'j': list(col).index(RCHmax), 'lay': l, 'hi':999, 'h0':999, 'RC':999, 'STO':999, 'outpathname':os.path.join(MM_ws,'_MM_0PzRCHmax.txt'), 'obs_h':[], 'obs_S':[]}
                                print 'row %d, col %d and day %d (%s)' % (row + 1, list(col).index(RCHmax) + 1, t, mpl.dates.num2date(inputDate[t] + 1.0).isoformat()[:10])
                                tRCHmax = t
                            except:
                                pass
            del cbc_RCH
            RCHmax = float(np.ceil(np.ma.max(RCHmax)))
            RCHmin = float(np.floor(np.ma.min(RCHmin)))
            # WEL
            if cMF.wel_yn == 1:
                cbc_WEL = h5_MF['WEL_d']
                cbcmax_d.append(np.ma.max(cbc_WEL))
                cbcmin_d.append(np.ma.min(cbc_WEL))
            # GHB
            if cMF.ghb_yn == 1:
                cbc_GHB = h5_MF['GHB_d']
                cbcmax_d.append(np.ma.max(cbc_GHB))
                cbcmin_d.append(np.ma.min(cbc_GHB))
                del cbc_GHB
            cbcmax_d = float(np.ceil(np.ma.max(cbcmax_d)))
            cbcmin_d = float(np.floor(np.ma.min(cbcmin_d)))
            # h
            h_MF_m = np.ma.masked_values(np.ma.masked_values(h5_MF['heads_d'], cMF.hnoflo, atol = 0.09), cMF.hdry, atol = 1E+25)
            hmaxMF = float(np.ceil(np.ma.max(h_MF_m[:,:,:,:].flatten())))
            hminMF = float(np.floor(np.ma.min(h_MF_m[:,:,:,:].flatten())))
            h5_MF.close()
        except:
            raise SystemExit('\nFATAL ERROR!\nInvalid MODFLOW HDF5 file. Run MARMITES and/or MODFLOW again.')
    else:
        DRNmax = cbcmax_d = 1
        DRNmin = cbcmin_d = -1
        hmaxMF = -9999.9
        hminMF = 9999.9

    if obs != None:
        x = 0
        for o in obs.keys():
            i = obs.get(o)['i']
            j = obs.get(o)['j']
            l = obs.get(o)['lay']
            obs_h = obs.get(o)['obs_h']
            hmaxMF_tmp = float(np.ceil(np.ma.max(h_MF_m[:,i,j,l].flatten())))
            hminMF_tmp = float(np.floor(np.ma.min(h_MF_m[:,i,j,l].flatten())))
            if obs_h != []:
                npa_m_tmp = np.ma.masked_values(obs_h, cMF.hnoflo, atol = 0.09)
                hmaxMM = float(np.ceil(np.ma.max(npa_m_tmp.flatten())))
                hminMM = float(np.floor(np.ma.min(npa_m_tmp.flatten())))
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

    # plot UNSAT/GW balance at the catchment scale
    tTgmin = -1
    if plt_out_obs == 1:
        if os.path.exists(h5_MM_fn):
            try:
                h5_MM = h5py.File(h5_MM_fn, 'r')
            except:
                raise SystemExit('\nFATAL ERROR!\nInvalid MARMITES HDF5 file. Run MARMITES and/or MODFLOW again.')
            # indexes of the HDF5 output arrays
            #    index = {'iRF':0, 'iPT':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iEXF':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSs':13, 'iETg':14, 'iETu':15, 'idSu':16, 'iHEADScorr':17, 'idtwt':18}
            #   index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'iRexf':4, 'idSu':5, 'iSu':6, 'iSAT':7, 'iMB_l':8}
            flxlbl   = ['RF', 'I', 'RFe', 'dSs', 'Ro', 'Es', 'dSu', 'EXF']
            flxlbl1  = ['Eu', 'Tu']
            flxlbl2  = ['ETu', 'Eg']
            flxlbl3  = ['Tg']
            flxlbl3a = ['ETg']
            flxlbl4  = ['Rp']
            sign =     [   1,  -1,     1,     1,   -1,   -1,    1,     1, -1, -1, -1, -1, -1, -1, -1]
            flxlst = []
            flx_Cat_TS = []
            flxmax_d = []
            flxmin_d = []
            for i in flxlbl:
                i = 'i'+i
                array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
                flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
                flxmax_d.append(np.ma.max(flx_tmp))
                flxmin_d.append(np.ma.min(flx_tmp))
                flxlst.append(facTim*(flx_tmp.sum())/sum(cMF.perlen)/sum(ncell_MM))
                array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
##                if i != 'iRo':
                flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
##                else:
##                flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1))
                del flx_tmp, array_tmp, array_tmp1
            for i in flxlbl1:
                flxlbl.append(i)
                i = 'i'+i
                flx_tmp1 = 0.0
                for l in range(_nslmax):
                    array_tmp = h5_MM['MM_S'][:,:,:,l,index_S.get(i)]
                    flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
                    flxmax_d.append(np.ma.max(flx_tmp))
                    flxmin_d.append(np.ma.min(flx_tmp))
                    flx_tmp1 += flx_tmp.sum()
                    array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
                    array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
                flxlst.append(facTim*flx_tmp1/sum(cMF.perlen)/sum(ncell_MM))
                flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
                del flx_tmp, flx_tmp1, array_tmp, array_tmp1
            for i in flxlbl2:
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
                flxlbl.append(i)
                i = 'i'+i
                if cMF.wel_yn == 1:
                    array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
                    flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
                    flxmax_d.append(np.ma.max(flx_tmp))
                    Tg_min = np.ma.min(flx_tmp)
                    flxmin_d.append(Tg_min)
                    flxlst.append(facTim*(flx_tmp.sum())/sum(cMF.perlen)/sum(ncell_MM))
                    if Tg_min < 0.0:
                        print '\nTg negative (%.2f) observed at:' % Tg_min
                        for row in range(cMF.nrow):
                            for t,col in enumerate(flx_tmp[:,row,:]):
                                try:
                                    print 'row %d, col %d and day %d' % (row + 1, list(col).index(Tg_min) + 1, t + 1)
                                    tTgmin = t
                                    if plt_out_obs == 1:
                                        obs['PzTgmin'] = {'x':999,'y':999, 'i': row, 'j': list(col).index(Tg_min), 'lay': 0, 'hi':999, 'h0':999, 'RC':999, 'STO':999, 'outpathname':os.path.join(MM_ws,'_MM_0PzTgmin.txt'), 'obs_h':[], 'obs_S':[]}
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
                flxlbl.append(i)
                i = 'i'+i
                array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
                flx_tmp = np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09)
                flxmax_d.append(np.ma.max(flx_tmp))
                flxmin_d.append(np.ma.min(flx_tmp))
                flxlst.append(facTim*flx_tmp.sum()/sum(cMF.perlen)/sum(ncell_MM))
                array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
                flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
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
            for i in ['Ss', 'PE', 'PT', 'inf']:
                i = 'i'+i
                array_tmp = h5_MM['MM'][:,:,:,index.get(i)]
                array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
                flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
                del array_tmp, array_tmp1
            plt_exportCATCH_fn = os.path.join(MM_ws, '_plt_0CATCHMENT.png')
            plt_titleCATCH = 'Time serie of fluxes averaged over the whole catchment'
            InMM = flxlst[0] + flxlst[3] + flxlst[6] + flxlst[7]
            OutMM = flxlst[1] + flxlst[4] + flxlst[5] + flxlst[10] + flxlst[14]
            MB_MM = 100*(InMM - OutMM)/((InMM+OutMM)/2)
            h5_MM.close()
            flxmax_d = float(np.ceil(np.ma.max(flxmax_d)))
            flxmin_d = float(np.floor(np.ma.min(flxmin_d)))
            flxmax = axefact*float(np.ceil(max(flxlst)))
            flxmin = axefact*float(np.floor(min(flxlst)))
            for l,(x,y) in enumerate(zip(flxlst, sign)):
                flxlst[l] = x*y
            del flxlbl1, flxlbl2, flxlbl3, flxlbl3a, sign
            if plt_out_obs == 1 and isinstance(cMF.h5_MF_fn, str):
                plt_export_fn = os.path.join(MM_ws, '_plt_0CATCHMENT_UNSATandGWbalances.png')
                # compute UZF_STO and store GW_RCH
                flxlbl.append('UZ_STO')
                rch_tmp = 0
                flxlst_tmp = []
                h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
                cbc_RCH = h5_MF['RCH_d']
                for l in range(cMF.nlay):
                    array_tmp = cbc_RCH[:,:,:,l]
                    array_tmp1 = np.sum(np.ma.masked_values(array_tmp, cMF.hnoflo, atol = 0.09), axis = 1)
                    flx_Cat_TS.append(np.sum(np.ma.masked_values(array_tmp1, cMF.hnoflo, atol = 0.09), axis = 1)/sum(ncell_MM))
                    rch_tmp1 = facTim*(array_tmp.sum())/sum(cMF.perlen)/sum(ncell_MM)
                    flxlst_tmp.append(rch_tmp1)
                    rch_tmp += rch_tmp1
                del array_tmp, array_tmp1, rch_tmp1, cbc_RCH
                flxlst.append(rch_tmp - inf)
                del rch_tmp, inf
                MMplot.plotTIMESERIES_CATCH(inputDate, flx_Cat_TS, plt_exportCATCH_fn, plt_titleCATCH, MF = 'y')
                InUZF = -flxlst[14] + flxlst[15]
                OutUZF = 0
                InMF = 0
                OutMF = 0
                for l in range(cMF.nlay):
                    # GW_STO
                    cbc_STO = h5_MF['STO_d']
                    flxlst.append(facTim*(cbc_STO[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))  # -1*
                    InMF += flxlst[-1]
                    del cbc_STO
                    # GW_RCH
                    flxlst.append(flxlst_tmp[l])
                    OutUZF += flxlst_tmp[l]
                    InMF += flxlst_tmp[l]
                    # EXF
                    cbc_EXF = h5_MF['EXF_d']
                    flxlst.append(facTim*(cbc_EXF[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))
                    OutMF += -flxlst[-1]
                    del cbc_EXF
                    if cMF.drn_yn == 1:
                        cbc_DRN = h5_MF['DRN_d']
                        if cMF.drncells[l]>0:
                            flxlst.append(facTim*(cbc_DRN[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))
                            OutMF += -flxlst[-1]
                        else:
                            flxlst.append(0.0)
                        del cbc_DRN
                    if cMF.wel_yn == 1:
                        cbc_WEL = h5_MF['WEL_d']
                        if ncell_MM[l]>0:
                            flxlst.append(facTim*(cbc_WEL[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))
                            OutMF += -flxlst[-1]
                        else:
                            flxlst.append(0.0)
                        del cbc_WEL
                    if cMF.ghb_yn == 1:
                        cbc_GHB = h5_MF['GHB_d']
                        if cMF.ghbcells[l] > 0:
                            flxlst.append(facTim*(cbc_GHB[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))
                            OutMF += -flxlst[-1]
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
                plt_title = 'MARMITES and MODFLOW water flux balance for the whole catchment\nMass balance error: MM = %1.2f%%, UZF = %1.2f%%, MF = %1.2f%%' % (MB_MM, MB_UZF, MB_MF)
            else:
                MMplot.plotTIMESERIES_CATCH(inputDate, flx_Cat_TS, plt_exportCATCH_fn, plt_titleCATCH, MF = 'n')
                plt_export_fn = os.path.join(MM_ws, '_plt_0CATCHMENT_UNSATbalance.png')
                plt_title = 'MARMITES water flux balance for the whole catchment\nMass balance error: MM = %1.2f%%' % (MB_MM)
            colors_flx = CreateColors.main(hi=0, hf=180, numbcolors = len(flxlbl))
            MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title, fluxmax = flxmax, fluxmin = flxmin, unit = plt_WB_unit)
            del flxlst
        # no MM, plot only MF balance if MF exists
        elif plt_out_obs == 1 and isinstance(cMF.h5_MF_fn, str):
            plt_export_fn = os.path.join(MM_ws, '_plt_0CATCHMENT_GWbalances.png')
            flxlbl = []
            flxlst = []
            InMF = 0
            OutMF = 0
            for l in range(cMF.nlay):
                h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
                # GW STO
                cbc_STO = h5_MF['STO_d']
                flxlst.append(facTim*(cbc_STO[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))  # *-1
                InMF += flxlst[-1]
                del cbc_STO
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
                        flxlst.append(facTim*(cbc_WEL[:,:,:,l].sum()/sum(cMF.perlen)/sum(ncell_MM)))
                        OutMF += flxlst[-1]
                    else:
                        flxlst.append(0.0)
                    del cbc_WEL
                if cMF.ghb_yn == 1:
                    cbc_GHB = h5_MF['GHB_d']
                    if cMF.ghbcells > 0:
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
            plt_title = 'MODFLOW water flux balance for the whole catchment\nMass balance error: MF = %1.2f%%' % (MB_MF)
            MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title, fluxmax = cbcmax_d, fluxmin = cbcmin_d, unit = plt_WB_unit)
            del flxlst

    # exporting MM time series results to ASCII files and plots at observations cells
    if plt_out_obs == 1 and os.path.exists(h5_MM_fn):
        h5_MM = h5py.File(h5_MM_fn, 'r')
        h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
        colors_nsl = CreateColors.main(hi=00, hf=180, numbcolors = (_nslmax+1))
        x = 0
        flxlst = []
        plt_exportBAL_fn = []
        plt_titleBAL = []
        for o in obs.keys():
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
            if obs_S != []:
                obs_S_tmp = obs_S
            else:
                obs_S_tmp = []
            if cMF.wel_yn == 1:
                cbc_WEL = -h5_MF['WEL_d'][:,i,j,0]
            else:
                cbc_WEL = 0
            # Export time series results at observations points as ASCII file
            cMF.MM_PROCESS.ExportResultsMM(i, j, inputDate, _nslmax, MM, index, MM_S, index_S, cbc_RCH[:,i,j,0], cbc_WEL, h_satflow, h_MF_m[:,i,j,l], obs_h_tmp, obs_S_tmp, outFileExport, o)
            del cbc_WEL
            outFileExport.close()
            # Export time series results at observations points as ASCII file for PEST
            # TODO reformulate the export format, it should be [date, SM_l1, SM_l2,...], i.e. the same format as the obs_SM and obs_heads files
            cMF.MM_PROCESS.ExportResultsPEST(i, j, inputDate, _nslmax, MM[:,index.get('iHEADScorr')], obs_h_tmp, obs_S_tmp, outPESTheads, outPESTsm, o, MM_S[:,:,index_S.get('iSu_pc')])
            # plot time series results as plot
            if plt_out_obs == 1:
                plt_title = o + '\ni = %d, j = %d, l = %d, x = %d, y = %d, %s' % (i+1, j+1, l+1, obs.get(o)['x'], obs.get(o)['y'], soilnam)
                plt_titleBAL.append(plt_title)
                # index = {'iRF':0, 'iPT':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iEXF':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSs':13, 'iETg':14, 'iETu':15, 'idSu':16, 'iHEADScorr':17, 'idtwt':18}
                # index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'iRexf':4, 'idSu':5, 'iSu':6, 'iSAT':7, 'iMB_l':8}
                plt_export_fn = os.path.join(MM_ws, '_plt_0'+ o + '.png')
                # def plotTIMESERIES(DateInput, P, PT, PE, Pe, dPOND, POND, Ro, Eu, Tu, Eg, Tg, S, dS, Spc, Rp, EXF, ETg, Es, MB, MB_l, dtwt, SAT, R, h_MF, h_MF_corr, h_SF, hobs, Sobs, Sm, Sr, hnoflo, plt_export_fn, plt_title, colors_nsl, hmax, hmin):
                MMplot.plotTIMESERIES(
                inputDate,
                MM[:,index.get('iRF')],
                MM[:,index.get('iPT')],
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
                hmax[x] + hdiff/2,
                hmin[x] - hdiff/2,
                o
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
                    -1*facTim*(MM[:,index.get('iI')].sum()/sum(cMF.perlen)),
                     facTim*(MM[:,index.get('iRFe')].sum()/sum(cMF.perlen)),
                     facTim*(MM[:,index.get('idSs')].sum()/sum(cMF.perlen)),
                     -1*facTim*(MM[:,index.get('iRo')].sum()/sum(cMF.perlen)),
                    -1*facTim*(MM[:,index.get('iEs')].sum()/sum(cMF.perlen)),
                     facTim*(MM[:,index.get('idSu')].sum()/sum(cMF.perlen)),
                     facTim*(MM[:,index.get('iEXF')].sum()/sum(cMF.perlen)),
                    -1*facTim*(MM_S[:,:,index_S.get('iEu')].sum()/sum(cMF.perlen)),
                    -1*facTim*(MM_S[:,:,index_S.get('iTu')].sum()/sum(cMF.perlen)),
                    -1*facTim*(MM[:,index.get('iETu')].sum()/sum(cMF.perlen)),
                    -1*facTim*(MM[:,index.get('iEg')].sum()/sum(cMF.perlen)),
                    -1*facTim*(MM[:,index.get('iTg')].sum()/sum(cMF.perlen)),
                    -1*facTim*(MM[:,index.get('iETg')].sum()/sum(cMF.perlen)),
                    -1*facTim*conv_fact*((cMF.perlen*h5_MM['finf'][:,i,j]).sum()/sum(cMF.perlen))
                    ])
                if plt_out_obs == 1 and isinstance(cMF.h5_MF_fn, str):
                    plt_exportBAL_fn.append(os.path.join(MM_ws, '_plt_0'+ o + '_UNSATandGWbalances.png'))
                    # compute UZF_STO and store GW_RCH
                    rch_tmp = 0
                    flxlst_tmp = []
                    for l in range(cMF.nlay):
                        rch_tmp1 = facTim*(cbc_RCH[:,i,j,l].sum()/sum(cMF.perlen))
                        flxlst_tmp.append(rch_tmp1)
                        rch_tmp += rch_tmp1
                    flxlst[-1].append(-rch_tmp + facTim*conv_fact*((cMF.perlen*h5_MM['finf'][:,i,j]).sum()/sum(cMF.perlen)))
                    del rch_tmp, rch_tmp1, cbc_RCH
                    for l in range(cMF.nlay):
                        # GW STO
                        cbc_STO = h5_MF['STO_d']
                        flxlst[-1].append(facTim*(cbc_STO[:,i,j,l].sum()/sum(cMF.perlen)))    # -1*
                        del cbc_STO
                        # GW_RCH
                        flxlst[-1].append(flxlst_tmp[l])
                        h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
                        # EXF
                        cbc_EXF = h5_MF['EXF_d']
                        flxlst[-1].append(facTim*(cbc_EXF[:,i,j,l].sum()/sum(cMF.perlen)))
                        del cbc_EXF
                        if cMF.drn_yn == 1:
                            cbc_DRN = h5_MF['DRN_d']
                            flxlst[-1].append(facTim*(cbc_DRN[:,i,j,l].sum()/sum(cMF.perlen)))
                            del cbc_DRN
                        if cMF.wel_yn == 1:
                            cbc_WEL = h5_MF['WEL_d']
                            flxlst[-1].append(facTim*(cbc_WEL[:,i,j,l].sum()/sum(cMF.perlen)))
                            del cbc_WEL
                        if cMF.ghb_yn == 1:
                            cbc_GHB = h5_MF['GHB_d']
                            flxlst[-1].append(facTim*(cbc_GHB[:,i,j,l].sum()/sum(cMF.perlen)))
                            del cbc_GHB
                    del flxlst_tmp
                else:
                    plt_exportBAL_fn.append(os.path.join(MM_ws, '_plt_0'+ o + '_UNSATbalance.png'))
        flxmax = axefact*float(np.ceil(np.asarray(flxlst).max()))
        flxmin = axefact*float(np.floor(np.asarray(flxlst).min()))
        for l, (lst, fn, title) in enumerate(zip(flxlst,plt_exportBAL_fn, plt_titleBAL)):
            MMplot.plotGWbudget(flxlst = lst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = fn, plt_title = title, fluxmax = flxmax, fluxmin = flxmin, unit = plt_WB_unit)
        del flxlst
        del h_satflow, MM, MM_S
        h5_MM.close()
        h5_MF.close()
        # output for PEST
        outPESTheads.close()
        outPESTsm.close()
        del obs, obs_h, obs_S

    if plt_out == 1:
        SP_lst = []
        Date_lst = []
        JD_lst = []
        SP = 0
        while SP < len(h_MF_m):
            SP_lst.append(SP)
            Date_lst.append(inputDate[SP])
            JD_lst.append(JD[SP])
            SP += plt_freq
        if tTgmin < 0:
            lst = [len(h_MF_m)-1, tRCHmax]
        else:
            lst = [len(h_MF_m)-1, tRCHmax, tTgmin]
        for e in lst:
            SP_lst.append(e)
            Date_lst.append(inputDate[e])
            JD_lst.append(JD[e])

    # plot MF output
    if plt_out == 1 and isinstance(cMF.h5_MF_fn, str):
        # plot heads (grid + contours), DRN, etc... at specified SP
        h5_MF = h5py.File(cMF.h5_MF_fn, 'r')
        # plot for selected time step
        t = 0
        for SP in SP_lst:
            # plot heads [m]
            V = []
            for L in range(cMF.nlay):
                V.append(h_MF_m[SP,:,:,L])
            if hmaxMF == hminMF:
                if hmaxMF < 10E-9:
                    hmaxMF = 1.0
                    hminMF = -1.0
                else:
                    hmaxMF *= 1.15
                    hminMF *= 0.85
                ctrs_tmp = False
            else:
                ctrs_tmp = ctrsMF
            MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation (m)', msg = 'DRY', plt_title = 'MF_HEADS', MM_ws = MM_ws, interval_type = 'arange', interval_diff = (hmaxMF - hminMF)/nrangeMF, contours = ctrs_tmp, Vmax = hmaxMF, Vmin = hminMF, ntick = ntick)
            # plot GW drainage [mm]
            if cMF.drn_yn == 1:
                V = []
                cbc_DRN = h5_MF['DRN_d']
                for L in range(cMF.nlay):
                    V.append(np.ma.masked_array(cbc_DRN[SP,:,:,L], mask[L])*(-1.0))
                DRNmax_tmp = np.ma.max(V)
                DRNmin_tmp = np.ma.min(V)
                if DRNmax_tmp == DRNmin_tmp:
                    if DRNmax_tmp < 10E-9:
                        DRNmax_tmp = 1.0
                        DRNmin_tmp = -1.0
                    else:
                        DRNmax_tmp *= 1.15
                        DRNmin_tmp *= 0.85
                    ctrs_tmp = False
                else:
                    ctrs_tmp = ctrsMF
                MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage (mm/day)', msg = '- no drainage', plt_title = 'MF_DRN', MM_ws = MM_ws, interval_type = 'linspace', interval_num = 5, Vmin = DRNmin_tmp, contours = ctrs_tmp, Vmax = DRNmax_tmp, fmt='%.3G', ntick = ntick)
            # plot GW RCH [mm]
            V = []
            cbc_RCH = h5_MF['RCH_d']
            for L in range(cMF.nlay):
                V.append(np.ma.masked_array(cbc_RCH[SP,:,:,L], mask[L]))
            RCHmax_tmp = np.ma.max(V)
            RCHmin_tmp = np.ma.min(V)
            if RCHmax_tmp == RCHmin_tmp:
                if RCHmax_tmp < 10E-9:
                    RCHmax_tmp = 1.0
                    RCHmin_tmp = -1.0
                else:
                    RCHmax_tmp *= 1.15
                    RCHmin_tmp *= 0.85
                ctrs_tmp = False
            else:
                ctrs_tmp = ctrsMF
            MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater recharge (mm/day)', msg = '- no flux', plt_title = 'MF_RCH', MM_ws = MM_ws, interval_type = 'linspace', interval_num = 5, Vmin = RCHmin_tmp, contours = ctrs_tmp, Vmax = RCHmax_tmp, fmt='%.3G', ntick = ntick)
            t += 1
            del V
        del t
        # plot average for the whole simulated period
        # plot GW drainage [mm]
        if cMF.drn_yn == 1:
            V = []
            for L in range(cMF.nlay):
                V.append(np.ma.masked_array(np.sum(cbc_DRN[:,:,:,L], axis = 0)/sum(cMF.perlen)*(-1.0), mask[L]))
            DRNmax_tmp = np.ma.max(V)
            DRNmin_tmp = np.ma.min(V)
            if DRNmax_tmp == DRNmin_tmp:
                if DRNmax_tmp < 10E-9:
                    DRNmax_tmp = 1.0
                    DRNmin_tmp = -1.0
                else:
                    DRNmax_tmp *= 1.15
                    DRNmin_tmp *= 0.85
                ctrs_tmp = False
            else:
                ctrs_tmp = ctrsMF
            MMplot.plotLAYER(SP = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage (mm/day)', msg = '- no drainage', plt_title = 'MF_average_DRN', MM_ws = MM_ws, interval_type = 'arange', interval_diff = (DRNmax_tmp - DRNmin_tmp)/nrangeMM, Vmax = DRNmax_tmp, Vmin = DRNmin_tmp, contours = ctrs_tmp, ntick = ntick)
        # plot GW RCH [mm]
        V = []
        for L in range(cMF.nlay):
            V.append(np.ma.masked_array(np.sum(cbc_RCH[:,:,:,L], axis = 0)/sum(cMF.perlen), mask[L]))
        RCHmax_tmp = np.ma.max(V)
        RCHmin_tmp = np.ma.min(V)
        if RCHmax_tmp == RCHmin_tmp:
            if RCHmax_tmp < 10E-9:
                RCHmax_tmp = 1.0
                RCHmin_tmp = -1.0
            else:
                RCHmax_tmp *= 1.15
                RCHmin_tmp *= 0.85
            ctrs_tmp = False
        else:
            ctrs_tmp = ctrsMF
        MMplot.plotLAYER(SP = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = cMF.nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater recharge (mm/day)', msg = '- no flux', plt_title = 'MF_average_RCH', MM_ws = MM_ws, interval_type = 'arange', interval_diff = (RCHmax_tmp - RCHmin_tmp)/nrangeMM, Vmax = RCHmax_tmp, Vmin = RCHmin_tmp, contours = ctrs_tmp, ntick = ntick)
        h5_MF.close()
        del V, cbc_RCH
        if cMF.drn_yn == 1:
            del cbc_DRN, DRNmax_tmp, DRNmin_tmp
        del RCHmax_tmp, RCHmin_tmp

    # plot MM output
    if plt_out == 1 and os.path.exists(h5_MM_fn):
        flxlbl = ['RF', 'RFe', 'I', 'EXF', 'dSs', 'Ro', 'Es', 'Eg', 'Tg', 'ETg', 'ETu', 'dSu']
        for i in flxlbl:
            # plot average for the whole simulated period
            i1 = 'i'+i
            h5_MM = h5py.File(h5_MM_fn, 'r')
            MM = h5_MM['MM'][:,:,:,index.get(i1)]
            h5_MM.close()
            V = [np.sum(np.ma.masked_values(MM, cMF.hnoflo, atol = 0.09), axis = 0)/sum(cMF.perlen)]
            V[0] = np.ma.masked_values(V[0], cMF.hnoflo, atol = 0.09)
            Vmax = np.ma.max(V[0]) #float(np.ceil(np.ma.max(V)))
            Vmin = np.ma.min(V[0]) #float(np.floor(np.ma.min(V)))
            if Vmax == Vmin:
                if Vmax < 10E-9:
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
            MMplot.plotLAYER(SP = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_average_' + i), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
            del V
            # plot for selected time step
            t = 0
            for SP in SP_lst:
                V = [np.ma.masked_values(MM[SP,:,:], cMF.hnoflo, atol = 0.09)]
                Vmax = np.ma.max(V[0]) #float(np.ceil(np.ma.max(V)))
                Vmin = np.ma.min(V[0]) #float(np.floor(np.ma.min(V)))
                if Vmax == Vmin:
                    if Vmax < 10E-9:
                        Vmax = 1.0
                        Vmin = -1.0
                    else:
                        Vmax *= 1.15
                        Vmin *= 0.85
                    ctrs_tmp = False
                else:
                    ctrs_tmp = ctrsMM
                MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_'+i), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
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
            Vmax = np.ma.max(V) #float(np.ceil(np.ma.max(V)))
            Vmin = np.ma.min(V) #float(np.floor(np.ma.min(V)))
            if Vmax == Vmin:
                if Vmax < 10E-9:
                    Vmax = 1.0
                    Vmin = -1.0
                else:
                    Vmax *= 1.15
                    Vmin *= 0.85
                ctrs_tmp = False
            else:
                ctrs_tmp = ctrsMM
            MMplot.plotLAYER(SP = 'NA', Date = 'NA', JD = 'NA', ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_'+ 'average_' + i), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
            del V
            # plot for selected time step
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
                Vmax = np.ma.max(V) #float(np.ceil(np.ma.max(V)))
                Vmin = np.ma.min(V) #float(np.floor(np.ma.min(V)))
                if Vmax == Vmin:
                    if Vmax < 10E-9:
                        Vmax = 1.0
                        Vmin = -1.0
                    else:
                        Vmax *= 1.15
                        Vmin *= 0.85
                    ctrs_tmp = False
                else:
                    ctrs_tmp = ctrsMM
                MMplot.plotLAYER(SP = SP, Date = Date_lst[t], JD = JD_lst[t], ncol = cMF.ncol, nrow = cMF.nrow, nlay = cMF.nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_'+i), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrs_tmp, ntick = ntick)
                t += 1
            del V, MM, t
        del SP_lst, flxlbl, i, i1, h_diff_surf
    del gridSOIL, inputDate
    del hmaxMF, hminMF, hmin, hdiff, cbcmax_d, cbcmin_d
    if cMF.drn_yn == 1:
        del DRNmax, DRNmin

timeendExport = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
durationExport=(timeendExport-timestartExport)
durationTotal = (timeendExport-timestart)

# final report of successful run
print '\n##############\nMARMITES executed successfully!\n%s' % mpl.dates.datetime.datetime.today().isoformat()[:19]
if MMunsat_yn > 0:
    print '\nLOOP %d/%d' % (LOOP-1, ccnum)
    for txt in msg_end_loop:
        print txt
print '\n%d stress periods, %d days' % (cMF.nper,sum(cMF.perlen))
print '%d rows x %d cols (%d cells) x %d layers' % (cMF.nrow, cMF.ncol, cMF.nrow*cMF.ncol, cMF.nlay)
print '%d MM active cells in total' % (sum(ncell_MM))
l = 1
for n in ncell:
    print  'LAYER %d' % l
    print '%d MF active cells' % (n)
    print '%d MM active cells' % (ncell_MM[l-1])
    l += 1
print ('\nApproximate run times:')
if MMsurf_yn > 0:
    print ('MARMITES surface: %s minute(s) and %.1f second(s)') % (str(int(durationMMsurf*24.0*60.0)), (durationMMsurf*24.0*60.0-int(durationMMsurf*24.0*60.0))*60)
if MMunsat_yn > 0:
    print ('MARMITES unsaturated zone: %s minute(s) and %.1f second(s)') % (str(int(durationMMunsat*24.0*60.0)), (durationMMunsat*24.0*60.0-int(durationMMunsat*24.0*60.0))*60)
if MF_yn == 1:
    print ('MODFLOW: %s minute(s) and %.1f second(s)') % (str(int(durationMF*24.0*60.0)), (durationMF*24.0*60.0-int(durationMF*24.0*60.0))*60)
print ('Export: %s minute(s) and %.1f second(s)') % (str(int(durationExport*24.0*60.0)), (durationExport*24.0*60.0-int(durationExport*24.0*60.0))*60)
print ('Total: %s minute(s) and %.1f second(s)') % (str(int(durationTotal*24.0*60.0)), (durationTotal*24.0*60.0-int(durationTotal*24.0*60.0))*60)
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