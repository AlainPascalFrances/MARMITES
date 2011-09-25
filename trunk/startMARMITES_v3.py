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
import numpy as np
import startMARMITESsurface as startMMsurf
import MARMITESunsat_v3 as MMunsat
import MARMITESprocess_v3 as MMproc
import ppMODFLOW_flopy_v3 as ppMF
import MARMITESplot_v3 as MMplot
import CreateColors


#####################################

# workspace (ws) definition
timestart = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
print '\n##############\nMARMITES started!\n%s\n##############' % mpl.dates.num2date(timestart).isoformat()[:19]

# read input file (called _input.ini in the MARMITES workspace
# the first character on the first line has to be the character used to comment
# the file can contain any comments as the user wish, but the sequence of the input has to be respected
MM_ws = r'E:\00code_ws\00_TESTS\MARMITESv3_r13c6l2'
MM_fn = '__inputMM.ini'

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
    obsCHECK = int(inputFile[l].strip())
    #run MARMITESsurface  (1 is YES, 0 is NO)
    l += 1
    MMsurf_yn = int(inputFile[l].strip())
    #run MARMITESunsat  (1 is YES, 0 is NO)
    l += 1
    MMunsat_yn = int(inputFile[l].strip())
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
except:
    print '\nType error in the input file %s' % (MM_fn)
    sys.exit()
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
        outMMsurf_fn = startMMsurf.MMsurf(MMsurf_ws, inputFile_TS_fn, inputFile_PAR_fn, outputFILE_fn, MM_ws, outMMsurf_fn)

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
        inputZON_TS_RF_fn = str(inputFile[l].strip())
        l += 1
        inputZON_TS_PET_fn = str(inputFile[l].strip())
        l += 1
        inputZON_TS_RFe_fn = str(inputFile[l].strip())
        l += 1
        inputZON_TS_PE_fn = str(inputFile[l].strip())
        l += 1
        inputZON_TS_E0_fn = str(inputFile[l].strip())
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
        print '\nType error in file [' + inputFile_fn + ']'
        sys.exit()
    del inputFile

    numDays = len(MMproc.readFile(MM_ws, inputDate_fn))

    # #############################
    # ###  READ MODFLOW CONFIG ####
    # #############################

    print'\n##############'
    print 'MODFLOW initialization'
    nrow, ncol, delr, delc, top, reggrid, nlay, nper, perlen, nstp, timedef, hnoflo, hdry, laytyp, lenuni, itmuni, ibound_fn = ppMF.ppMFini(MF_ws, MF_ini_fn, out = 'MM', numDays = numDays)

    # ####   SUMMARY OF MODFLOW READINGS   ####
    # active cells in layer 1_____________________ibound[0]
    # elevation___________________________________top
    # aquifer type of layer 1_____________________laytyp[0]
    # numero total de time step___________________sum(nstp)
    # code for dry cell___________________________hdry
    # cell size___________________________________delr[0]

    # compute cbc conversion factor from volume to mm
    if lenuni == 1:
        conv_fact = 304.8
    elif lenuni == 2:
        conv_fact = 1000.0
    elif lenuni == 3:
        conv_fact = 10.0
    else:
        print 'FATAL ERROR!\nDefine the length unit in the MODFLOW ini file!\n (see USGS Open-File Report 00-92)'
        sys.exit()
        # TODO if lenuni!=2 apply conversion factor to delr, delc, etc...
    if laytyp[0]==0:
        print 'FATAL ERROR!\nThe first layer cannot be confined type!\nChange your parameter laytyp in the MODFLOW lpf package.\n(see USGS Open-File Report 00-92)'
        sys.exit()
    if itmuni != 4:
        print 'FATAL ERROR! Time unit is not in days!'
        sys.exit()

    # #############################
    # ### MF time processing
    # #############################
    # if required by user, compute nper, perlen,etc based on RF analysis in the METEO zones
    if timedef>=0:
        if isinstance(nper, str):
            try:
                perlenmax = int(nper.split()[1].strip())
            except:
                print '\nError in nper format of the MODFLOW ini file!\n'
                sys.exit()
        MFtime_fn = 'MFtimediscret.txt'
        nper, perlen, nstp, inputZON_TS_RF_fn, inputZON_TS_PET_fn, inputZON_TS_RFe_fn, inputZON_TS_PE_fn, inputZON_TS_E0_fn = ppMF.ppMFtime(MM_ws, MF_ws, MFtime_fn, nper, perlen, nstp, timedef, inputDate_fn, inputZON_TS_RF_fn, inputZON_TS_PET_fn, inputZON_TS_RFe_fn, inputZON_TS_PE_fn, inputZON_TS_E0_fn, NMETEO, NVEG, NSOIL)
    else:
        MFtime_fn = None

    print'\n##############'
    print 'MARMITESunsat initialization'

    # MARMITES INITIALIZATION
    MM_PROCESS = MMproc.PROCESS(MM_ws                = MM_ws,
                            MF_ws                    = MF_ws,
                            nrow                     = nrow,
                            ncol                     = ncol,
                            xllcorner                = xllcorner,
                            yllcorner                = yllcorner,
                            cellsizeMF               = delr[0],
                            nstp                     = nstp,
                            hnoflo                   = hnoflo
                            )
    MM_UNSAT = MMunsat.UNSAT(hnoflo = hnoflo)
    MM_SATFLOW = MMunsat.SATFLOW()

    ibound = np.zeros((nrow,ncol), dtype = int)
    ibound_path = os.path.join(MF_ws, ibound_fn[0])
    MM_PROCESS.convASCIIraster2array(ibound_path, ibound)
    ncell = (ibound[:,:]!=0).sum()

    # READ input ESRI ASCII rasters # missing gridIRR_fn
    print "\nImporting ESRI ASCII files to initialize MARMITES..."
    gridMETEO = MM_PROCESS.inputEsriAscii(grid_fn = gridMETEO_fn, datatype = int)

    gridSOIL = MM_PROCESS.inputEsriAscii(grid_fn = gridSOIL_fn, datatype = int)

    gridSOILthick = MM_PROCESS.inputEsriAscii(grid_fn = gridSOILthick_fn,
     datatype = float)

    gridSshmax = MM_PROCESS.inputEsriAscii(grid_fn = gridSshmax_fn,
     datatype = float)

    gridSsw = MM_PROCESS.inputEsriAscii(grid_fn = gridSsw_fn,
     datatype = float)

    ##gridIRR = MM_PROCESS.inputEsriAscii(grid_fn                  = gridIRR_fn)

    # READ input time series and parameters   # missing IRR_fn
    gridVEGarea, RFzonesTS, E0zonesTS, PETvegzonesTS, RFevegzonesTS, PEsoilzonesTS, inputDate = MM_PROCESS.inputTS(NMETEO = NMETEO,
                                    NVEG                     = NVEG,
                                    NSOIL                    = NSOIL,
                                    inputDate_fn             = inputDate_fn,
                                    inputZON_TS_RF_fn        = inputZON_TS_RF_fn,
                                    inputZON_TS_PET_fn       = inputZON_TS_PET_fn,
                                    inputZON_TS_RFe_fn       = inputZON_TS_RFe_fn,
                                    inputZON_TS_PE_fn        = inputZON_TS_PE_fn,
                                    inputZON_TS_E0_fn        = inputZON_TS_E0_fn
     ) # IRR_fn

    # SOIL PARAMETERS
    _nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _Su_ini, _Ks = MM_PROCESS.inputSoilParam(MM_ws = MM_ws, SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
    _nslmax = max(_nsl)

    try:
        top = float(top[0])
        top_array = np.ones((nrow,ncol))*top
    except:
        if isinstance(top[0], str):
            top_path = os.path.join(MF_ws, top[0])
            top_array = np.zeros((nrow,ncol))
            top_array = MM_PROCESS.convASCIIraster2array(top_path, top_array)

    # READ observations time series (heads and soil moisture)
    if obsCHECK==1:
        print "\nReading observations time series (hydraulic heads and soil moisture)..."
        obs, outpathname, obs_h, obs_S = MM_PROCESS.inputObs(
                                MM_ws = MM_ws,
                                inputObs_fn = inputObs_fn,
                                inputObsHEADS_fn = inputObsHEADS_fn,
                                inputObsSM_fn = inputObsSM_fn,
                                inputDate   = inputDate,
                                _nslmax     = _nslmax,
                                nlay        = nlay
                                )
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
    else:
        print "\nNo reading of observations time series (hydraulic heads and soil moisture) required."
        obs = None

    # compute thickness, top and bottom elevation of each soil layer
    TopAquif = top_array*1000.0 # conversion from m to mm
    # topography elevation
    TopSoil = TopAquif + gridSOILthick*1000.0

    # indexes of the HDF5 output arrays
    index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iEXF':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSs':13, 'iETg':14, 'iETu':15, 'idSu':16, 'iHEADScorr':17, 'idtwt':18, 'iuzthick':19}
    index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'iRexf':4, 'idSu':5, 'iSu':6, 'iSAT':7, 'iMB_l':8}

    # #############################
    # ### 1st MODFLOW RUN with initial user-input recharge
    # #############################
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
    h5_MM_fn = os.path.join(MM_ws,'_h5_MM.h5')
    h5_MF_fn = ppMF.ppMF(MM_ws, xllcorner, yllcorner, MF_ws, MF_ini_fn, finf_MM = (h5_MM_fn, 'finf'), wel_MM = (h5_MM_fn, 'ETg'), report = report, verbose = verbose, chunks = chunks, MFtime_fn = MFtime_fn, numDays = numDays, obs = obs)

    h5_MF = h5py.File(h5_MF_fn)

    # heads format is : timestep, nrow, ncol, nlay
    # cbc format is: (kstp), kper, textprocess, nrow, ncol, nlay
    cbc_nam = []
    cbc_uzf_nam = []
    for c in h5_MF['cbc_nam']:
        cbc_nam.append(c.strip())
    for c in h5_MF['cbc_uzf_nam']:
        cbc_uzf_nam.append(c.strip())
    imfSTO = cbc_nam.index('STORAGE')
    imfDRN = cbc_nam.index('DRAINS')
    imfWEL = cbc_nam.index('WELLS')
    imfEXF = cbc_uzf_nam.index('SURFACE LEAKAGE')
    imfRCH = cbc_uzf_nam.index('UZF RECHARGE')

    timeendMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    durationMF +=  timeendMF-timestartMF

except StandardError, e:  #Exception
    print '\nERROR! Abnormal MM run interruption in the initialization!\nError description:'
    traceback.print_exc(file=sys.stdout)
    h5_MF_fn = None
#    traceback.print_exc(limit=1, file=sys.stdout)
    timeend = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    duration += (timeend-timestart) -durationMF
    sys.exit()

# #############################
# 2nd phase : MM/MF loop #####
# #############################
duration = 0.0
h_diff_surf = None
try:
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
            # ###########################
            # ###  MARMITES INPUT #######
            # ###########################

            print'\n##############'
            print 'MARMITESunsat RUN'

            # SOIL PARAMETERS
            _nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _Su_ini, _Ks = MM_PROCESS.inputSoilParam(MM_ws = MM_ws, SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
            _nslmax = max(_nsl)

            # ###############
            # Create HDF5 arrays to store MARMITES output
            h5_MM = h5py.File(h5_MM_fn, 'w')
            # arrays for fluxes independent of the soil layering
            h5_MM.create_dataset(name = 'iMM', data = np.asarray(index.items()))
            if chunks == 1:
                h5_MM.create_dataset(name = 'MM', shape = (sum(perlen),nrow,ncol,len(index)), dtype = np.float, chunks = (1,nrow,ncol,len(index)),  compression = 'gzip', compression_opts = 5, shuffle = True)
            else:
                h5_MM.create_dataset(name = 'MM', shape = (sum(perlen),nrow,ncol,len(index)), dtype = np.float)
            # arrays for fluxes in each soil layer
            h5_MM.create_dataset(name = 'iMM_S', data = np.asarray(index_S.items()))
            if chunks == 1:
                h5_MM.create_dataset(name = 'MM_S', shape = (sum(perlen),nrow,ncol,_nslmax,len(index_S)), dtype = np.float, chunks = (1,nrow,ncol,_nslmax,len(index_S)),  compression = 'gzip', compression_opts = 5, shuffle = True)
            else:
                h5_MM.create_dataset(name = 'MM_S', shape = (sum(perlen),nrow,ncol,_nslmax,len(index_S)), dtype = np.float)
            # arrays to compute net recharge to be exported to MF
            if chunks == 1:
                h5_MM.create_dataset(name = 'finf', shape = (nper,nrow,ncol), dtype = np.float, chunks = (1,nrow,ncol),  compression = 'gzip', compression_opts = 5, shuffle = True)
                h5_MM.create_dataset(name = 'ETg', shape = (nper,nrow,ncol), dtype = np.float, chunks = (1,nrow,ncol),  compression = 'gzip', compression_opts = 5, shuffle = True)
            else:
                h5_MM.create_dataset(name = 'finf', shape = (nper,nrow,ncol), dtype = np.float)
                h5_MM.create_dataset(name = 'ETg', shape = (nper,nrow,ncol), dtype = np.float)
            # ###############
            # # main loop: calculation of soil water balance in each cell-grid for each time step inside each stress period
    #        t0=0
            print '\nComputing...'

            # initial values of SP
            Su_ini_tmp_array = np.zeros([nrow,ncol,_nslmax])
            Rp_ini_tmp_array = np.zeros([nrow,ncol,_nslmax])
            Ss_ini_tmp_array = np.zeros([nrow,ncol])
            Ss_ini_tmp_array = np.zeros([nrow,ncol])
            for n in range(nper):
                tstart_MM = 0
                for t in range(n):
                    tstart_MM += perlen[t]
                tend_MM = tstart_MM + perlen[n]
                tstart_MF = 0
                for t in range(n):
                    tstart_MF += nstp[t]
                tend_MF = tstart_MF + nstp[n]
                MM = np.zeros([perlen[n],nrow,ncol,len(index)], dtype=float)
                MM_S = np.zeros([perlen[n],nrow,ncol,_nslmax,len(index_S)], dtype=float)
                MM_finf_MF = np.zeros([nrow,ncol], dtype=float)
                MM_wel_MF = np.zeros([nrow,ncol], dtype=float)
                h_MF_per = h5_MF['heads'][tstart_MF:tend_MF,:,:,0]
                exf_MF_per = h5_MF['cbc_uzf'][tstart_MF:tend_MF,imfEXF,:,:,0]
                # loop into the grid
                for i in range(nrow):
                    for j in range(ncol):
                        SOILzone_tmp = gridSOIL[i,j]-1
                        METEOzone_tmp = gridMETEO[i,j]-1
                        if ibound[i,j]!=0.0:
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
                                                         TopSoilLay   = TopSoilLay,
                                                         BotSoilLay   = BotSoilLay,
                                                         Tl      = Tl,
                                                         Sm      = _Sm[SOILzone_tmp],
                                                         Sfc     = _Sfc[SOILzone_tmp],
                                                         Sr      = _Sr[SOILzone_tmp],
                                                         Su_ini  = Su_ini_tmp,
                                                         Ss_ini  = Ss_ini_tmp,
                                                         Rp_ini  = Rp_ini_tmp_array[i,j,:],
                                                         Ks      = _Ks[SOILzone_tmp],
                                                         Ss_max  = 1000*1.12*gridSshmax[i,j]*gridSsw[i,j]/delr[j],
                                                         Ss_ratio= 1.12*gridSsw[i,j]/delr[j],
                                                         HEADS   = h_MF_tmp,
                                                         EXF     = exf_MF_tmp,
                                                         RF      = RFzonesTS[METEOzone_tmp][tstart_MF:tend_MF],
                                                         E0      = E0zonesTS[METEOzone_tmp][tstart_MF:tend_MF],
                                                         PETveg  = PETvegzonesTS_tmp,
                                                         RFeveg  = RFevegzonesTS_tmp,
                                                         PEsoil  = PEsoilzonesTS_tmp,
                                                         VEGarea = VEGarea_tmp,
                                                         Zr      = Zr,
                                                         nstp    = nstp[n],
                                                         perlen  = perlen[n],
                                                         dti     = dti,
                                                         hdry    = hdry,
                                                         kTu_min = kTu_min,
                                                         kTu_n = kTu_n)
                            if (float(perlen[n])/float(nstp[n]))!=1.0:
                                for stp in range(nstp[n]):
                                    ts = float(perlen[n])/float(nstp[n])
                                    tstart =    stp*ts
                                    tend =   (1+stp)*ts
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
                            Su_ini_tmp_array[i,j,:]  = MM_S[nstp[n]-1,i,j,:,index_S.get('iSu')]
                            Rp_ini_tmp_array[i,j,:]  = MM_S[nstp[n]-1,i,j,:,index_S.get('iRp')]
                            Ss_ini_tmp_array[i,j]    = MM[nstp[n]-1,i,j,index.get('iSs')]
                        else:
                            if perlen[n]>1:
                                MM[:,i,j,:] = hnoflo
                                MM_S[:,i,j,:,:] = hnoflo
                            else:
                                MM[:,i,j,:] = hnoflo
                                MM_S[:,i,j,:,:] = hnoflo
                            MM_finf_MF[i,j] = hnoflo
                            MM_wel_MF[i,j] = hnoflo
                dti = float(perlen[n])/float(nstp[n])
                h5_MM['MM'][tstart_MM:tend_MM,:,:,:]     = MM[:,:,:,:]
                h5_MM['MM_S'][tstart_MM:tend_MM,:,:,:,:] = MM_S[:,:,:,:,:]
                h5_MM['finf'][n,:,:] = MM_finf_MF
                h5_MM['ETg'][n,:,:] = MM_wel_MF
            del MM, MM_S, MM_finf_MF, MM_wel_MF, exf_MF_per, Su_ini_tmp_array, Rp_ini_tmp_array, Ss_ini_tmp_array, dti
            h5_MM.close()

            # CHECK MM amd MF CONVERG.
            h_MF_average = np.ma.average(np.ma.masked_values(h_MF_per, hnoflo, atol = 0.09))
            h_diff.append(h_MF_average - h_pSP)
            h_diff_surf = np.ma.masked_values(h_MF_per, hnoflo, atol = 0.09) - np.ma.masked_values(h_pSP_all, hnoflo, atol = 0.09)
            h_diff_all_max = np.nanmax(h_diff_surf)
            h_diff_all_min = np.nanmin(h_diff_surf)
            if abs(h_diff_all_max)>abs(h_diff_all_min):
                h_diff_all.append(h_diff_all_max)
            else:
                h_diff_all.append(h_diff_all_min)
            LOOP += 1
            LOOPlst.append(LOOP)
            h_pSP = h_MF_average
            h_pSP_all = h_MF_per
            del h_MF_per
            if np.absolute(h_diff[LOOP])>0.0:
                h_diff_log.append(np.log10(np.absolute(h_diff[LOOP])))
                h_diff_all_log.append(np.log10(np.absolute(h_diff_all[LOOP])))
            else:
                h_diff_log.append(np.log10(convcrit))
                h_diff_all_log.append(np.log10(convcrit))
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
            h5_MF_fn = ppMF.ppMF(MM_ws, xllcorner, yllcorner, MF_ws, MF_ini_fn, finf_MM = (h5_MM_fn, 'finf'), wel_MM = (h5_MM_fn, 'ETg'), report = report, verbose = verbose, chunks = chunks, MFtime_fn = MFtime_fn, numDays = numDays, obs = obs)

            h5_MF = h5py.File(h5_MF_fn)

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

            ax2=fig.add_subplot(3,1,3, sharex = ax1)
            plt.setp(ax2.get_xticklabels(), fontsize=8)
            plt.setp(ax2.get_yticklabels(), fontsize=8)
            ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
            plt.ylabel('log(abs(h_diff)) [log(m)]', fontsize=10, horizontalalignment = 'center')
            plt.grid(True)
            plt.xlabel('trial', fontsize=10)
            plt.plot(LOOPlst[2:], h_diff_log[2:], linestyle='-', marker='o', markersize=5, c = 'orange', markerfacecolor='orange', markeredgecolor='red')
            plt.plot(LOOPlst[2:], h_diff_all_log[2:], linestyle='-', marker='o', markersize=5, c = 'green', markerfacecolor='green', markeredgecolor='blue')

            ax3=fig.add_subplot(3,1,2, sharex = ax1)
            plt.setp(ax3.get_xticklabels(), fontsize=8)
            plt.setp(ax3.get_yticklabels(), fontsize=8)
            plt.plot(LOOPlst[2:], h_diff[2:], linestyle='-', marker='o', markersize=5, c = 'orange', markerfacecolor='orange', markeredgecolor='red')
            plt.plot(LOOPlst[2:], h_diff_all[2:], linestyle='-', marker='o', markersize=5, c = 'green', markerfacecolor='green', markeredgecolor='blue')
            ax3.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.3G'))
            plt.ylabel('h_diff [m]', fontsize=10, horizontalalignment = 'center')
        #        plt.ylabel.Text.position(0.5, -0.5)
            plt.grid(True)

            ax2.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%2d'))
            ax3.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%2d'))

            plt.xlim(1,LOOP)
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
        h5_MF_fn = ppMF.ppMF(MM_ws, xllcorner, yllcorner, MF_ws, MF_ini_fn, finf_MM = (h5_MM_fn, 'finf'), wel_MM = (h5_MM_fn, 'ETg'), report = report, verbose = verbose, chunks = chunks, MFtime_fn = MFtime_fn, numDays = numDays, obs = obs)

        timeendMF = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
        durationMF += (timeendMF-timestartMF)

except StandardError, e:  #Exception
    print '\nERROR! Abnormal MM run interruption in the MM/MF loop!\nError description:'
    h5_MF_fn = None
    traceback.print_exc(file=sys.stdout)
    timeend = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    duration += (timeend-timestart) -durationMF
    sys.exit()

# #############################
# 3rd phase : export results #####
# #############################

try:
    del gridVEGarea, RFzonesTS, E0zonesTS, PETvegzonesTS, RFevegzonesTS, PEsoilzonesTS
    del gridMETEO, gridSOILthick, gridSshmax, gridSsw

    # #############################
    # ###  OUTPUT EXPORT   ########
    # #############################

    timestartExport = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    print '\n##############\nMARMITES exporting...'

    # reading MF output
    if isinstance(h5_MF_fn, str):
        h5_MF = h5py.File(h5_MF_fn)
        cbc_DRN = np.zeros((sum(perlen), nrow, ncol, nlay))
        cbc_STO = np.zeros((sum(perlen), nrow, ncol, nlay))
        cbc_RCH = np.zeros((sum(perlen), nrow, ncol, nlay))
        cbc_WEL = np.zeros((sum(perlen), nrow, ncol, nlay))
        # cbc format is : (kstp), kper, textprocess, nrow, ncol, nlay
        t = 0
        if timedef>=0:
            h_MF = np.zeros((sum(perlen), nrow, ncol, nlay))
            for n in range(nper):
                if perlen[n]!=1.0:
                    for x in range(perlen[n]):
                        h_MF[t,:,:,:] = h5_MF['heads'][n,:,:,:]
                        cbc_DRN_tmp = h5_MF['cbc'][n,imfDRN,:,:,:]
                        cbc_STO_tmp = h5_MF['cbc'][n,imfSTO,:,:,:]
                        cbc_RCH_tmp = h5_MF['cbc_uzf'][n,imfRCH,:,:,:]
                        cbc_WEL_tmp = h5_MF['cbc'][n,imfWEL,:,:,:]
                        if reggrid == 1:
                            cbc_DRN[t,:,:,:] = conv_fact*cbc_DRN_tmp[:,:,:]/(delr[0]*delc[0])
                            cbc_STO[t,:,:,:] = conv_fact*cbc_STO_tmp[:,:,:]/(delr[0]*delc[0])
                            cbc_RCH[t,:,:,:] = conv_fact*cbc_RCH_tmp[:,:,:]/(delr[0]*delc[0])
                            cbc_WEL[t,:,:,:] = conv_fact*cbc_WEL_tmp[:,:,:]/(delr[0]*delc[0])
                        else:
                            for i in range(nrow):
                                for j in range(ncol):
                                    cbc_DRN[t,i,j,:] = conv_fact*cbc_DRN_tmp[:,i,j,:]/(delr[j]*delc[i])
                                    cbc_STO[t,i,j,:] = conv_fact*cbc_STO_tmp[:,i,j,:]/(delr[j]*delc[i])
                                    cbc_RCH[t,i,j,:] = conv_fact*cbc_RCH_tmp[:,i,j,:]/(delr[j]*delc[i])
                                    cbc_WEL[t,i,j,:] = conv_fact*cbc_WEL_tmp[:,i,j,:]/(delr[j]*delc[i])
                        t += 1
                else:
                    h_MF[t,:,:,:] = h5_MF['heads'][n,:,:,:]
                    cbc_DRN_tmp = h5_MF['cbc'][n,imfDRN,:,:,:]
                    cbc_STO_tmp = h5_MF['cbc'][n,imfSTO,:,:,:]
                    cbc_RCH_tmp = h5_MF['cbc_uzf'][n,imfRCH,:,:,:]
                    cbc_WEL_tmp = h5_MF['cbc'][n,imfWEL,:,:,:]
                    if reggrid == 1:
                        cbc_DRN[t,:,:,:] = conv_fact*cbc_DRN_tmp[:,:,:]/(delr[0]*delc[0])
                        cbc_STO[t,:,:,:] = conv_fact*cbc_STO_tmp[:,:,:]/(delr[0]*delc[0])
                        cbc_RCH[t,:,:,:] = conv_fact*cbc_RCH_tmp[:,:,:]/(delr[0]*delc[0])
                        cbc_WEL[t,:,:,:] = conv_fact*cbc_WEL_tmp[:,:,:]/(delr[0]*delc[0])
                    else:
                        for i in range(nrow):
                            for j in range(ncol):
                                cbc_DRN[t,i,j,:] = conv_fact*cbc_DRN_tmp[:,i,j,:]/(delr[j]*delc[i])
                                cbc_STO[t,i,j,:] = conv_fact*cbc_STO_tmp[:,i,j,:]/(delr[j]*delc[i])
                                cbc_RCH[t,i,j,:] = conv_fact*cbc_RCH_tmp[:,i,j,:]/(delr[j]*delc[i])
                                cbc_WEL[t,i,j,:] = conv_fact*cbc_WEL_tmp[:,i,j,:]/(delr[j]*delc[i])
                    t += 1
            del cbc_DRN_tmp, cbc_STO_tmp, cbc_RCH_tmp, cbc_WEL_tmp
        else:
            h_MF = h5_MF['heads']
            cbc_DRN_tmp = h5_MF['cbc'][:,imfDRN,:,:,:]
            cbc_STO_tmp = h5_MF['cbc'][:,imfSTO,:,:,:]
            cbc_RCH_tmp = h5_MF['cbc_uzf'][:,imfRCH,:,:,:]
            cbc_WEL_tmp = h5_MF['cbc'][:,imfWEL,:,:,:]
            if reggrid == 1:
                cbc_DRN[:,:,:,:] = conv_fact*cbc_DRN_tmp[:,:,:]/(delr[0]*delc[0])
                cbc_STO[:,:,:,:] = conv_fact*cbc_STO_tmp[:,:,:]/(delr[0]*delc[0])
                cbc_RCH[:,:,:,:] = conv_fact*cbc_RCH_tmp[:,:,:]/(delr[0]*delc[0])
                cbc_WEL[:,:,:,:] = conv_fact*cbc_WEL_tmp[:,:,:]/(delr[0]*delc[0])
            else:
                cbc_tmp = h5_MF['cbc'][:,:,i,j,:]
                for i in range(nrow):
                    for j in range(ncol):
                        cbc_DRN[:,i,j,:] = conv_fact*cbc_DRN_tmp[i,j,:]/(delr[0]*delc[0])
                        cbc_STO[:,i,j,:] = conv_fact*cbc_STO_tmp[i,j,:]/(delr[0]*delc[0])
                        cbc_RCH[:,i,j,:] = conv_fact*cbc_RCH_tmp[i,j,:]/(delr[0]*delc[0])
                        cbc_WEL[:,i,j,:] = conv_fact*cbc_WEL_tmp[i,j,:]/(delr[0]*delc[0])
                del cbc_DRN_tmp, cbc_STO_tmp, cbc_RCH_tmp, cbc_WEL_tmp
        h_MF_m = np.ma.masked_values(h_MF, hnoflo, atol = 0.09)
        del h_MF
        h5_MF.close()
        top_array_m = np.ma.masked_values(top_array, hnoflo, atol = 0.09)
        del top_array
        index_cbc_uzf = [imfRCH]
        index_cbc = [imfSTO, imfDRN, imfWEL]
    else:
        h_MF_m = np.zeros((sum(perlen), nrow, ncol, nlay))
        cbc_DRN = cbc_STO = cbc_RCH = cbc_WEL = np.zeros((sum(perlen), nrow, ncol, nlay))
        imfDRN = imfSTO = imfRCH = imfWEL = 0
        top_array_m = np.zeros((nrow, ncol))
        del top_array

    # computing range for plotting
    h_MF_m = np.ma.masked_values(h_MF_m, hdry, atol = 1E+25)
    hmax = []
    hmin = []
    hmaxMM = []
    hminMM = []
    hdiff = []
    cbcmax = []
    cbcmin = []
    if obs != None:
        for o in range(len(obs.keys())):
            i = obs.get(obs.keys()[o])['i']
            j = obs.get(obs.keys()[o])['j']
            l = obs.get(obs.keys()[o])['lay']
            hmaxMF = float(np.ceil(np.nanmax(h_MF_m[:,i,j,l].flatten())))
            hminMF = float(np.floor(np.nanmin(h_MF_m[:,i,j,l].flatten())))
            if obs_h[o] != []:
                npa_m_tmp = np.ma.masked_values(obs_h[o], hnoflo, atol = 0.09)
                hmaxMM = float(np.ceil(np.nanmax(npa_m_tmp.flatten())))
                hminMM = float(np.floor(np.nanmin(npa_m_tmp.flatten())))
                del npa_m_tmp
            else:
                hmaxMM = -999
                hminMM = 999
            hmax.append(np.nanmax((hmaxMF, hmaxMM)))
            hmin.append(np.nanmin((hminMF, hminMM)))
            hdiff.append(hmax[o]-hmin[o])
        hdiff = np.nanmax(hdiff)
    else:
        for o in range(len(obs.keys())):
            hmax.append(999.9)
            hmin.append(-999.9)
        hdiff = 2000
    try:
        hmaxMF = float(np.ceil(np.nanmax(h_MF_m[:,:,:,:].flatten())))
        hminMF = float(np.floor(np.nanmin(h_MF_m[:,:,:,:].flatten())))
    except:
        hmaxMF = 999.9
        hminMF = -999.9
    if isinstance(h5_MF_fn, str):
        DRNmax = np.nanmax(-cbc_DRN)
        cbcmax.append(DRNmax)
        DRNmax = float(np.ceil(np.nanmax(DRNmax)))
        DRNmin = np.nanmin(-cbc_DRN)
        cbcmin.append(DRNmin)
        DRNmin = float(np.floor(np.nanmin(DRNmin)))
        cbcmax.append(np.nanmax(-cbc_STO))
        cbcmax.append(np.nanmax(-cbc_RCH))
        cbcmax.append(np.nanmax(-cbc_WEL))
        cbcmin.append(np.nanmin(-cbc_STO))
        cbcmin.append(np.nanmin(-cbc_RCH))
        cbcmin.append(np.nanmin(-cbc_WEL))
        cbcmax = float(np.ceil(np.nanmax(cbcmax)))
        cbcmin = float(np.floor(np.nanmin(cbcmin)))
    else:
        DRNmax = cbcmax = 1
        DRNmin = cbcmin = -1

    print '\nExporting ASCII files and plots...'
    # plot UNSAT/GW balance at the catchment scale
    if plot_out == 1:
        if os.path.exists(h5_MM_fn):
            h5_MM = h5py.File(h5_MM_fn, 'r')
            # indexes of the HDF5 output arrays
            #    index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iEXF':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSs':13, 'iETg':14, 'iETu':15, 'idSu':16, 'iHEADScorr':17, 'idtwt':18}
            #   index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'iRexf':4, 'idSu':5, 'iSu':6, 'iSAT':7, 'iMB_l':8}
            flxlbl = ['RF', 'I', 'EXF', 'dSs', 'Ro', 'Es', 'dSu']
            flxlbl1 = ['Eu', 'Tu']
            flxlbl1a = ['Rp']
            flxlbl2 = ['ETu', 'Eg', 'Tg', 'ETg']
            sign = [1,-1,1,1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1]
            flxlst = []
            flxmax = []
            flxmin = []
            for i in flxlbl:
                i = 'i'+i
                flxlst.append(np.ma.masked_values(h5_MM['MM'][:,:,:,index.get(i)], hnoflo, atol = 0.09).sum()/sum(perlen)/ncell)
            for i in flxlbl1:
                tmp = 0.0
                flxlbl.append(i)
                i = 'i'+i
                for l in range(_nslmax):
                    tmp += np.ma.masked_values(h5_MM['MM_S'][:,:,:,l,index_S.get(i)], hnoflo, atol = 0.09).sum()
                flxlst.append(tmp/sum(perlen)/ncell)
                del tmp
            for i in flxlbl1a:
                flxlbl.append(i)
                i = 'i'+i
                flxlst.append(np.ma.masked_values(h5_MM['MM_S'][:,:,:,-1,index_S.get(i)], hnoflo, atol = 0.09).sum()/sum(perlen)/ncell)
            for i in flxlbl2:
                flxlbl.append(i)
                i = 'i'+i
                flxlst.append(np.ma.masked_values(h5_MM['MM'][:,:,:,index.get(i)], hnoflo, atol = 0.09).sum()/sum(perlen)/ncell)
            h5_MM.close()
            for l, (x,y) in enumerate(zip(flxlst, sign)):
                flxlst[l] = x*y
            del flxlbl1, flxlbl1a, flxlbl2, sign
            flxmax = float(np.ceil(np.nanmax(flxlst)))
            flxmin = float(np.floor(np.nanmin(flxlst)))
            flxmax = 1.15*max(cbcmax, flxmax)
            flxmin = 1.15*min(cbcmin, flxmin)
            if isinstance(h5_MF_fn, str):
                plt_export_fn = os.path.join(MM_ws, '_plt_0UNSATandGWbalances.png')
                for l in range(nlay):
                    flxlst.append(cbc_RCH[:,:,:,l].sum()/sum(perlen)/ncell)
                    flxlst.append(cbc_STO[:,:,:,l].sum()/sum(perlen)/ncell)
                    flxlst.append(cbc_DRN[:,:,:,l].sum()/sum(perlen)/ncell)
                    flxlst.append(cbc_WEL[:,:,:,l].sum()/sum(perlen)/ncell)
                    for x in range(len(index_cbc_uzf)):
                        flxlbl.append('GW_' + cbc_uzf_nam[index_cbc_uzf[x]] + '_L' + str(l+1))
                    for x in range(len(index_cbc)):
                        flxlbl.append('GW_' + cbc_nam[index_cbc[x]] + '_L' + str(l+1))
                # MB_tmp = sum(flxlst) - (flxlst[11] + flxlst[13] + flxlst[16] + flxlst[17] + flxlst[19] + flxlst[20])
                MB_tmp = -999
                plt_title = 'MARMITES and MODFLOW water flux balance for the whole catchment\nsum of fluxes: %.4f' % MB_tmp
            else:
                plt_export_fn = os.path.join(MM_ws, '_plt_0UNSATbalance.png')
                # MB_tmp = sum(flxlst) - (flxlst[11] + flxlst[13])
                MB_tmp = -999
                plt_title = 'MARMITES water flux balance for the whole catchment\nsum of fluxes: %.4f' % MB_tmp
            colors_flx = CreateColors.main(hi=0, hf=180, numbcolors = len(flxlbl))
            MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title, fluxmax = flxmax, fluxmin = flxmin)
            del flxlst
        # no MM, plot only MF balance if MF exists
        elif isinstance(h5_MF_fn, str):
            plt_export_fn = os.path.join(MM_ws, '_plt_GWbalances.png')
            # MB_tmp = sum(flxlst) - (flxlst[11] + flxlst[13] + flxlst[16] + flxlst[17] + flxlst[19] + flxlst[20])
            MB_tmp = -999
            plt_title = 'MODFLOW water flux balance for the whole catchment\nsum of fluxes: %.4f' % MB_tmp
            flxlbl = []
            flxlst = []
            for l in range(nlay):
                flxlst.append(cbc_RCH[:,:,:,l].sum()/sum(perlen)/ncell)
                flxlst.append(cbc_STO[:,:,:,l].sum()/sum(perlen)/ncell)
                flxlst.append(cbc_DRN[:,:,:,l].sum()/sum(perlen)/ncell)
                flxlst.append(cbc_WEL[:,:,:,l].sum()/sum(perlen)/ncell)
                for x in range(len(index_cbc_uzf)):
                    flxlbl.append('GW_' + cbc_uzf_nam[index_cbc_uzf[x]] + '_L' + str(l+1))
                for x in range(len(index_cbc)):
                    flxlbl.append('GW_' + cbc_nam[index_cbc[x]] + '_L' + str(l+1))
            colors_flx = CreateColors.main(hi=0, hf=180, numbcolors = len(flxlbl))
            MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title, fluxmax = cbcmax, fluxmin = cbcmin)
            del flxlst

    # exporting MM to ASCII files and plots at observations cells
    if obsCHECK == 1 and os.path.exists(h5_MM_fn):
        h5_MM = h5py.File(h5_MM_fn, 'r')
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
            h_satflow = MM_SATFLOW.run(cbc_RCH[:,i,j,0], float(obs.get(obs.keys()[o])['hi']),float(obs.get(obs.keys()[o])['h0']),float(obs.get(obs.keys()[o])['RC']),float(obs.get(obs.keys()[o])['STO']))
            # export ASCII file at piezometers location
            #TODO extract heads at piezo location and not center of cell
            if obs_h[o] != []:
                obs_h_tmp = obs_h[o][0,:]
            else:
                obs_h_tmp = []
            if obs_S[o] != []:
                obs_S_tmp = obs_h[o][0,:]
            else:
                obs_S_tmp = []
            MM_PROCESS.ExportResultsMM(i, j, inputDate, _nslmax, MM, index, MM_S, index_S, cbc_RCH[:,i,j,0], -cbc_WEL[:,i,j,0], h_satflow, h_MF_m[:,i,j,l], obs_h_tmp, obs_S_tmp, outFileExport[o], obs.keys()[o])
            outFileExport[o].close()
            MM_PROCESS.ExportResultsPEST(i, j, inputDate, _nslmax, MM[:,index.get('iHEADScorr')], obs_h_tmp, obs_S_tmp, outPESTheads, outPESTsm, obs.keys()[o], MM_S[:,:,index_S.get('iSu_pc')])
            # plot
            if plot_out == 1:
                plt_title = obs.keys()[o]
                # index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iEXF':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSs':13, 'iETg':14, 'iETu':15, 'idSu':16, 'iHEADScorr':17, 'idtwt':18}
                # index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'iRexf':4, 'idSu':5, 'iSu':6, 'iSAT':7, 'iMB_l':8}
                plt_export_fn = os.path.join(MM_ws, '_plt_0'+ obs.keys()[o] + '.png')
                # def allPLOT(DateInput, P, PET, PE, Pe, dPOND, POND, Ro, Eu, Tu, Eg, Tg, S, dS, Spc, Rp, EXF, ETg, Es, MB, MB_l, dtwt, SAT, R, h_MF, h_MF_corr, h_SF, hobs, Sobs, Sm, Sr, hnoflo, plt_export_fn, plt_title, colors_nsl, hmax, hmin):
                MMplot.allPLOT(
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
                hnoflo,
                plt_export_fn,
                plt_title,
                colors_nsl,
                hmin[o] + hdiff,
                hmin[o]
                )
                # plot water balance at each obs. cell
                plt_export_fn = os.path.join(MM_ws, '_plt_0'+ obs.keys()[o] + '_UNSATbalance.png')
                # flxlbl = ['RF', 'I', 'EXF', 'dSs', 'Ro', 'Es', 'dSu']
                # flxlbl1 = ['Eu', 'Tu']
                # flxlbl1a = ['Rp','Rexf']
                # flxlbl2 = ['ETu', 'Eg', 'Tg', 'ETg']
                flxlst =[MM[:,index.get('iRF')].sum()/sum(perlen),
                    -MM[:,index.get('iI')].sum()/sum(perlen),
                     MM[:,index.get('iEXF')].sum()/sum(perlen),
                     MM[:,index.get('idSs')].sum()/sum(perlen),
                    -MM[:,index.get('iRo')].sum()/sum(perlen),
                    -MM[:,index.get('iEs')].sum()/sum(perlen),
                     MM[:,index.get('idSu')].sum()/sum(perlen),
                    -MM_S[:,:,index_S.get('iEu')].sum()/sum(perlen),
                    -MM_S[:,:,index_S.get('iTu')].sum()/sum(perlen),
                    -MM_S[:,-1,index_S.get('iRp')].sum()/sum(perlen),
                    -MM[:,index.get('iETu')].sum()/sum(perlen),
                    -MM[:,index.get('iEg')].sum()/sum(perlen),
                    -MM[:,index.get('iTg')].sum()/sum(perlen),
                    -MM[:,index.get('iETg')].sum()/sum(perlen)]
                #MB_tmp = sum(flxlst) - (flxlst[11] + flxlst[13])
                MB_tmp = -999
                if isinstance(h5_MF_fn, str):
                    plt_export_fn = os.path.join(MM_ws, '_plt_0'+ obs.keys()[o] + '_UNSATandGWbalances.png')
                    for l in range(nlay):
                        flxlst.append(cbc_RCH[:,i,j,l].sum()/sum(perlen))
                        flxlst.append(cbc_STO[:,i,j,l].sum()/sum(perlen))
                        flxlst.append(cbc_DRN[:,i,j,l].sum()/sum(perlen))
                        flxlst.append(cbc_WEL[:,i,j,l].sum()/sum(perlen))
                    #MB_tmp -= flxlst[16] + flxlst[17] + flxlst[19] + flxlst[20]
                    MB_tmp = -999
                MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title + '\nsum of fluxes: %.4f' % MB_tmp, fluxmax = flxmax, fluxmin = flxmin)
                del flxlst
        del h_satflow, MM, MM_S
        h5_MM.close()
        # output for PEST
        outPESTheads.close()
        outPESTsm.close()
        del obs, obs_h, obs_S
    del MM_PROCESS, cbc_STO, cbc_RCH, cbc_WEL

    # plot MM output
    if plot_out == 1 and os.path.exists(h5_MM_fn):
        #h5_MM = h5py.File(h5_MM_fn, 'r')
        flxlbl = ['RF', 'RFe', 'I', 'EXF', 'dSs', 'Ro', 'Es', 'Eg', 'Tg', 'ETg', 'ETu', 'dSu']
        TSlst = []
        TS = 0
        while TS < sum(perlen):
            TSlst.append(TS)
            TS += plot_freq
        TSlst.append(sum(perlen)-1)
        # plot at time interval
        for i in flxlbl:
            # plot average for the whole simulated period
            i1 = 'i'+i
            h5_MM = h5py.File(h5_MM_fn, 'r')
            MM = h5_MM['MM'][:,:,:,index.get(i1)]
            h5_MM.close()
            V = [np.ma.masked_values(MM, hnoflo, atol = 0.09)]
            V[0] = np.add.accumulate(V[0], axis = 0)[sum(perlen)-1,:,:]/sum(perlen)
            Vmax = np.nanmax(V[0]) #float(np.ceil(np.nanmax(V)))
            Vmin = np.nanmin(V[0]) #float(np.floor(np.nanmin(V)))
            if i == 'dSu':
                i_lbl = '$\Delta$Su' #'$\Delta$S$_{u}'
            elif i == 'dSs':
                i_lbl = '$\Delta$Ss'
            else:
                i_lbl = i
            if Vmax!=0.0 or Vmin!=0.0:
                MMplot.plotLAYER(TS = 99998, ncol = ncol, nrow = nrow, nlay = nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i_lbl + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_' + 'average_' + i), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrsMM, ntick = ntick)
            del V
            for TS in TSlst:
                V = [np.ma.masked_values(MM[TS,:,:], hnoflo, atol = 0.09)]
                if Vmax!=0.0 or Vmin!=0.0:
                    MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_'+i), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrsMM, ntick = ntick)
            del V, MM
        flxlbl = ['Eu', 'Tu','Rp']
        for i in flxlbl:
            # plot average for the whole simulated period
            i1 = 'i'+i
            if i != 'Rp':
                V = [np.zeros([nrow,ncol])]
                for l in range(_nslmax):
                    h5_MM = h5py.File(h5_MM_fn, 'r')
                    MM = h5_MM['MM_S'][:,:,:,l,index_S.get(i1)]
                    h5_MM.close()
                    V1 = np.ma.masked_values(MM, hnoflo, atol = 0.09)
                    V1 = np.add.accumulate(V1, axis = 0)[sum(perlen)-1,:,:]/sum(perlen)
                    V += V1
                del V1
            else:
                h5_MM = h5py.File(h5_MM_fn, 'r')
                MM = h5_MM['MM_S'][:,:,:,-1,index_S.get(i1)]
                h5_MM.close()
                V = [np.ma.masked_values(MM, hnoflo, atol = 0.09)]
                V[0] = np.add.accumulate(V[0], axis = 0)[sum(perlen)-1,:,:]/sum(perlen)
                i = 'Rp_botlayer'
            Vmax = np.nanmax(V[0]) #float(np.ceil(np.nanmax(V)))
            Vmin = np.nanmin(V[0]) #float(np.floor(np.nanmin(V)))
            if Vmax!=0.0 or Vmin!=0.0:
                MMplot.plotLAYER(TS = 99998, ncol = ncol, nrow = nrow, nlay = nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_'+ 'average_' + i), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrsMM, ntick = ntick)
            del V
            for TS in TSlst:
                h5_MM = h5py.File(h5_MM_fn, 'r')
                MM = h5_MM['MM_S'][TS,:,:,:,index_S.get(i1)]
                h5_MM.close()
                if i1 != 'iRp':
                    V = [np.ma.masked_values(MM[:,:,:], hnoflo, atol = 0.09)]
                    V[0] = np.add.accumulate(V[0], axis = 2)[:,:,_nslmax-1]
                else:
                    V = [np.ma.masked_values(MM[:,:,-1], hnoflo, atol = 0.09)]
                if Vmax!=0.0 or Vmin!=0.0:
                    MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = 1, V = V,  cmap = plt.cm.Blues, CBlabel = (i + ' (mm/day)'), msg = 'no flux', plt_title = ('MM_'+i), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrsMM, ntick = ntick)
            del V, MM
        if h_diff_surf != None:
            flxlbl = ['h_diff_surf']
            Vmax = np.nanmax(h_diff_surf) #float(np.ceil(np.nanmax(V)))
            Vmin = np.nanmin(h_diff_surf) #float(np.floor(np.nanmin(V)))
            if Vmax!=0.0 or Vmin!=0.0:
                MMplot.plotLAYER(TS = 99998, ncol = ncol, nrow = nrow, nlay = 1, nplot = 1, V = h_diff_surf,  cmap = plt.cm.Blues, CBlabel = ('(m)'), msg = 'no value', plt_title = ('_HEADSmaxdiff_ConvLoop'), MM_ws = MM_ws, interval_type = 'arange', interval_diff = (Vmax - Vmin)/nrangeMM, Vmax = Vmax, Vmin = Vmin, contours = ctrsMM, ntick = ntick)
        del TSlst, flxlbl, i, i1, h_diff_surf

    # plot MF output
    if plot_out == 1 and isinstance(h5_MF_fn, str):
        # plot heads (grid + contours), DRN, etc... at specified TS
        TSlst = []
        TS = 0
        while TS < len(h_MF_m):
            TSlst.append(TS)
            TS += plot_freq
        TSlst.append(len(h_MF_m)-1)
        for TS in TSlst:
            # plot heads [m]
            V=[]
            for L in range(nlay):
                V.append(h_MF_m[TS,:,:,L])
            MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation (m)', msg = 'DRY', plt_title = 'HEADS', MM_ws = MM_ws, interval_type = 'arange', interval_diff = (hmaxMF - hminMF)/nrangeMF, Vmax = hmaxMF, Vmin = hminMF, ntick = ntick)
            del V
            # plot GW drainage [mm]
            V = []
            for L in range(nlay):
                V.append(-cbc_DRN[TS,:,:,L])
            MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage (mm/day)', msg = '- no drainage', plt_title = 'DRN', MM_ws = MM_ws, interval_type = 'linspace', interval_num = 5, Vmin = DRNmin, contours = ctrsMF, Vmax = DRNmax, fmt='%.3G', ntick = ntick)
            del V
        del TSlst

    del h_MF_m, cbc_DRN
    del top_array_m, gridSOIL, inputDate
    del hmaxMF, hminMF, hmax, hmin, hdiff, DRNmax, DRNmin, cbcmax, cbcmin

    timeendExport = mpl.dates.datestr2num(mpl.dates.datetime.datetime.today().isoformat())
    durationExport=(timeendExport-timestartExport)

    # final report of successful run
    print ('\n##############\nMARMITES executed successfully!\n%s\n') % mpl.dates.num2date(timestart).isoformat()[:19]
    print ('%s stress periods\n%s days\n%sx%s cells (rows x cols)') % (str(int(nper)),str(int(sum(perlen))),str(nrow),str(ncol))
    print ('\nMARMITES run time: %s minute(s) and %.1f second(s)') % (str(int(duration*24.0*60.0)), (duration*24.0*60.0-int(duration*24.0*60.0))*60)
    print ('MODFLOW run time: %s minute(s) and %.1f second(s)') % (str(int(durationMF*24.0*60.0)), (durationMF*24.0*60.0-int(durationMF*24.0*60.0))*60)
    print ('Export run time: %s minute(s) and %.1f second(s)') % (str(int(durationExport*24.0*60.0)), (durationExport*24.0*60.0-int(durationExport*24.0*60.0))*60)
    print ('\nOutput written in folder: \n%s\n##############\n') % MM_ws

except StandardError, e:  #Exception
    print ('\n##############\nMARMITES ERROR!')
    print '\nAbnormal MM run interruption in the export phase!\nError description:'
    traceback.print_exc(file=sys.stdout)
#    traceback.print_exc(limit=1, file=sys.stdout)
    try:
        h5_MM.close()
        h5_MF.close()
    except:
        pass
    sys.exit()

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