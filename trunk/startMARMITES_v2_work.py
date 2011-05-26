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

""" See info in MARMITESunsat.py"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.2"
__date__ = "November 2010"

import pylab
import sys
import traceback
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESsurf')
import startMARMITESsurface as startMMsurf
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESunsat_v2')
import MARMITESunsat_v2 as MMunsat
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESutilities')
import MARMITESprocess as MMproc
sys.path.append(r'E:\00code\MARMITES\trunk\ppMF_FloPy')
import ppMODFLOW_flopy as ppMF
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESutilities\MARMITESplot')
import MARMITESplot as MMplot
import CreateColors
import StringIO
import h5py

#####################################

# workspace (ws) definition
timestart = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
print '\n##############\nMARMITES started!\n%s\n##############' % pylab.num2date(timestart).isoformat()[:19]

# read input file (called _input.ini in the MARMITES workspace
# the first character on the first line has to be the character used to comment
# the file can contain any comments as the user wish, but the sequence of the input has to be respected
MM_ws = r'E:\00code\00ws\zz_TESTS\MARMITESv2_r13c6l2'
MM_fn = '_inputMM.ini'

inputFile = MMproc.readFile(MM_ws,MM_fn)

l=0
try:
    # # ECHO ON/OFF (1 - interpreter verbose BUT NO report, 0 - NO interpreter verbose BUT report)
    verbose = int(inputFile[l].strip())
    l += 1
    # output plot (1 is YES, 0 is NO)
    plot_out  = int(inputFile[l].strip())
    l += 1
    #run MARMITESunsat  (1 is YES, 0 is NO)
    MMunsat_yn = int(inputFile[l].strip())
    l += 1
    #run MARMITESsurface  (1 is YES, 0 is NO)
    MMsurf_yn = int(inputFile[l].strip())
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
    plot_freq =  int(inputFile[l].strip())
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
    gridPONDhmax_fn =  inputFile[l].strip()
    l += 1
    gridPONDw_fn =  inputFile[l].strip()
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
    report_fn = os.path.join(MM_ws,'_report_MMrun.txt')
    print '\nECHO OFF (no screen output).\nSee MM run report in file:\n%s\n' % report_fn
    s = sys.stdout
    report = open(report_fn, 'w')
    sys.stdout = report
    print '\n##############\nMARMITES started!\n%s\n##############' % pylab.num2date(timestart).isoformat()[:19]
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

# MAIN CODE #
try:
    # #############################
    # ###  READ MODFLOW CONFIG ####
    # #############################

    print'\n##############'
    print 'Importing MODFLOW configuration file'
    nrow, ncol, delr, delc, nlay, nper, perlen, nstp, hnoflo, hdry, laytyp, lenuni, itmuni = ppMF.ppMFini(MF_ws, MF_ini_fn, out = 'MM')

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
        # TODO if lenuni<>2 apply conversion factor to delr, delc, etc...
    if laytyp[0]==0:
        print 'FATAL ERROR!\nThe first layer cannot be confined type!\nChange your parameter laytyp in the MODFLOW lpf package.\n(see USGS Open-File Report 00-92)'
        sys.exit()
    if itmuni <> 4:
        print 'FATAL ERROR! Time unit is not in days!'
        sys.exit()

    # #############################
    # ###  MARMITES SURFACES  #####
    # #############################

    print'\n##############'
    print 'MARMITESsurf RUN'

    if MMsurf_yn>0:
        outMMsurf_fn = startMMsurf.MMsurf(MMsurf_ws, inputFile_TS_fn, inputFile_PAR_fn, outputFILE_fn, MM_ws, outMMsurf_fn)

    inputFile = MMproc.readFile(MM_ws,outMMsurf_fn)

    l=0
    TRANS_vdw = []
    Zr = []
    k_Tu_slp = []
    k_Tu_inter = []
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
            k_Tu_slp.append(float(line[v]))
        l += 1
        line = inputFile[l].split()
        for v in range(NVEG):
            k_Tu_inter.append(float(line[v]))
        l += 1
        line = inputFile[l].split()
        for v in range(NSOIL):
            TRANS_sdw.append(int(line[v]))
    except:
        print '\nType error in file [' + inputFile_fn + ']'
        sys.exit()
    del inputFile

    # #############################
    # ### 1st MODFLOW RUN with initial user-input recharge
    # #############################

    # if required by user, compute nper, perlen,etc based on RF analysis in the METEO zones
    if isinstance(nper, str):
        try:
            MFtime_fn = nper.split()[0].strip()
            perlenmax = int(nper.split()[1].strip())
        except:
            print '\nError in nper format of the MODFLOW ini file!\n'
            sys.exit()
        nper, perlen, nstp, inputZON_TS_RF_fn, inputZON_TS_PET_fn, inputZON_TS_RFe_fn, inputZON_TS_PE_fn, inputZON_TS_E0_fn = ppMF.ppMFtime(MM_ws, MF_ws, MFtime_fn, perlenmax, inputDate_fn, inputZON_TS_RF_fn, inputZON_TS_PET_fn, inputZON_TS_RFe_fn, inputZON_TS_PE_fn, inputZON_TS_E0_fn, NMETEO, NVEG, NSOIL)
    else:
        MFtime_fn = None

    durationMF = 0.0
    timestartMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
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
    if MMunsat_yn > 0:
        top_array, ibound, h5_MF_fn = ppMF.ppMF(MF_ws, MF_ini_fn, rch_input = (h5_MM_fn, 'rch'), wel_input = (h5_MM_fn, 'ETg'), report = report, verbose = verbose, chunks = chunks)
    else:
        top_array, ibound, h5_MF_fn = ppMF.ppMF(MF_ws, MF_ini_fn, report = report, verbose = verbose, chunks = chunks)

    h5_MF = h5py.File(h5_MF_fn)

    # cbc format is: (kstp), kper, textprocess, ncol, nrow, nlay
    cbc_nam = []
    for c in h5_MF['cbc_nam']:
        cbc_nam.append(c.strip())
    iDRN = cbc_nam.index('DRAINS')
    iSTO = cbc_nam.index('STORAGE')
    iRCH = cbc_nam.index('RECHARGE')
    iWEL = cbc_nam.index('WELLS')

    timeendMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
    durationMF +=  timeendMF-timestartMF

    # ####   SUMMARY OF MODFLOW READINGS   ####
    # active cells in layer 1_____________________ibound[0]
    # elevation___________________________________top
    # heads in layer n____________________________heads[timestep, row, col, layer]
    # aquifer type of layer 1_____________________laytyp[0]
    # numero total de time step___________________sum(nstp)
    # code for dry cell___________________________hdry
    # cell size___________________________________delr[0]

    print'\n##############'
    print 'MARMITESunsat initialization'

    # MARMITES INITIALIZATION
    MM_PROCESS = MMproc.PROCESS(MM_ws               = MM_ws,
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

    # READ input ESRI ASCII rasters # missing gridIRR_fn
    print "\nImporting ESRI ASCII files for MARMITES initialization..."
    gridMETEO = MM_PROCESS.inputEsriAscii(grid_fn = gridMETEO_fn, datatype = int)

    gridSOIL = MM_PROCESS.inputEsriAscii(grid_fn = gridSOIL_fn, datatype = int)

    gridSOILthick = MM_PROCESS.inputEsriAscii(grid_fn = gridSOILthick_fn,
     datatype = float)

    gridPONDhmax = MM_PROCESS.inputEsriAscii(grid_fn = gridPONDhmax_fn,
     datatype = float)

    gridPONDw = MM_PROCESS.inputEsriAscii(grid_fn = gridPONDw_fn,
     datatype = float)

    # TODO compute correction factor for pionding and convert to mm
    # TODO check unit consistency (input should be m and be converted to mm)

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
                                    inputZON_TS_E0_fn        = inputZON_TS_E0_fn,
                                    MFtime_fn                = MFtime_fn
     ) # IRR_fn

    # SOIL PARAMETERS
    _nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _Si, _Ks = MM_PROCESS.inputSoilParam(SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
    _nslmax = max(_nsl)

    # indexes of the HDF5 output arrays
    index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iPOND':4, 'iRo':5, 'iSEEPAGE':6, 'iEs':7, 'iMB':8, 'iINTER':9, 'iE0':10, 'iEg':11, 'iTg':12, 'iR':13, 'iRn':14, 'idPOND':15, 'iETg':16}
    index_S = {'iEu':0, 'iTu':1,'iSpc':2, 'iRp':3, 'idS':4, 'iS':5}

    duration = 0.0
    if MMunsat_yn > 0:

        h_pSP = 0
        LOOP = 0
        LOOPlst = [LOOP]
        h_diff = [1000]
        h_diff_log = [1]
        loopdry = 0

        plt_ConvLoop_fn = os.path.join(MM_ws, '_plt_ConvLoop_MM_MF.png')

        # #############################
        # ###  CONVERGENCE LOOP   #####
        # #############################

        while abs(h_diff[LOOP]) > convcrit:
            print '\n##############\nCONVERGENCE LOOP #%s\n##############' % str(LOOP)
            h_MFsum = 0.0
            # ###########################
            # ###  MARMITES INPUT #######
            # ###########################

            print'\n##############'
            print 'MARMITESunsat RUN'

            # SOIL PARAMETERS
            _nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _Si, _Ks = MM_PROCESS.inputSoilParam(SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
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
                h5_MM.create_dataset(name = 'rch', shape = (sum(nstp),nrow,ncol), dtype = np.float, chunks = (1,nrow,ncol),  compression = 'gzip', compression_opts = 5, shuffle = True)
                h5_MM.create_dataset(name = 'ETg', shape = (sum(nstp),nrow,ncol), dtype = np.float, chunks = (1,nrow,ncol),  compression = 'gzip', compression_opts = 5, shuffle = True)
            else:
                h5_MM.create_dataset(name = 'rch', shape = (sum(nstp),nrow,ncol), dtype = np.float)
                h5_MM.create_dataset(name = 'ETg', shape = (sum(nstp),nrow,ncol), dtype = np.float)
            # ###############
            # # main loop: calculation of soil water balance in each cell-grid for each time step inside each stress period
    #        t0=0
            print '\nComputing...'

            # initial values of SP
            Si_tmp_array = np.zeros([nrow,ncol,_nslmax], dtype=float)
            Rpi_tmp_array = np.zeros([nrow,ncol,_nslmax], dtype=float)
            PONDi_tmp_array = np.zeros([nrow,ncol], dtype=float)
            DRNi_tmp_array = np.zeros([nrow,ncol], dtype=float)
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
                MM_rch = np.zeros([nstp[n],nrow,ncol], dtype=float)
                MM_ETg = np.zeros([nstp[n],nrow,ncol], dtype=float)
                h_MF_per = h5_MF['heads'][tstart_MF:tend_MF,:,:,0]
                drn_MF_per = h5_MF['cbc'][tstart_MF:tend_MF,iDRN,:,:,0]
                ncell = 0
                # loop into the grid
                for i in range(nrow):
                    for j in range(ncol):
                        SOILzone_tmp = gridSOIL[i,j]-1
                        METEOzone_tmp = gridMETEO[i,j]-1
                        if ibound[i,j,0]<>0.0:
                            ncell += 1
                            nsl_tmp   = _nsl[SOILzone_tmp]
                            st_tmp    = _st[SOILzone_tmp]
                            slprop_tmp= _slprop[SOILzone_tmp]
                            Sm_tmp    = _Sm[SOILzone_tmp]
                            Sfc_tmp   = _Sfc[SOILzone_tmp]
                            Sr_tmp    = _Sr[SOILzone_tmp]
                            Si_tmp = []
                            if n==0:
                                for l in range(nsl_tmp):
                                    Si_tmp.append(_Si[SOILzone_tmp][l])
                                PONDi  = 0.0
                                DRNi = 0.0
                            else:
                                Si_tmp = Si_tmp_array[i,j,:]
                                PONDi  = PONDi_tmp_array[i,j]
                                DRNi  = DRNi_tmp_array[i,j]
                            Rpi_tmp = Rpi_tmp_array[i,j,:]
                            Ks_tmp    = _Ks[SOILzone_tmp]
                            PONDm_tmp = 1000*1.12*gridPONDhmax[i,j]*gridPONDw[i,j]/delr[j]
                            PONDratio = 1.12*gridPONDw[i,j]/delr[j]
                            D_tmp = gridSOILthick[i,j]
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
                            h_MF_tmp = h_MF_per[0:nstp[n],i,j]
                            drn_MF_tmp = -conv_fact*drn_MF_per[0:nstp[n],i,j]/(delr[j]*delc[i])
                            # for the while loop
                            h_MF_L0_m = np.ma.masked_values(h_MF_tmp, hdry, atol = 1E+25)
                            h_MF_L0_m_avg = np.nansum(h_MF_L0_m)/nstp[n]
                            if isinstance(h_MF_L0_m_avg, float):
                                h_MFsum += h_MF_L0_m_avg
                            # cal functions for reservoirs calculations
                            MM1_temp, MM2_temp = MM_UNSAT.run(
                                                         i, j, n,
                                                         nsl     = nsl_tmp,
                                                         st      = st_tmp,
                                                         slprop  = slprop_tmp,
                                                         Sm      = Sm_tmp,
                                                         Sfc     = Sfc_tmp,
                                                         Sr      = Sr_tmp,
                                                         Si      = Si_tmp,
                                                         PONDi   = PONDi,
                                                         Rpi     = Rpi_tmp,
                                                         D       = D_tmp,
                                                         Ks      = Ks_tmp,
                                                         PONDm   = PONDm_tmp,
                                                         PONDratio = PONDratio,
                                                         ELEV    = top_array[i,j],
                                                         HEADS   = h_MF_tmp,
                                                         DRN     = drn_MF_tmp,
                                                         DRNi    = DRNi,
                                                         RF      = RFzonesTS[METEOzone_tmp][tstart_MF:tend_MF],
                                                         E0      = E0zonesTS[METEOzone_tmp][tstart_MF:tend_MF],
                                                         PETveg  = PETvegzonesTS_tmp,
                                                         RFeveg  = RFevegzonesTS_tmp,
                                                         PEsoil  = PEsoilzonesTS_tmp,
                                                         VEGarea = VEGarea_tmp,
                                                         Zr      = Zr,
                                                         nstp  = nstp[n],
                                                         hdry    = hdry,
                                                         k_Tu_slp = k_Tu_slp,
                                                         k_Tu_inter = k_Tu_inter)
                            if perlen[n]>1:
                                for p in range(perlen[n]):
                                    for k in range(len(index)):
                                        MM[p,i,j,k] = MM1_temp[:,k]
                                    for k in range(len(index_S)):
                                        for l in range(nsl_tmp):
                                            MM_S[p,i,j,l,k] = MM2_temp[:,l,k]
                                del MM1_temp, MM2_temp
                            else:
                                for k in range(len(index)):
                                    MM[:,i,j,k] = MM1_temp[:,k]
                                for k in range(len(index_S)):
                                    for l in range(nsl_tmp):
                                        MM_S[:,i,j,l,k] = MM2_temp[:,l,k]
                                del MM1_temp, MM2_temp
                            divfact = float(perlen[n])
                            MM_rch[:,i,j] = MM[:,i,j,index.get('iR')].sum()/(divfact*1000)
                            MM_ETg[:,i,j] = (delr[j]*delc[i])*MM[:,i,j,index.get('iETg')].sum()/(divfact*1000)
                            # setting initial conditions for the next SP
                            for l in range(nsl_tmp):
                                Si_tmp_array[i,j,l] = MM_S[nstp[n]-1,i,j,l,index_S.get('iSpc')]
                                Rpi_tmp_array[i,j,l] = MM_S[nstp[n]-1,i,j,l,index_S.get('iRp')]
                                PONDi_tmp_array[i,j] = MM[nstp[n]-1,i,j,index.get('iPOND')]
                                DRNi_tmp_array[i,j]  = drn_MF_tmp[-1]
                        else:
                            if perlen[n]>1:
                                for p in range(perlen[n]):
                                    for k in range(len(index)):
                                        MM[p,i,j,k] = hnoflo
                                    for k in range(len(index_S)):
                                        for l in range(nsl_tmp):
                                            MM_S[p,i,j,l,k] = hnoflo
                            else:
                                for k in range(len(index)):
                                        MM[:,i,j,k] = hnoflo
                                for k in range(len(index_S)):
                                    for l in range(nsl_tmp):
                                        MM_S[:,i,j,l,k] = hnoflo
                            MM_rch[:,i,j] = hnoflo
                            MM_ETg[:,i,j] = hnoflo
                    h5_MM['MM'][tstart_MM:tend_MM,:,:,:] = MM[:,:,:,:]
                    h5_MM['MM_S'][tstart_MM:tend_MM,:,:,:,:] = MM_S[:,:,:,:,:]
                    h5_MM['rch'][tstart_MF:tend_MF,:,:] = MM_rch
                    h5_MM['ETg'][tstart_MF:tend_MF,:,:] = MM_ETg
                del h_MF_per, drn_MF_per
                del MM, MM_S, MM_rch, MM_ETg
    #            t0 += int(nstp[n])

             #   print '\nSTRESS PERIOD %i/%i DONE!' % (n+1, nper)
            h5_MM.close()

            # #############################
            # ###  PRINT CONVERG. PLOT ###
            # #############################
            h_MFsum = h_MFsum/ncell/nper
            h_diff.append(h_MFsum - h_pSP)
            LOOP += 1
            LOOPlst.append(LOOP)
            h_pSP = h_MFsum
            if pylab.absolute(h_diff[LOOP])>0.0:
                h_diff_log.append(pylab.log10(pylab.absolute(h_diff[LOOP])))
            else:
                h_diff_log.append(pylab.log10(convcrit))
            if LOOP <2:
                print "\nInitial average heads:\n%.3f m" % h_diff[LOOP]
            else:
                print "\nHeads diff. from previous conv. loop:\n%.3f m" % h_diff[LOOP]
            if h_MFsum == 0.0:
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

            # #############################
            # ### MODFLOW RUN with MM-computed recharge
            # #############################
            h5_MF.close()
            durationMF = 0.0
            timestartMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
            print'\n##############'
            print 'MODFLOW RUN (MARMITES fluxes)'
            if verbose == 0:
                print '\n--------------'
                sys.stdout = s
                report.close()
                s = sys.stdout
                report = open(report_fn, 'a')
                sys.stdout = report
            top_array, ibound, h5_MF_fn = ppMF.ppMF(MF_ws, MF_ini_fn, rch_input = (h5_MM_fn, 'rch'), wel_input = (h5_MM_fn, 'ETg'), report = report, verbose = verbose, chunks = chunks)

            h5_MF = h5py.File(h5_MF_fn)

            timeendMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
            durationMF += (timeendMF-timestartMF)

        # #############################
        # ###  END CONVERGENCE LOOP ###
        # #############################

        # export loop plot
        print'\n##############'
        print 'Exporting plot of the convergence loop...'
        fig = plt.figure()
        if LOOP>0:
            ax1=fig.add_subplot(3,1,1)
            plt.setp(ax1.get_xticklabels(), fontsize=8)
            plt.setp(ax1.get_yticklabels(), fontsize=8)
            ax1.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.3G'))
            plt.ylabel('h_diff\n[m]', fontsize=10, horizontalalignment = 'center')
            plt.grid(True)
            plt.plot(LOOPlst[1:], h_diff[1:], linestyle='-', marker='o', markersize=5, c = 'orange', markerfacecolor='orange', markeredgecolor='red')
        if LOOP>1:
            ax2=fig.add_subplot(3,1,3, sharex = ax1)
            plt.setp(ax2.get_xticklabels(), fontsize=8)
            plt.setp(ax2.get_yticklabels(), fontsize=8)
            ax2.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.3G'))
            plt.ylabel('log(abs(h_diff))\n[log(m)]', fontsize=10, horizontalalignment = 'center')
            plt.grid(True)
            plt.xlabel('trial', fontsize=10)
            plt.plot(LOOPlst[2:], h_diff_log[2:], linestyle='-', marker='o', markersize=5, c = 'orange', markerfacecolor='orange', markeredgecolor='red')
            ax3=fig.add_subplot(3,1,2, sharex = ax1)
            plt.setp(ax3.get_xticklabels(), fontsize=8)
            plt.setp(ax3.get_yticklabels(), fontsize=8)
            plt.plot(LOOPlst[2:], h_diff[2:], linestyle='-', marker='o', markersize=5, c = 'orange', markerfacecolor='orange', markeredgecolor='red')
            ax3.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.3G'))
            plt.ylabel('h_diff\n[m]', fontsize=10, horizontalalignment = 'center')
        #        plt.ylabel.Text.position(0.5, -0.5)
            plt.grid(True)
            ax2.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%2d'))
            ax3.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%2d'))
        if LOOP>1:
            plt.xlim(1,LOOP)
            ax1.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%2d'))
            ax1.xaxis.set_ticks(LOOPlst[1:])
        plt.savefig(plt_ConvLoop_fn)
        plt.cla()
        plt.clf()
        plt.close('all')
        del fig, LOOPlst, h_diff, h_diff_log

        timeend = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
        duration += (timeend-timestart) -durationMF

        # #############################
        # ### MODFLOW RUN with MM-computed recharge
        # #############################
        h5_MF.close()
        durationMF = 0.0
        timestartMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
        print'\n##############'
        print 'MODFLOW RUN (MARMITES fluxes after conv. loop)'
        if verbose == 0:
            print '\n--------------'
            sys.stdout = s
            report.close()
            s = sys.stdout
            report = open(report_fn, 'a')
            sys.stdout = report

        top_array, ibound, h5_MF_fn = ppMF.ppMF(MF_ws, MF_ini_fn, rch_input = (h5_MM_fn, 'rch'), wel_input = (h5_MM_fn, 'ETg'), report = report, verbose = verbose, chunks = chunks)

        timeendMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
        durationMF += (timeendMF-timestartMF)

    del gridVEGarea, RFzonesTS, E0zonesTS, PETvegzonesTS, RFevegzonesTS, PEsoilzonesTS
    del gridMETEO, gridSOILthick, gridPONDhmax, gridPONDw

    # #############################
    # ###  OUTPUT EXPORT   ########
    # #############################

    timestartExport = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
    print '\n##############\nMARMITES exporting...'

    # reading MF output
    h5_MF = h5py.File(h5_MF_fn)

    h_MF = np.zeros((sum(perlen), nrow, ncol, nlay))
    cbc_MF = np.zeros((sum(perlen), len(cbc_nam), nrow, ncol, nlay))
    t = 0
    ncell = 0
    for n in range(nper):
        if perlen[n]>1:
            for x in range(perlen[n]):
                h_MF[t,:,:,:] = h5_MF['heads'][n,:,:,:]
                for i in range(nrow):
                    for j in range(ncol):
                        cbc_MF[t,:,i,j,:] = conv_fact*h5_MF['cbc'][n,:,i,j,:]/(delr[j]*delc[i])
                        if ibound[i,j,0]<>0.0:
                            ncell += 1
                t += 1
        else:
            h_MF[t,:,:,:] = h5_MF['heads'][n,:,:,:]
            for i in range(nrow):
                for j in range(ncol):
                    cbc_MF[t,:,i,j,:] = conv_fact*h5_MF['cbc'][n,:,i,j,:]/(delr[j]*delc[i])
                    if ibound[i,j,0]<>0.0:
                            ncell += 1
            t += 1
    ncell /= sum(perlen)
    h_MF_m = np.ma.masked_values(h_MF, hnoflo, atol = 0.09)
    del h_MF
    top_array_m = np.ma.masked_values(top_array, hnoflo, atol = 0.09)
    h5_MF.close()

    # reading MM output
    h5_MM = h5py.File(h5_MM_fn, 'r')
    MM = np.ma.masked_values(h5_MM['MM'], hnoflo, atol = 0.09)
    MM_S = np.ma.masked_values(h5_MM['MM_S'], hnoflo, atol = 0.09)
    h5_MM.close()

    # READ observations time series (heads and soil moisture)
    obsCHECK = 1  # 0: no obs, 1: obs
    if obsCHECK==1:
        print "\nReading observations time series (hydraulic heads and soil moisture)..."
        obs, outpathname, obs_h, obs_S = MM_PROCESS.inputObs(
                                inputObs_fn = inputObs_fn,
                                inputObsHEADS_fn = inputObsHEADS_fn,
                                inputObsSM_fn = inputObsSM_fn,
                                inputDate   = inputDate,
                                _nslmax     = _nslmax,
                                nlay        = nlay,
                                MFtime_fn   = MFtime_fn,
                                )
        # Write first output in a txt file
        outFileExport = []
        for o in range(len(obs.keys())):
            outFileExport.append(open(outpathname[o], 'w'))
            S_str=''
            Spc_str=''
            dS_str=''
            Rp_str=''
            Eu_str=''
            Tu_str=''
            Smeasout = ''
            for l in range(_nslmax):
                S_str = S_str + 'S_l' + str(l+1) + ','
                Spc_str = Spc_str + 'Spc_l' + str(l+1) + ','
                dS_str = dS_str + 'dS_l' + str(l+1) + ','
                Eu_str = Eu_str + 'Eu_l' + str(l+1) + ','
                Tu_str = Tu_str + 'Tu_l' + str(l+1) + ','
                Rp_str = Rp_str + 'Rp_l' + str(l+1) + ','
                Smeasout = Smeasout + 'Smeas_' + str(l+1) + ','
            header='Date,RF,E0,PET,PE,RFe,Inter,'+Eu_str+Tu_str+'Eg,Tg,ETg,WEL_MF,Es,'+S_str+Spc_str+dS_str+'dPOND,POND,Ro,SEEPAGE,DRN_MF,'+Rp_str+'R,Rn,R_MF,hSATFLOW,hMF,hmeas,' + Smeasout + 'MB\n'
            outFileExport[o].write(header)
        outPESTheads_fn      = 'PESTheads.dat'
        outPESTsm_fn         = 'PESTsm.dat'
        outPESTheads=open(os.path.join(MM_ws,outPESTheads_fn), 'w')
        outPESTsm=open(os.path.join(MM_ws,outPESTsm_fn), 'w')
    else:
        print "\nNo reading of observations time series (hydraulic heads and soil moisture) required."

    print '\nExporting ASCII files and plots...'
    h_MF_m = np.ma.masked_values(h_MF_m, hdry, atol = 1E+25)
    hmax = []
    hmin = []
    DRNmax = []
    DRNmin = []
    cbcmax = []
    cbcmin = []
    for L in range(nlay):
        hmax.append(np.nanmax(h_MF_m[:,:,:,:].flatten()))
        hmin.append(np.nanmin(h_MF_m[:,:,:,:].flatten()))
        DRNmax.append(np.nanmax(-cbc_MF[:,iDRN,:,:,:]).flatten())
        DRNmin.append(np.nanmin(-cbc_MF[:,iDRN,:,:,:]).flatten())
        cbcmax.append(np.nanmax(-cbc_MF[:,:,:,:,:]).flatten())
        cbcmin.append(np.nanmin(-cbc_MF[:,:,:,:,:]).flatten())
    for o in range(len(obs.keys())):
            npa_m_tmp = np.ma.masked_values(obs_h[o], hnoflo, atol = 0.09)
            hmax.append(np.nanmax(npa_m_tmp.flatten()))
            hmin.append(np.nanmin(npa_m_tmp.flatten()))
    hmax = float(np.ceil(np.nanmax(hmax)))
    hmin = float(np.floor(np.nanmin(hmin)))
    DRNmax = float(np.ceil(np.nanmax(DRNmax)))
    DRNmin = float(np.floor(np.nanmin(DRNmin)))
    cbcmax = float(np.ceil(np.nanmax(cbcmax)))
    cbcmin = float(np.floor(np.nanmin(cbcmin)))
    # plot water budget
    index_cbc = [iRCH, iSTO, iDRN, iWEL]
    if MMunsat_yn > 0:
        plt_export_fn = os.path.join(MM_ws, '_plt_UNSATandGWbudgets.png')
        plt_title = 'MARMITES and MODFLOW water flux budget for the whole catchment'
        flxlbl = ['RF', 'INTER', 'SEEPAGE', 'dS', 'dPOND', 'Ro', 'Es', 'Eu', 'Tu', 'R', 'Rn', 'Eg', 'Tg', 'ETg']
        flxlst = [MM[:,:,:,index.get('iRF')].sum()/sum(perlen)/ncell,
            -MM[:,:,:,index.get('iINTER')].sum()/sum(perlen)/ncell,
            MM[:,:,:,index.get('iSEEPAGE')].sum()/sum(perlen)/ncell,
            MM_S[:,:,:,:,index_S.get('idS')].sum()/sum(perlen)/ncell,
            MM[:,:,:,index.get('idPOND')].sum()/sum(perlen)/ncell,
            -MM[:,:,:,index.get('iRo')].sum()/sum(perlen)/ncell,
            -MM[:,:,:,index.get('iEs')].sum()/sum(perlen)/ncell,
            -MM_S[:,:,:,:,index_S.get('iEu')].sum()/sum(perlen)/ncell,
            -MM_S[:,:,:,:,index_S.get('iTu')].sum()/sum(perlen)/ncell,
            -MM[:,:,:,index.get('iR')].sum()/sum(perlen)/ncell,
            -MM[:,:,:,index.get('iRn')].sum()/sum(perlen)/ncell,
            -MM[:,:,:,index.get('iEg')].sum()/sum(perlen)/ncell,
            -MM[:,:,:,index.get('iTg')].sum()/sum(perlen)/ncell,
            MM[:,:,:,index.get('iETg')].sum()/sum(perlen)/ncell]
    else:
        plt_export_fn = os.path.join(MM_ws, '_plt_GWbudgets.png')
        plt_title = 'MODFLOW water flux budget for the whole catchment'
        flxlbl = []
        flxlst = []
    for l in range(nlay):
        for x in range(len(index_cbc)):
            flxlst.append(cbc_MF[:,index_cbc[x],:,:,l].sum()/sum(perlen)/ncell)
            flxlbl.append('GW_' + cbc_nam[index_cbc[x]] + '_L' + str(l+1))
    colors_flx = CreateColors.main(hi=0, hf=180, numbcolors = len(flxlbl))
    MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title, fluxmax = cbcmax, fluxmin = cbcmin)
    del flxlst

    # exporting ASCII files and plots at observations cells
    if obsCHECK == 1:
        colors_nsl = CreateColors.main(hi=00, hf=180, numbcolors = (_nslmax+1))
        for o in range(len(obs.keys())):
            i = obs.get(obs.keys()[o])['i']
            j = obs.get(obs.keys()[o])['j']
            l = obs.get(obs.keys()[o])['lay']
            # SATFLOW
            if MMunsat_yn > 0:
                h_satflow = MM_SATFLOW.run(MM[:,i,j,index.get('iR')], float(obs.get(obs.keys()[o])['hi']),float(obs.get(obs.keys()[o])['h0']),float(obs.get(obs.keys()[o])['RC']),float(obs.get(obs.keys()[o])['STO']))
            # export ASCII file at piezometers location
            #TODO extract heads at piezo location and not center of cell
            if MMunsat_yn > 0:
                MM_PROCESS.ExportResultsMM(i, j, inputDate, _nslmax, MM, index, MM_S, index_S, -cbc_MF[:,iDRN,i,j,0], cbc_MF[:,iRCH,i,j,0], -cbc_MF[:,iWEL,i,j,0], h_satflow, h_MF_m[:,i,j,l], obs_h[o][0,:], obs_S[o], outFileExport[o], obs.keys()[o])
            outFileExport[o].close()
            MM_PROCESS.ExportResultsPEST(i, j, inputDate, _nslmax, h_MF_m[:,i,j,l], obs_h[o][0,:], obs_S[o], outPESTheads, outPESTsm, obs.keys()[o], MM_S, index_S)
            # plot
            if plot_out == 1:
                plt_title = obs.keys()[o]
                if MMunsat_yn > 0:
                    # DateInput, P, PET, Pe, POND, dPOND, Ro, ETa, S, Rp, R, h, hmeas, Smeas, Sm, Sr):
                    # index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iPOND':4, 'iRo':5, 'iSEEPAGE':6, 'iEs':7, 'iMB':8, 'iINTER':9, 'iE0':10, 'iEg':11, 'iTg':12, 'iR':13, 'iRn':14, 'idPOND':15, 'iETg':16}
                    plt_export_fn = os.path.join(MM_ws, '_plt_'+ obs.keys()[o] + '.png')
                    # def allPLOT(DateInput, P, PET, PE, Pe, dPOND, POND, Ro, Eu, Tu, Eg, Tg, S, dS, Spc, Rp, SEEPAGE, R, Rn, Es, MB, h_MF, h_SF, hmeas, Smeas, Sm, Sr, hnoflo, plt_export_fn, plt_title, colors_nsl, hmax, hmin):
                    MMplot.allPLOT(
                    inputDate,
                    MM[:,i,j,index.get('iRF')],
                    MM[:,i,j,index.get('iPET')],
                    MM[:,i,j,index.get('iPE')],
                    MM[:,i,j,index.get('iRFe')],
                    MM[:,i,j,index.get('idPOND')],
                    MM[:,i,j,index.get('iPOND')],
                    MM[:,i,j,index.get('iRo')],
                    MM_S[:,i,j,0:_nsl[gridSOIL[i,j]-1],index_S.get('iEu')],
                    MM_S[:,i,j,0:_nsl[gridSOIL[i,j]-1],index_S.get('iTu')],
                    MM[:,i,j,index.get('iEg')],
                    MM[:,i,j,index.get('iTg')],
                    MM_S[:,i,j,0:_nsl[gridSOIL[i,j]-1],index_S.get('iS')],
                    MM_S[:,i,j,0:_nsl[gridSOIL[i,j]-1],index_S.get('idS')],
                    MM_S[:,i,j,0:_nsl[gridSOIL[i,j]-1],index_S.get('iSpc')],
                    MM_S[:,i,j,0:_nsl[gridSOIL[i,j]-1],index_S.get('iRp')],
                    MM[:,i,j,index.get('iSEEPAGE')],
                    MM[:,i,j,index.get('iR')],
                    MM[:,i,j,index.get('iRn')],
                    MM[:,i,j,index.get('iEs')],
                    MM[:,i,j,index.get('iMB')],
                    h_MF_m[:,i,j,0], h_satflow, obs_h[o][0,:], obs_S[o],
                    _Sm[gridSOIL[i,j]-1],
                    _Sr[gridSOIL[i,j]-1],
                    hnoflo,
                    plt_export_fn,
                    plt_title,
                    colors_nsl,
                    hmax,
                    hmin
                    )
                # plot water budget at each obs. cell
                if MMunsat_yn > 0:
                    plt_export_fn = os.path.join(MM_ws, '_plt_'+ obs.keys()[o] + 'UNSATandGWbudgets.png')
                    flxlst =[MM[:,i,j,index.get('iRF')].sum()/sum(perlen),
                        -MM[:,i,j,index.get('iINTER')].sum()/sum(perlen),
                         MM[:,i,j,index.get('iSEEPAGE')].sum()/sum(perlen),
                         MM_S[:,i,j,:,index_S.get('idS')].sum()/sum(perlen),
                         MM[:,i,j,index.get('idPOND')].sum()/sum(perlen),
                        -MM[:,i,j,index.get('iRo')].sum()/sum(perlen),
                        -MM[:,i,j,index.get('iEs')].sum()/sum(perlen),
                        -MM_S[:,i,j,:,index_S.get('iEu')].sum()/sum(perlen),
                        -MM_S[:,i,j,:,index_S.get('iTu')].sum()/sum(perlen),
                        -MM[:,i,j,index.get('iR')].sum()/sum(perlen),
                        -MM[:,i,j,index.get('iRn')].sum()/sum(perlen),
                        -MM[:,i,j,index.get('iEg')].sum()/sum(perlen),
                        -MM[:,i,j,index.get('iTg')].sum()/sum(perlen),
                         MM[:,i,j,index.get('iETg')].sum()/sum(perlen)]
                else:
                    plt_export_fn = os.path.join(MM_ws, '_plt_'+ obs.keys()[o] + 'GWbudgets.png')
                    flxlst = []
                for l in range(nlay):
                    for x in range(len(index_cbc)):
                        flxlst.append(cbc_MF[:,index_cbc[x],i,j,l].sum()/sum(perlen))
                MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title, fluxmax = cbcmax, fluxmin = cbcmin)
                del flxlst
            if MMunsat_yn > 0:
                del h_satflow
        # output for PEST
        outPESTheads.close()
        outPESTsm.close()
        del obs, obs_h, obs_S

    del MM_PROCESS

    if plot_out == 1:
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
            MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation (m)', msg = 'DRY', plt_title = 'HEADS', MM_ws = MM_ws, interval_type = 'arange', interval_diff = 0.5, Vmax = hmax, Vmin = hmin)
            del V
            # plot diff between drain elevation and heads elevation [m]
            DrnHeadsLtop = top_array_m - h_MF_m[TS,:,:,0]
            DrnHeadsLtop_m = np.ma.masked_greater(DrnHeadsLtop,0.0)
            V = [DrnHeadsLtop_m]
            diffMin = 0
            diffMax = np.nanmin(DrnHeadsLtop_m)
            MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = 1, V = V,  cmap = plt.cm.RdYlGn, CBlabel = 'diff. between DRN elev and hyd. heads elev. (m)', msg = ' - no drainage', plt_title = 'HEADSDRNdiff', MM_ws = MM_ws, interval_type = 'linspace', interval_num = 5, Vmax = diffMin, Vmin = diffMax, fmt='%.3G')
            # plot GW drainage [mm]
            V = []
            for L in range(nlay):
                V.append(-cbc_MF[TS,iDRN,:,:,L])
            MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage (mm/time step)', msg = '- no drainage', plt_title = 'DRN', MM_ws = MM_ws, interval_type = 'linspace', interval_num = 5, Vmin = DRNmin, Vmax = DRNmax, fmt='%.3G')
            del DrnHeadsLtop, DrnHeadsLtop_m, V

    del cbc_MF, h_MF_m, MM, MM_S
    del top_array_m, gridSOIL, inputDate
    del hmax, hmin, DRNmax, DRNmin, cbcmax, cbcmin

    timeendExport = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
    durationExport=(timeendExport-timestartExport)

    # final report of successful run

    print ('\n##############\nMARMITES executed successfully!')
    print ('%s stress periods\n%s days\n%sx%s cells (rows x cols)') % (str(int(nper)),str(int(sum(perlen))),str(nrow),str(ncol))
    print ('\nMARMITES run time: %s minute(s) and %.1f second(s)') % (str(int(duration*24.0*60.0)), (duration*24.0*60.0-int(duration*24.0*60.0))*60)
    print ('MODFLOW run time: %s minute(s) and %.1f second(s)') % (str(int(durationMF*24.0*60.0)), (durationMF*24.0*60.0-int(durationMF*24.0*60.0))*60)
    print ('Export run time: %s minute(s) and %.1f second(s)') % (str(int(durationExport*24.0*60.0)), (durationExport*24.0*60.0-int(durationExport*24.0*60.0))*60)
    print ('\nOutput written in folder: \n%s\n##############\n') % MM_ws

    del nrow, ncol, delr, delc, nlay, nstp, nper, hnoflo, hdry, ibound, laytyp, h5_MF_fn, top_array, lenuni

except StandardError, e:  #Exception
    print '\nERROR! Abnormal MM run interruption!\nError description:'
    traceback.print_exc(file=sys.stdout)
#    traceback.print_exc(limit=1, file=sys.stdout)

if verbose == 0:
    sys.stdout = s
    report.close()
    print 'MARMITES terminated!\n%s\n' % pylab.datetime.datetime.today().isoformat()[:19]
    del s
    try:
        h5_MM.close()
        h5_MF.close()
    except:
        pass