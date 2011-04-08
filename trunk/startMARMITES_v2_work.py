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

__author__ = "Alain P. Francés <frances08512@itc.nl>"
__version__ = "0.2"
__date__ = "November 2010"

import pylab
import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESsurf')
import startMARMITESsurface as startMMsurf
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESunsat_v2')
import MARMITESunsat_v2 as MMunsat
sys.path.append(r'E:\00code\MARMITES\trunk\ppMF_FloPy')
import ppMODFLOW_flopy as ppMF
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESplot')
import MARMITESplot as MMplot
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESsurf')
import CreateColors
#####################################

#try:

# workspace (ws) definition
timestart = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
print'\n##############'
print 'MARMITES starteeeeeeeeeeed!!!'

messagemanual="Please read the manual!\n(that by the way still doesn't exist...)"

# read input file (called _input.ini in the MARMITES workspace
# the first character on the first line has to be the character used to comment
# the file can contain any comments as the user wish, but the sequence of the input has to be respected
inputFile_fn = r'E:\00code\00ws\zz_TESTS\MARMITESv2_r13c6l2\_inputMM.ini'
inputFile = []
if os.path.exists(inputFile_fn):
    fin = open(inputFile_fn, 'r')
else:
    print "Input file doesn't exist, verify name and path!"
    sys.exit()
line = fin.readline().split()
delimChar = line[0]
try:
    for line in fin:
        line_tmp = line.split(delimChar)
        if not line_tmp == []:
            if (not line_tmp[0] == '') and (not line_tmp[0] == '\n'):
                inputFile.append(line_tmp[0])
        else:
            raise NameError('InputFileFormat')
except NameError:
    print 'Error in the input file, check format!\n%s' % (messagemanual)
    sys.exit()
except e:
    print "Unexpected error in the input file:\n", sys.exc_info()[0]
    print messagemanual
    sys.exit()
l=0
try:
    # Define MARMITES ws folders
    MM_ws = inputFile[l]
    l = l+1
    #run MARMITESsurface  1 is YES, 0 is NO
    MMsurf_yn = int(inputFile[l].strip())
    l = l+1
    # Define MARMITESsurface folder
    inputFOLDER_fn = inputFile[l].strip()
    l = l+1
    # METEO TIME SERIES file name
    inputFile_TS_fn = inputFile[l].strip()
    l = l+1
    # METEO/VEGETATION/SOIL/WATER PARAMETERS file name
    inputFile_PAR_fn = inputFile[l].strip()
    l = l+1
    # ouputprefix
    outputFILE_fn = inputFile[l].strip()
    l = l+1
    # ZONEVEGSOILfile
    MMsurf_fn = inputFile[l].strip()
    l = l+1
    # Define MODFLOW ws folders
    MF_ws = inputFile[l]
    l = l+1
    #GRID (ll means lower left)
    xllcorner = float(inputFile[l].strip())
    l = l+1
    yllcorner = float(inputFile[l].strip())
    l = l+1
    gridMETEO_fn = inputFile[l].strip()
    l = l+1
    gridSOIL_fn = inputFile[l].strip()
    l = l+1
    gridSOILthick_fn = inputFile[l].strip()
    l = l+1
    gridIRR_fn = inputFile[l].strip()
    l = l+1
    gridSUSTm_fn =  inputFile[l].strip()
    l = l+1
    SOILparam_fn = inputFile[l].strip()
    l = l+1
    IRR_fn = inputFile[l].strip()
    l = l+1
    inputObs_fn = inputFile[l].strip()
except:
    print '\nType error in the input\n%s' % (messagemanual)
    sys.exit()
fin.close()
print ('\nMARMITES workspace:\n%s\n\nMARMITESsurf workspace:\n%s\n\nMODFLOW workspace:\n%s' % (MM_ws, inputFOLDER_fn, MF_ws))

# #############################
# ###  MARMITES SURFACES  #####
# #############################

print'\n##############'
print 'MARMITESsurf running...'

if MMsurf_yn>0:
    MMsurf_fn = startMMsurf.MMsurf(inputFOLDER_fn, inputFile_TS_fn, inputFile_PAR_fn, outputFILE_fn, MM_ws)

inputFile = []
inputFile_fn = os.path.join(MM_ws,MMsurf_fn)
if os.path.exists(inputFile_fn):
    fin = open(inputFile_fn, 'r')
else:
    print "File [" + inputFile_fn + "] doesn't exist, verify name and path!"
    sys.exit()
line = fin.readline().split()
delimChar = line[0]
try:
    for line in fin:
        line_tmp = line.split(delimChar)
        if not line_tmp == []:
            if (not line_tmp[0] == '') and (not line_tmp[0] == '\n'):
                inputFile.append(line_tmp[0])
        else:
            raise NameError('InputFileFormat')
except NameError:
    print 'Error in file [' + inputFile_fn + '], check format!\n%s' % (messagemanual)
    sys.exit()
except e:
    print "Unexpected error in file [" + inputFile_fn + "]\n", sys.exc_info()[0]
    print messagemanual
    sys.exit()
l=0
TRANS_vdw = []
Zr = []
k_Tu_d = []
k_Tu_w = []
TRANS_sdw = []
try:
    NMETEO = int(inputFile[l].strip())
    l = l+1
    NVEG = int(inputFile[l].strip())
    l = l+1
    NSOIL = int(inputFile[l].strip())
    l = l +1
    inputDate_fn = str(inputFile[l].strip())
    l = l +1
    inputZON_TS_RF_fn = str(inputFile[l].strip())
    l = l +1
    inputZON_TS_PET_fn = str(inputFile[l].strip())
    l = l +1
    inputZON_TS_RFe_fn = str(inputFile[l].strip())
    l = l +1
    inputZON_TS_PE_fn = str(inputFile[l].strip())
    l = l +1
    inputZON_TS_E0_fn = str(inputFile[l].strip())
    l = l +1
    line = inputFile[l].split()
    for v in range(NVEG):
        TRANS_vdw.append(int(line[v]))
    l = l +1
    line = inputFile[l].split()
    for v in range(NVEG):
        Zr.append(float(line[v]))
    l = l +1
    line = inputFile[l].split()
    for v in range(NVEG):
        k_Tu_d.append(float(line[v]))
    l = l +1
    line = inputFile[l].split()
    for v in range(NVEG):
        k_Tu_w.append(float(line[v]))
    l = l +1
    line = inputFile[l].split()
    for v in range(NSOIL):
        k_Tu_w.append(int(line[v]))
except:
    print 'Type error in file [' + inputFile_fn + ']%s' % (messagemanual)
    sys.exit()
fin.close()

h_pSP = 0
TRY = 0
TRYlst = [TRY]
h_diff = [10]
h_diff_log = [1]
convcrit = 0.01
ccnum = 5 # convergence cycle number

plotCONVERGENCE_export_fn = os.path.join(MM_ws, '00conv.png')
duration = 0.0
durationMF = 0.0

# #############################
# ###  CONVERGENCE LOOP   #####
# #############################

while abs(h_diff[TRY]) > convcrit:
    h_MFsum = 0

    # #############################
    # ###  MODFLOW FILES   #####
    # #############################

    timestartMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
    print'\n##############'
    print 'MODFLOW initialization'
    SP_d, nrow, ncol, delr, delc, nlay, perlen, nper, top, hnoflo, hdry, ibound, AqType, h_MF, cbc, cbc_nam_tmp, top_array, inputFileMF_fn, lenuni = ppMF.ppMF(MF_ws, MM_ws)

    h_MF_m = np.ma.masked_values(h_MF, hnoflo, atol = 0.09)
    top_array_m = np.ma.masked_values(top_array, hnoflo, atol = 0.09)
    cbc_nam = []
    for c in cbc_nam_tmp:
        cbc_nam.append(c.strip())
    iDRN = cbc_nam.index('DRAINS')
    iSTO = cbc_nam.index('STORAGE')
    # convert cbc from volume to mm
    if lenuni == 1:
        conv_fact = 304.8
    elif lenuni == 2:
        conv_fact = 1000.0
    elif lenuni == 3:
        conv_fact = 10.0
    else:
        print 'ERROR! Define the length unit in the MODFLOW ini file!\n (see USGS Open-File Report 00-92)'
        sys.exit()

    rch_fn = 'rch_per'
    fin = open(inputFileMF_fn, 'r')
    num_line = 0
    inputFile=[]
    for line in fin:
        inputFile.append(line)
        num_line = num_line+1
    fin = open(inputFileMF_fn, 'w')
    for l in range(num_line):
        if l == num_line-1:
            fin.write(rch_fn)
        else:
            fin.write(inputFile[l])
    fin.close()

    timeendMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
    durationMF = durationMF + (timeendMF-timestartMF)

    # ####   SUMMARY OF MODFLOW READINGS   ####
    # active cells in layer 1_____________________ibound[0]
    # elevation___________________________________top
    # heads in layer n____________________________heads[timestep, row, col, layer]
    # aquifer type of layer 1_____________________AqType
    # numero total de time step___________________sum(perlen)
    # code for dry cell___________________________hdry
    # cell size___________________________________delr[0]
    #TODO use the PEST utilities for time extrapolation and use time step > 1 day in MODFLOW

    # ###########################
    # ###  MARMITES INPUT #######
    # ###########################

    print'\n##############'
    print 'MARMITESunsat initialization...'

    # MARMITES INITIALIZATION
    MM_PROCESS = MMunsat.PROCESS(MM_ws                  = MM_ws,
                            MF_ws                    = MF_ws,
                            nrow                     = nrow,
                            ncol                     = ncol,
                            xllcorner                = xllcorner,
                            yllcorner                = yllcorner,
                            cellsizeMF               = delr[0],
                            perlen                   = perlen,
                            hnoflo                   = hnoflo
                            )
    MM_UNSAT = MMunsat.UNSAT(hnoflo = hnoflo)
    MM_SATFLOW = MMunsat.SATFLOW()

    # READ input ESRI ASCII rasters # missing gridIRR_fn
    print "\nImporting ESRI ASCII files..."
    gridMETEO = MM_PROCESS.inputEsriAscii(grid_fn = gridMETEO_fn, datatype = int)

    gridSOIL = MM_PROCESS.inputEsriAscii(grid_fn = gridSOIL_fn, datatype = int)

    gridSOILthick = MM_PROCESS.inputEsriAscii(grid_fn = gridSOILthick_fn,
     datatype = float)

    gridSUSTm = MM_PROCESS.inputEsriAscii(grid_fn = gridSUSTm_fn,
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

    _nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _Si, _Ks = MM_PROCESS.inputSoilParam(
                SOILparam_fn             = SOILparam_fn,
                NSOIL                    = NSOIL
                )
    _nslmax = max(_nsl)

    # READ observations time series (heads and soil moisture)
    obsCHECK = 1  # 0: no obs, 1: obs
    if obsCHECK==1:
        print "\nReading observations time series (hydraulic heads and soil moisture), please wait..."
        obs, outpathname, obs_h, obs_S = MM_PROCESS.inputObs(
                                inputObs_fn = inputObs_fn,
                                inputDate   = inputDate
                                )
        # Write first output in a txt file
        outFileExport = []
        for o in range(len(obs.keys())):
            outFileExport.append(open(outpathname[o], 'w'))
            S_str=''
            Rp_str=''
            Eu_str=''
            Tu_str=''
            for l in range(_nslmax):
                S_str = S_str + 'S_l' + str(l+1) + ','
                Eu_str = Eu_str + 'Eu_l' + str(l+1) + ','
                Tu_str = Tu_str + 'Tu_l' + str(l+1) + ','
                Rp_str = Rp_str + 'Rp_l' + str(l+1) + ','
                header='Date,RF,E0,PET,PE,RFe,Inter,'+Eu_str+Tu_str+'Eg,Tg,Es,'+S_str+'SUST,SUSTcount,Qs,Qscount,'+Rp_str+'R,hSATFLOW,hMF,hmeas,Smeas,SSATpart,FLOODcount,MB\n'
            outFileExport[o].write(header)
        outPESTheads_fn      = 'PESTheads.dat'
        outPESTsm_fn         = 'PESTsm.dat'
        outPESTheads=open(os.path.join(MM_ws,outPESTheads_fn), 'w')
        outPESTsm=open(os.path.join(MM_ws,outPESTsm_fn), 'w')
    else:
        print "\nNo reading of observations time series (hydraulic heads and soil moisture) required."

    # ###############
    # OUTPUT: average flux for the whole simulated period [L/T]
    outFileRF_PERall  = MM_PROCESS.outputEAgrd(outFile_fn = 'outRF.asc')
    outFilePET_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outPET.asc')
    outFilePE_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outPE.asc')
    outFileRFe_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outRFe.asc')
    outFileSUST_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outSUST.asc')
    outFileQs_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outQs.asc')
    outFileEs_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outEs.asc')
    outFileEu_PERall = []
    outFileTu_PERall = []
    outFileS_PERall = []
    outFileRp_PERall = []
    for l in range(_nslmax):
        outFileEu_PERall.append(MM_PROCESS.outputEAgrd(outFile_fn         = 'outEu_l'+str(l+1)+'.asc'))
        outFileTu_PERall.append(MM_PROCESS.outputEAgrd(outFile_fn         = 'outTu_l'+str(l+1)+'.asc'))
        outFileS_PERall.append(MM_PROCESS.outputEAgrd(outFile_fn         = 'outS_l'+str(l+1)+'.asc'))
        outFileRp_PERall.append(MM_PROCESS.outputEAgrd(outFile_fn         = 'outRp_l'+str(l+1)+'.asc'))
    outFileEg_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outEg.asc')
    outFileTg_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outTg.asc')
    outFileR_PERall = MM_PROCESS.outputEAgrd(outFile_fn  = 'outR.asc')
    outFileRn_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outRn.asc')
    outFileFLOOD_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outFLOOD.asc')
    outFileSATpart_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outSATpart.asc')
    outFileMB_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outMB.asc')
    outFileINTER_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outINTER.asc')
    outFilePOND_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outPOND.asc')
    outFileRunoff_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outRunoff.asc')

    # ###############
    # Create arrays to store output
    # arrays for fluxes independent of the soil layering
    # WARNING: R has to be the last one
    index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iSUST':4, 'iQs':5, 'iFLOOD':6, 'iEs':7, 'iSATpart':8, 'iMB':9, 'iINTER':10, 'iPOND':11, 'iRunoff':12, 'iE0':13, 'iEg':14, 'iTg':15, 'iR':16, 'iRn':17}
    # for the whole simulated period
    resavg_PERall = np.zeros([nrow,ncol,len(index)], dtype=float)
    # for each SP
    res_PERall = np.zeros([nrow,ncol,len(index),sum(perlen)], dtype=float)
    # arrays for fluxes in each soil layer
    index_S = {'iEu':0, 'iTu':1,'iS':2, 'iRp':3}
    # for the whole simulated period
    resavg_PERall_S = np.zeros([nrow,ncol,len(index_S),_nslmax], dtype=float)
    # for each SP
    res_PERall_S = np.zeros([nrow,ncol,len(index_S),_nslmax,sum(perlen)], dtype=float)
    # to compute net recharge to be exported to MF
    Rn4MF=np.zeros([nper,nrow,ncol], dtype=float)

    # ###############
    # # main loop: calculation of soil water balance in each cell-grid for each time step inside each stress period
    t0=0
    print'\n##############'
    print 'MARMITESunsat computing...'

    # initial values of SP
    Si_tmp_array = np.zeros([nrow,ncol,_nslmax], dtype=float)
    Rpi_tmp_array = np.zeros([nrow,ncol,_nslmax], dtype=float)

    for n in range(nper):
        # ###############
        # OUTPUT: average flux for each SP [L/T]
        gridoutRF = np.zeros([nrow,ncol], dtype=float)
        gridoutPET = np.zeros([nrow,ncol], dtype=float)
        gridoutPE = np.zeros([nrow,ncol], dtype=float)
        gridoutRFe = np.zeros([nrow,ncol], dtype=float)
        gridoutSUST = np.zeros([nrow,ncol], dtype=float)
        gridoutQs = np.zeros([nrow,ncol], dtype=float)
        gridoutS = []
        gridoutEu = []
        gridoutTu = []
        gridoutRp = []
        for l in range(_nslmax):
            gridoutS.append(np.zeros([nrow,ncol], dtype=float))
            gridoutEu.append(np.zeros([nrow,ncol], dtype=float))
            gridoutTu.append(np.zeros([nrow,ncol], dtype=float))
            gridoutRp.append(np.zeros([nrow,ncol], dtype=float))
        gridoutEs = np.zeros([nrow,ncol], dtype=float)
        gridoutEg = np.zeros([nrow,ncol], dtype=float)
        gridoutTg = np.zeros([nrow,ncol], dtype=float)
        gridoutRn = np.zeros([nrow,ncol], dtype=float)
        gridoutR = np.zeros([nrow,ncol], dtype=float)
        gridoutFLOOD = np.zeros([nrow,ncol], dtype=float)
        gridoutSATpart = np.zeros([nrow,ncol], dtype=float)
        gridoutMB = np.zeros([nrow,ncol], dtype=float)
        gridoutINTER = np.zeros([nrow,ncol], dtype=float)
        gridoutPOND = np.zeros([nrow,ncol], dtype=float)
        gridoutRunoff = np.zeros([nrow,ncol], dtype=float)
        outFileRPER  = MM_PROCESS.outputEAgrd(outFile_fn = rch_fn + str(n+1) + '.asc', outFolder  = MM_PROCESS.MF_ws)
        if SP_d == 0:
            outFileRF = MM_PROCESS.outputEAgrd(outFile_fn = 'outRF_PER' + str(n+1) + '.asc')
            outFilePET = MM_PROCESS.outputEAgrd(outFile_fn = 'outPET_PER' + str(n+1) + '.asc')
            outFilePE = MM_PROCESS.outputEAgrd(outFile_fn = 'outPE_PER' + str(n+1) + '.asc')
            outFileRFe = MM_PROCESS.outputEAgrd(outFile_fn = 'outRFe_PER' + str(n+1) + '.asc')
            outFileSUST = MM_PROCESS.outputEAgrd(outFile_fn = 'outSUST_PER' + str(n+1) + '.asc')
            outFileQs = MM_PROCESS.outputEAgrd(outFile_fn = 'outQs_PER' + str(n+1) + '.asc')
            outFileEu = []
            outFileTu = []
            outFileS = []
            outFileRp = []
            for l in range(_nslmax):
                outFileEu.append(MM_PROCESS.outputEAgrd(outFile_fn = 'outEu_l'+str(l+1)+'_PER' + str(n+1) + '.asc'))
                outFileTu.append(MM_PROCESS.outputEAgrd(outFile_fn = 'outTu_l'+str(l+1)+'_PER' + str(n+1) + '.asc'))
                outFileS.append(MM_PROCESS.outputEAgrd(outFile_fn = 'outS_l'+str(l+1)+'_PER' + str(n+1) + '.asc'))
                outFileRp.append(MM_PROCESS.outputEAgrd(outFile_fn = 'outRp_l'+str(l+1)+'_PER' + str(n+1) + '.asc'))
            outFileEs = MM_PROCESS.outputEAgrd(outFile_fn = 'outEs_PER' + str(n+1) + '.asc')
            outFileEg = MM_PROCESS.outputEAgrd(outFile_fn = 'outEg_PER' + str(n+1) + '.asc')
            outFileTg = MM_PROCESS.outputEAgrd(outFile_fn = 'outTg_PER' + str(n+1) + '.asc')
            outFileRn = MM_PROCESS.outputEAgrd(outFile_fn = 'outRn_PER' + str(n+1) + '.asc')
            outFileR = MM_PROCESS.outputEAgrd(outFile_fn = 'outR_PER' + str(n+1) + '.asc')
            outFileFLOOD,  = MM_PROCESS.outputEAgrd(outFile_fn = 'outFLOOD_PER' + str(n+1) + '.asc')
            outFileSATpart, = MM_PROCESS.outputEAgrd(outFile_fn = 'outSATpart_PER' + str(n+1) + '.asc')
            outFileMB, = MM_PROCESS.outputEAgrd(outFile_fn = 'outMB_PER' + str(n+1) + '.asc')
            outFileINTER, = MM_PROCESS.outputEAgrd(outFile_fn = 'outINTER_PER' + str(n+1) + '.asc')
            outFilePOND, = MM_PROCESS.outputEAgrd(outFile_fn = 'outPOND_PER' + str(n+1) + '.asc')
            outFileRunoff, = MM_PROCESS.outputEAgrd(outFile_fn = 'outRunoff_PER' + str(n+1) + '.asc')
            outFileHEADS_MF = MM_PROCESS.outputEAgrd(outFile_fn = 'outHEADS_PER' + str(n+1) + '.asc')

        tstart = 0
        for t in range(n):
            tstart = tstart + perlen[t]
        tend = tstart + perlen[n]
        results=np.zeros([nrow,ncol,len(index),perlen[n]], dtype=float)
        results_S=np.zeros([nrow,ncol,len(index_S),_nslmax,perlen[n]], dtype=float)
        ncell = 0
        for i in range(nrow):
            for j in range(ncol):
                SOILzone_tmp = gridSOIL[i,j]-1
                METEOzone_tmp = gridMETEO[i,j]-1
                if n == 0:
                    cbc[:,:,i,j,:] = conv_fact*cbc[:,:,i,j,:]/(delr[j]*delc[i])
                if ibound[i,j,0]<>0:
                    ncell = ncell + 1
                    nsl_tmp   = _nsl[SOILzone_tmp]
                    st_tmp    = _st[SOILzone_tmp]
                    slprop_tmp= _slprop[SOILzone_tmp]
                    Sm_tmp    = _Sm[SOILzone_tmp]
                    Sfc_tmp   = _Sfc[SOILzone_tmp]
                    Sr_tmp    = _Sr[SOILzone_tmp]
                    if n==0:
                        Si_tmp    = _Si[SOILzone_tmp]
                    else:
                        Si_tmp = Si_tmp_array[i,j,:]
                    Rpi_tmp = Rpi_tmp_array[i,j,:]
                    Ks_tmp    = _Ks[SOILzone_tmp]
                    SUSTm_tmp = gridSUSTm[i,j]
                    D_tmp = gridSOILthick[i,j]
                    PEsoilzonesTS_tmp = PEsoilzonesTS[METEOzone_tmp,SOILzone_tmp,tstart:tend]
                    PEsoilzonesTS_tmp = np.asarray(PEsoilzonesTS_tmp)
                    PETvegzonesTS_tmp = []
                    RFevegzonesTS_tmp = []
                    for z in range(NVEG):
                        PETvegzonesTS_tmp.append(PETvegzonesTS[METEOzone_tmp,z,tstart:tend])
                        RFevegzonesTS_tmp.append(RFevegzonesTS[METEOzone_tmp,z,tstart:tend])
                    PETvegzonesTS_tmp = np.asarray(PETvegzonesTS_tmp)
                    RFevegzonesTS_tmp = np.asarray(RFevegzonesTS_tmp)
                    VEGarea_tmp=np.zeros([NVEG], dtype=np.float)
                    for v in range(NVEG):
                        VEGarea_tmp[v]=gridVEGarea[v,i,j]
                    if n==0:
                        h_MF_tmp = np.zeros([perlen[n]], dtype = float)
                        h_MF_tmp[0] = h_MF[0,i,j,0]
                        h_MF_tmp[1:perlen[n]] = h_MF[tstart:tend-1,i,j,0]
                        cbc_tmp = np.zeros([perlen[n]], dtype = float)
                        cbc_tmp[0] = -cbc[0,iDRN,i,j,0]
                        cbc_tmp[1:perlen[n]] = -cbc[tstart:tend-1,iDRN,i,j,0]
                    else:
                        h_MF_tmp = h_MF[tstart-1:tend-1,i,j,0]
                        cbc_tmp = -cbc[tstart-1:tend-1,iDRN,i,j,0]
                    # cal functions for reservoirs calculations
                    results1_temp, results2_temp = MM_UNSAT.run(
                                                 i, j,
                                                 nsl   = nsl_tmp,
                                                 st    = st_tmp,
                                                 slprop= slprop_tmp,
                                                 Sm    = Sm_tmp,
                                                 Sfc   = Sfc_tmp,
                                                 Sr    = Sr_tmp,
                                                 Si    = Si_tmp,
                                                 Rpi   = Rpi_tmp,
                                                 D     = D_tmp,
                                                 Ks    = Ks_tmp,
                                                 SUSTm = SUSTm_tmp,
                                                 ELEV    = top[i,j],
                                                 HEADS   = h_MF_tmp,
                                                 DRN     = cbc_tmp,
                                                 RF      = RFzonesTS[METEOzone_tmp][tstart:tend],
                                                 E0      = E0zonesTS[METEOzone_tmp][tstart:tend],
                                                 PETveg  = PETvegzonesTS_tmp,
                                                 RFeveg  = RFevegzonesTS_tmp,
                                                 PEsoil  = PEsoilzonesTS_tmp,
                                                 VEGarea = VEGarea_tmp,
                                                 Zr      = Zr,
                                                 perlen  = perlen[n],
                                                 AqType  = AqType,
                                                 hdry    = hdry)
                    for k in range(len(index)-2):
                        results[i,j,k,:] = results1_temp[k,:]
                        res_PERall[i,j,k,tstart:tend] = results1_temp[k,:]
                    for k in range(len(index_S)):
                        for l in range(nsl_tmp):
                           results_S[i,j,k,l,:] = results2_temp[k,l,:]
                           res_PERall_S[i,j,k,l,tstart:tend] = results2_temp[k,l,:]
                    results[i,j,index.get('iR'),:] = results_S[i,j,index_S.get('iRp'),nsl_tmp-1,:]*1.0
                    results[i,j,index.get('iRn'),:] = results[i,j,index.get('iR'),:] - results[i,j,index.get('iEg'),:] - results[i,j,index.get('iTg'),:]
                    res_PERall[i,j,index.get('iR'),tstart:tend] = results[i,j,index.get('iR'),:]*1.0
                    res_PERall[i,j,index.get('iRn'),tstart:tend] = results[i,j,index.get('iR'),:] - results[i,j,index.get('iEg'),:] - results[i,j,index.get('iTg'),:]
                    # for the while loop
                    h_MFsum= h_MFsum + h_MF_m[:,i,j,0].sum()/nper
                     # output cumulative (1) or averaged by time (sum(perlen))
                    divfact=float(perlen[n])
                    gridoutRF[i,j]      = results[i,j,index.get('iRF'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iRF')] = resavg_PERall[i,j,index.get('iRF')] + gridoutRF[i,j]
                    gridoutPET[i,j]     = results[i,j,index.get('iPET'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iPET')] = resavg_PERall[i,j,index.get('iPET')] + gridoutPET[i,j]
                    gridoutPE[i,j]     = results[i,j,index.get('iPE'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iPE')] = resavg_PERall[i,j,index.get('iPE')] + gridoutPE[i,j]
                    gridoutRFe[i,j]     = results[i,j,index.get('iRFe'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iRFe')] = resavg_PERall[i,j,index.get('iRFe')] + gridoutRFe[i,j]
                    gridoutSUST[i,j]    = results[i,j,index.get('iSUST'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iSUST')] = resavg_PERall[i,j,index.get('iSUST')] + gridoutSUST[i,j]
                    gridoutQs[i,j]      = results[i,j,index.get('iQs'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iQs')] = resavg_PERall[i,j,index.get('iQs')] + gridoutQs[i,j]
                    gridoutEs[i,j]     = results[i,j,index.get('iEs'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iEs')] = resavg_PERall[i,j,index.get('iEs')] + gridoutEs[i,j]
                    gridoutEg[i,j]       = results[i,j,index.get('iEg'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iEg')] = resavg_PERall[i,j,index.get('iEg')] + gridoutEg[i,j]
                    gridoutTg[i,j]       = results[i,j,index.get('iTg'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iTg')] = resavg_PERall[i,j,index.get('iTg')] + gridoutTg[i,j]
                    gridoutRn[i,j]       = results[i,j,index.get('iRn'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iRn')] = resavg_PERall[i,j,index.get('iRn')] + gridoutRn[i,j]
                    gridoutR[i,j]       = results[i,j,index.get('iR'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iR')] = resavg_PERall[i,j,index.get('iR')] + gridoutR[i,j]
                    gridoutFLOOD[i,j]   = results[i,j,index.get('iFLOOD'),:].sum()
                    resavg_PERall[i,j,index.get('iFLOOD')] = resavg_PERall[i,j,index.get('iFLOOD')] + gridoutFLOOD[i,j]
                    gridoutSATpart[i,j] = results[i,j,index.get('iSATpart'),:].sum()
                    resavg_PERall[i,j,index.get('iSATpart')] = resavg_PERall[i,j,index.get('iSATpart')] + gridoutSATpart[i,j]
                    gridoutMB[i,j]      = results[i,j,index.get('iMB'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iMB')] = resavg_PERall[i,j,index.get('iMB')] + gridoutMB[i,j]
                    gridoutINTER[i,j]   = results[i,j,index.get('iINTER'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iINTER')] = resavg_PERall[i,j,index.get('iINTER')] + gridoutINTER[i,j]
                    gridoutPOND[i,j]   = results[i,j,index.get('iPOND'),:].sum()
                    resavg_PERall[i,j,index.get('iPOND')] = resavg_PERall[i,j,index.get('iPOND')] + gridoutPOND[i,j]
                    gridoutRunoff[i,j]   = results[i,j,index.get('iRunoff'),:].sum()
                    resavg_PERall[i,j,index.get('iRunoff')] = resavg_PERall[i,j,index.get('iRunoff')] + gridoutRunoff[i,j]

                    for l in range(nsl_tmp):
                        gridoutEu[l][i,j] = results_S[i,j,index_S.get('iEu'),l,:].sum()/divfact
                        resavg_PERall_S[i,j,index_S.get('iEu'),l] = resavg_PERall_S[i,j,index_S.get('iEu'),l] + gridoutEu[l][i,j]
                        gridoutTu[l][i,j] = results_S[i,j,index_S.get('iTu'),l,:].sum()/divfact
                        resavg_PERall_S[i,j,index_S.get('iTu'),l] = resavg_PERall_S[i,j,index_S.get('iTu'),l] + gridoutTu[l][i,j]
                        gridoutS[l][i,j] = results_S[i,j,index_S.get('iS'),l,:].sum()/divfact
                        resavg_PERall_S[i,j,index_S.get('iS'),l] = resavg_PERall_S[i,j,index_S.get('iS'),l] + gridoutS[l][i,j]
                        gridoutRp[l][i,j] = results_S[i,j,index_S.get('iRp'),l,:].sum()/divfact
                        resavg_PERall_S[i,j,index_S.get('iRp'),l] = resavg_PERall_S[i,j,index_S.get('iRp'),l] + gridoutRp[l][i,j]

                    Rn4MF[n,i,j] = results[i,j,index.get('iRn'),:].sum()/(divfact*1000)

                    # setting initial conditions for the next SP
                    for l in range(nsl_tmp):
                        Si_tmp_array[i,j,l] = results_S[i,j,index_S.get('iS'),l,perlen[n]-1]
                        Rpi_tmp_array[i,j,l] = results_S[i,j,index_S.get('iRp'),l,perlen[n]-1]
                    if SP_d==0:
                        outFileRF.write('%.6f'%(gridoutRF[i,j]) + ' ')
                        outFilePET.write('%.6f'%(gridoutPET[i,j]) + ' ')
                        outFilePE.write('%.6f'%(gridoutPE[i,j]) + ' ')
                        outFileRFe.write('%.6f'%(gridoutRFe[i,j]) + ' ')
                        outFileSUST.write('%.6f'%(gridoutSUST[i,j]) + ' ')
                        outFileQs.write('%.6f'%(gridoutQs[i,j]) + ' ')
                        for l in range(_nslmax):
                            outFileEu[l].write('%.6f'%(gridoutEu[l][i,j]) + ' ')
                            outFileTu[l].write('%.6f'%(gridoutTu[l][i,j]) + ' ')
                            outFileS[l].write('%.6f'%(gridoutS[l][i,j]) + ' ')
                            outFileRp[l].write('%.6f'%(gridoutRp[l][i,j]) + ' ')
                        outFileEs.write('%.6f'%(gridoutEs[i,j]) + ' ')
                        outFileEg.write('%.6f'%(gridoutEg[i,j]) + ' ')
                        outFileTg.write('%.6f'%(gridoutTg[i,j]) + ' ')
                        outFileRn.write('%.6f'%(gridoutRn[i,j]) + ' ')
                        outFileR.write('%.6f'%(gridoutR[i,j]) + ' ')
                        outFileFLOOD.write('%.6f'%(gridoutFLOOD[i,j]) + ' ')
                        outFileSATpart.write('%.6f'%(gridoutSATpart[i,j]) + ' ')
                        outFileMB.write('%.6f'%(gridoutMB[i,j]) + ' ')
                        outFileINTER.write('%.6f'%(gridoutINTER[i,j]) + ' ')
                        outFilePOND.write('%.6f'%(gridoutPOND[i,j]) + ' ')
                        outFileRunoff.write('%.6f'%(gridoutRunoff[i,j]) + ' ')
                        outFileHEADS_MF.write('%.6f'%(h_MF[tend-1,i,j,0]) + ' ')
                    if Rn4MF[n,i,j]>1e-10:
                        outFileRPER.write(str(Rn4MF[n,i,j])+' ')
                    else:
                        outFileRPER.write(str(0.0)+' ')
                else:
                    results[i,j,:,:]    = hnoflo
                    results_S[i,j,:,:,:]    = hnoflo
                    resavg_PERall[i,j,:] = hnoflo
                    resavg_PERall_S[i,j,:,:] = hnoflo
                    gridoutRF[i,j]      = hnoflo
                    gridoutPET[i,j]     = hnoflo
                    gridoutPE[i,j]      = hnoflo
                    gridoutRFe[i,j]     = hnoflo
                    gridoutSUST[i,j]    = hnoflo
                    gridoutQs[i,j]      = hnoflo
                    for l in range(nsl_tmp):
                        gridoutEu[l][i,j]     = hnoflo
                        gridoutTu[l][i,j]     = hnoflo
                        gridoutS[l][i,j]       = hnoflo
                        gridoutRp[l][i,j]      = hnoflo
                    gridoutEs[i,j]      = hnoflo
                    gridoutEg[i,j]      = hnoflo
                    gridoutTg[i,j]      = hnoflo
                    gridoutRn[i,j]      = hnoflo
                    gridoutR[i,j]       = hnoflo
                    gridoutFLOOD[i,j]   = hnoflo
                    gridoutSATpart[i,j] = hnoflo
                    gridoutMB[i,j]      = hnoflo
                    gridoutINTER[i,j]   = hnoflo
                    gridoutFLOOD[i,j]   = hnoflo
                    gridoutRunoff[i,j]  = hnoflo
                    Rn4MF[n,i,j]        = hnoflo
                    outFileRPER.write(str(hnoflo)+' ')
                    if SP_d==0:
                        outFileRF.write('%.6f'%(gridoutRF[i,j]) + ' ')
                        outFilePET.write('%.6f'%(gridoutPET[i,j]) + ' ')
                        outFilePE.write('%.6f'%(gridoutPE[i,j]) + ' ')
                        outFileRFe.write('%.6f'%(gridoutRFe[i,j]) + ' ')
                        outFileSUST.write('%.6f'%(gridoutSUST[i,j]) + ' ')
                        outFileQs.write('%.6f'%(gridoutQs[i,j]) + ' ')
                        for l in range(_nslmax):
                            outFileEu[l].write('%.6f'%(gridoutEu[l][i,j]) + ' ')
                            outFileTu[l].write('%.6f'%(gridoutTu[l][i,j]) + ' ')
                            outFileS[l].write('%.6f'%(gridoutS[l][i,j]) + ' ')
                            outFileRp[l].write('%.6f'%(gridoutRp[l][i,j]) + ' ')
                        outFileEs.write('%.6f'%(gridoutEs[i,j]) + ' ')
                        outFileEg.write('%.6f'%(gridoutEg[i,j]) + ' ')
                        outFileTg.write('%.6f'%(gridoutTg[i,j]) + ' ')
                        outFileRn.write('%.6f'%(gridoutRn[i,j]) + ' ')
                        outFileR.write('%.6f'%(gridoutR[i,j]) + ' ')
                        outFileFLOOD.write('%.6f'%(gridoutFLOOD[i,j]) + ' ')
                        outFileSATpart.write('%.6f'%(gridoutSATpart[i,j]) + ' ')
                        outFileMB.write('%.6f'%(gridoutMB[i,j]) + ' ')
                        outFileINTER.write('%.6f'%(gridoutINTER[i,j]) + ' ')
                        outFilePOND.write('%.6f'%(gridoutPOND[i,j]) + ' ')
                        outFileRunoff.write('%.6f'%(gridoutRunoff[i,j]) + ' ')
                        outFileHEADS_MF.write('%.6f'%(h_MF[tend-1,i,j,0]) + ' ')

            outFileRPER.write('\n')
            if SP_d==0:
                outFileRF.write('\n')
                outFilePET.write('\n')
                outFilePE.write('\n')
                outFileRFe.write('\n')
                outFileSUST.write('\n')
                outFileQs.write('\n')
                for l in range(_nslmax):
                    outFileEu[l].write('\n')
                    outFileTu[l].write('\n')
                    outFileS[l].write('\n')
                    outFileRp[l].write('\n')
                outFileEs.write('\n')
                outFileEg.write('\n')
                outFileTg.write('\n')
                outFileRn.write('\n')
                outFileR.write('\n')
                outFileFLOOD.write('\n')
                outFileSATpart.write('\n')
                outFileMB.write('\n')
                outFileINTER.write('\n')
                outFilePOND.write('\n')
                outFileRunoff.write('\n')
                outFileHEADS_MF.write('\n')
        # close export ASCII files
        outFileRPER.close()
        if SP_d==0:
            outFileRF.close()
            outFilePET.close()
            outFilePE.close()
            outFileRFe.close()
            outFileSUST.close()
            outFileQs.close()
            for l in range(_nslmax):
                outFileEu[l].close()
                outFileTu[l].close()
                outFileS[l].close()
                outFileRp[l].close()
            outFileEs.close()
            outFileEg.close()
            outFileTg.close()
            outFileRn.close()
            outFileR.close()
            outFileFLOOD.close()
            outFileSATpart.close()
            outFileMB.close()
            outFileINTER.close()
            outFilePOND.close()
            outFileRunoff.close()
            outFileHEADS_MF.close()

        t0=t0+int(perlen[n])
     #   print '\nSTRESS PERIOD %i/%i DONE!' % (n+1, nper)
    # #############################
    # ###  PRINT CONVERG. PLOT ###
    # #############################
    h_MFsum = h_MFsum/ncell
    h_diff.append(h_MFsum - h_pSP)
    TRY = TRY + 1
    TRYlst.append(TRY)
    h_pSP = h_MFsum
    if pylab.absolute(h_diff[TRY])>0.0:
        h_diff_log.append(pylab.log10(pylab.absolute(h_diff[TRY])))
    else:
        h_diff_log.append(pylab.log10(convcrit))
    if TRY>0:
        fig = plt.figure()
        ax1=fig.add_subplot(2,1,1)
        plt.setp(ax1.get_xticklabels(), fontsize=8)
        plt.setp(ax1.get_yticklabels(), fontsize=8)
        plt.ylabel('log(h_diff) [log(m)]', fontsize=10)
        plt.grid(True)
        plt.plot(TRYlst[1:], h_diff_log[1:], linestyle='-', marker='o', markersize=5, c = 'orange', markerfacecolor='orange', markeredgecolor='red')
        ax2=fig.add_subplot(2,1,2, sharex = ax1)
        plt.setp(ax2.get_xticklabels(), fontsize=8)
        plt.setp(ax2.get_yticklabels(), fontsize=8)
        plt.plot(TRYlst[1:], h_diff[1:], linestyle='-', marker='o', markersize=5, c = 'orange', markerfacecolor='orange', markeredgecolor='red')
        plt.ylabel('h_diff [m]', fontsize=10)
        plt.xlabel('trial', fontsize=10)
        plt.grid(True)
        if TRY>1:
            plt.xlim(1,TRY)
        ax1.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%2d'))
        ax1.xaxis.set_ticks(TRYlst[1:])
        ax2.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%2d'))
        plt.savefig(plotCONVERGENCE_export_fn)
        plt.close()
        del fig
    print "\nTRY %s\nh_diff = %.2f" % (str(TRY), h_diff[TRY])
    if TRY>ccnum:
        print'\nNo convergence between MARMITES and MODFLOW, change parameters.'
        break

plt.close()

# #############################
# ###  END CONVERGENCE LOOP ###
# #############################

# #############################
# ###  OUTPUT EXPORT   ########
# #############################

# write fluxes for the total time into grid
for i in range(nrow):
    for j in range(ncol):
        if ibound[i,j,0]<>0:
            outFileRF_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iRF')]/nper) + ' ')
            outFilePET_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iPET')]/nper) + ' ')
            outFilePE_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iPE')]/nper) + ' ')
            outFileRFe_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iRFe')]/nper) + ' ')
            outFileSUST_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iSUST')]/nper) + ' ')
            outFileQs_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iQs')]/nper) + ' ')
            for l in range(_nslmax):
                outFileEu_PERall[l].write('%.6f'%(resavg_PERall_S[i,j,index_S.get('iEu'),l]/nper) + ' ')
                outFileTu_PERall[l].write('%.6f'%(resavg_PERall_S[i,j,index_S.get('iTu'),l]/nper) + ' ')
                outFileS_PERall[l].write('%.6f'%(resavg_PERall_S[i,j,index_S.get('iS'),l]/nper) + ' ')
                outFileRp_PERall[l].write('%.6f'%(resavg_PERall_S[i,j,index_S.get('iRp'),l]/nper) + ' ')
            outFileEs_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iEs')]/nper) + ' ')
            outFileEg_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iEg')]/nper) + ' ')
            outFileTg_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iTg')]/nper) + ' ')
            outFileRn_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iRn')]/nper) + ' ')
            outFileR_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iR')]/nper) + ' ')
            outFileFLOOD_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iFLOOD')]) + ' ')
            outFileSATpart_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iSATpart')]) + ' ')
            outFileMB_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iMB')]) + ' ')
            outFileINTER_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iINTER')]/nper) + ' ')
            outFilePOND_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iPOND')]) + ' ')
            outFileRunoff_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iRunoff')]) + ' ')
        else:
            outFileRF_PERall.write('%.6f'%(hnoflo) + ' ')
            outFilePET_PERall.write('%.6f'%(hnoflo) + ' ')
            outFilePE_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileRFe_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileSUST_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileQs_PERall.write('%.6f'%(hnoflo) + ' ')
            for l in range(_nslmax):
                outFileEu_PERall[l].write('%.6f'%(hnoflo) + ' ')
                outFileTu_PERall[l].write('%.6f'%(hnoflo) + ' ')
                outFileS_PERall[l].write('%.6f'%(hnoflo) + ' ')
                outFileRp_PERall[l].write('%.6f'%(hnoflo) + ' ')
            outFileEs_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileEg_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileTg_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileRn_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileR_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileFLOOD_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileSATpart_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileMB_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileINTER_PERall.write('%.6f'%(hnoflo) + ' ')
            outFilePOND_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileRunoff_PERall.write('%.6f'%(hnoflo) + ' ')
    outFileRF_PERall.write('\n')
    outFilePET_PERall.write('\n')
    outFilePE_PERall.write('\n')
    outFileRFe_PERall.write('\n')
    outFileSUST_PERall.write('\n')
    outFileQs_PERall.write('\n')
    for l in range(_nslmax):
        outFileEu_PERall[l].write('\n')
        outFileTu_PERall[l].write('\n')
        outFileS_PERall[l].write('\n')
        outFileRp_PERall[l].write('\n')
    outFileEs_PERall.write('\n')
    outFileEg_PERall.write('\n')
    outFileTg_PERall.write('\n')
    outFileRn_PERall.write('\n')
    outFileR_PERall.write('\n')
    outFileFLOOD_PERall.write('\n')
    outFileSATpart_PERall.write('\n')
    outFileMB_PERall.write('\n')
    outFileINTER_PERall.write('\n')
    outFilePOND_PERall.write('\n')
    outFileRunoff_PERall.write('\n')
outFileRF_PERall.close()
outFilePET_PERall.close()
outFilePE_PERall.close()
outFileRFe_PERall.close()
outFileSUST_PERall.close()
outFileQs_PERall.close()
for l in range(_nslmax):
    outFileEu_PERall[l].close()
    outFileTu_PERall[l].close()
    outFileS_PERall[l].close()
    outFileRp_PERall[l].close()
outFileEs_PERall.close()
outFileEg_PERall.close()
outFileTg_PERall.close()
outFileRn_PERall.close()
outFileR_PERall.close()
outFileFLOOD_PERall.close()
outFileSATpart_PERall.close()
outFileMB_PERall.close()
outFileINTER_PERall.close()
outFilePOND_PERall.close()
outFileRunoff_PERall.close()

# final report of successful run

timeend = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
duration = duration + (timeend-timestart) -durationMF

timestartExport = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
print '\n##############\nMARMITES exporting...'

h_MF_m = np.ma.masked_values(h_MF_m, hdry, atol = 1E+25)
hmax = []
hmin = []
DRNmax = []
DRNmin = []
for L in range(nlay):
    hmax.append(np.nanmax(h_MF_m[:,:,:,:].flatten()))
    hmin.append(np.nanmin(h_MF_m[:,:,:,:].flatten()))
    DRNmax.append(np.nanmax(-cbc[:,iDRN,:,:,:]).flatten())
    DRNmin.append(np.nanmin(-cbc[:,iDRN,:,:,:]).flatten())
for o in range(len(obs.keys())):
    npa_m_tmp = np.ma.masked_values(obs_h[o], hnoflo, atol = 0.09)
    hmax.append(np.nanmax(npa_m_tmp.flatten()))
    hmin.append(np.nanmin(npa_m_tmp.flatten()))
hmax = float(np.ceil(np.nanmax(hmax)))
hmin = float(np.floor(np.nanmin(hmin)))
DRNmax = float(np.ceil(np.nanmax(DRNmax)))
DRNmin = float(np.floor(np.nanmin(DRNmin)))
# export data for observation cells and show calib graphs
if obsCHECK == 1:
    colors_nsl = CreateColors.main(hi=30, hf=50, numbcolors = (_nslmax+1))
    for o in range(len(obs.keys())):
        i = obs.get(obs.keys()[o])['i']
        j = obs.get(obs.keys()[o])['j']
        # SATFLOW
        h_satflow=MM_SATFLOW.run(res_PERall[i,j,index.get('iR'),:], float(obs.get(obs.keys()[o])['hi']),float(obs.get(obs.keys()[o])['h0']),float(obs.get(obs.keys()[o])['RC']),float(obs.get(obs.keys()[o])['STO']))
        # correct heads from MF (1 day delay)
        h_MF_tmp = np.zeros([sum(perlen)], dtype = float)
        h_MF_tmp[0] = h_MF[0,i,j,0]
        h_MF_tmp[1:sum(perlen)] = h_MF[0:sum(perlen)-1,i,j,0]
        # export ASCII file at piezometers location
        #TODO extract heads at piezo location and not center of cell
        MM_PROCESS.ExportResults(i, j, inputDate, _nslmax, res_PERall, index, res_PERall_S, index_S, h_satflow, h_MF_tmp, obs_h[o], obs_S[o], outFileExport[o], outPESTheads, outPESTsm, obs.keys()[o])
        # plot
        # DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, Sm, Sr):
        # {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iSUST':4, 'iQs':5, 'iFLOOD':6, 'iEs':7, 'iSATpart':8, 'iMB':9, 'iINTER':10, 'iPOND':11, 'iRunoff':12, 'iE0':13, 'iEg':14, 'iTg':15, 'iR':16, 'iRn':17}
        plot_export_fn = os.path.join(MM_ws, obs.keys()[o] + '.png')
        MMplot.allPLOT(
        inputDate,
        res_PERall[i,j,index.get('iRF'),:],
        res_PERall[i,j,index.get('iPET'),:],
        res_PERall[i,j,index.get('iPE'),:],
        res_PERall[i,j,index.get('iRFe'),:],
        res_PERall[i,j,index.get('iSUST'),:],
        res_PERall[i,j,index.get('iQs'),:],
        res_PERall_S[i,j,index_S.get('iEu'),0:_nsl[gridSOIL[i,j]-1],:],
        res_PERall_S[i,j,index_S.get('iTu'),0:_nsl[gridSOIL[i,j]-1],:],
        res_PERall[i,j,index.get('iEg'),:],
        res_PERall[i,j,index.get('iTg'),:],
        res_PERall_S[i,j,index_S.get('iS'),0:_nsl[gridSOIL[i,j]-1],:],
        res_PERall_S[i,j,index_S.get('iRp'),0:_nsl[gridSOIL[i,j]-1],:],
        res_PERall[i,j,index.get('iR'),:],
        res_PERall[i,j,index.get('iEs'),:],
        h_MF_tmp, h_satflow, obs_h[o], obs_S[o],
        _Sm[gridSOIL[i,j]-1],
        _Sr[gridSOIL[i,j]-1],
        hnoflo,
        plot_export_fn,
        colors_nsl,
        hmax,
        hmin
        )
        plot_exportMB_fn = os.path.join(MM_ws, obs.keys()[o] + '_MB.png')
        MMplot.plotMBerror(
        inputDate,
        res_PERall[i,j,index.get('iMB'),:],
        res_PERall[i,j,index.get('iFLOOD'),:],
        res_PERall[i,j,index.get('iSATpart'),:],
        res_PERall[i,j,index.get('iPOND'),:],
        res_PERall[i,j,index.get('iRunoff'),:],
        plot_exportMB_fn
        )
        outFileExport[o].close()
    # output for PEST
    outPESTheads.close()
    outPESTsm.close()

    # plot heads (grid + contours), DRN, etc... at specified TS
    TSlst = []
    TS = 0
    while TS < len(h_MF):
        TSlst.append(TS)
        TS = TS + 30
    TSlst.append(len(h_MF)-1)
    for TS in TSlst:
        # plot heads [m]
        V=[]
        for L in range(nlay):
            V.append(h_MF_m[TS,:,:,L])
        MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation (m)', msg = 'DRY', plttitle = 'HEADS', MM_ws = MM_ws, interval_type = 'arange', interval_diff = 0.5, Vmax = hmax, Vmin = hmin)
        # plot diff between drain elevation and heads elevation [m]
        DrnHeadsLtop = top_array_m - h_MF_m[TS,:,:,0]
        DrnHeadsLtop_m = np.ma.masked_greater(DrnHeadsLtop,0.0)
        V = [DrnHeadsLtop_m]
        diffMin = 0
        diffMax = np.nanmin(DrnHeadsLtop_m)
        MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = 1, V = V,  cmap = plt.cm.RdYlGn, CBlabel = 'diff. between DRN elev and hyd. heads elev. (m)', msg = ' - no drainage', plttitle = 'HEADSDRNdiff', MM_ws = MM_ws, interval_type = 'linspace', interval_num = 5, Vmax = diffMin, Vmin = diffMax, fmt='%.2G')
        # plot GW drainage [mm]
        V = []
        for L in range(nlay):
            V.append(-cbc[TS,iDRN,:,:,L])
        MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage (mm)', msg = '- no drainage', plttitle = 'DRN', MM_ws = MM_ws, interval_type = 'linspace', interval_num = 5, Vmin = DRNmin, Vmax = DRNmax, fmt='%.2G')

timeendExport = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
durationExport=(timeendExport-timestartExport)

print ('\n##############\nMARMITES executed successfully!')
print ('%s time steps\n%sx%s cells (rows x cols)') % (int(sum(perlen)),str(nrow),str(ncol))
print ('\nMARMITES run time: %s minute(s) and %.1f second(s)') % (str(int(duration*24.0*60.0)), (duration*24.0*60.0-int(duration*24.0*60.0))*60)
print ('MODFLOW run time: %s minute(s) and %.1f second(s)') % (str(int(durationMF*24.0*60.0)), (durationMF*24.0*60.0-int(durationMF*24.0*60.0))*60)
print ('Export run time: %s minute(s) and %.1f second(s)') % (str(int(durationExport*24.0*60.0)), (durationExport*24.0*60.0-int(durationExport*24.0*60.0))*60)
print ('\nOutput written in folder: \n%s\n##############\n') % MM_ws

##except (ValueError, TypeError, KeyboardInterrupt, IOError), e:
##    print e
##    raise e
#os.system('pause')

