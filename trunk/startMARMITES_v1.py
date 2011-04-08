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
__version__ = "0.1"
__date__ = "November 2010"

import pylab
import sys
import os
#import time
#import wx
#import struct
#import types
import numpy as np
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESsurf')
import startMARMITESsurface as startMMsurf
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESunsat_v1')
import MARMITESunsat_v1 as MARMunsat
sys.path.append(r'E:\00code\MARMITES\trunk\ppMF_FloPy')
import ppMODFLOW_flopy as ppMF
#####################################

#try:

# workspace (ws) definition
timestart = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
print '\nMARMITES starteeeeeeeeeeed!!!\n##############'

messagemanual="Please read the manual!\n(that by the way still doesn't exist...)"

# read input file (called _input.ini in the MARMITES workspace
# the first character on the first line has to be the character used to comment
# the file can contain any comments as the user wish, but the sequence of the input has to be respected
inputFile_fn = r'E:\00code\00ws\zz_TESTS\MARMITES_r13c6l2\_input.ini'
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
    MARM_ws = inputFile[l]
    l = l+1
    #run MARMITESsurface  1 is YES, 0 is NO
    MARMsurf_yn = int(inputFile[l].strip())
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
    MARMsurf_fn = inputFile[l].strip()
    l = l+1
    #run MODFLOWpreprocessing  1 is YES, 0 is NO
    MFpp_yn = int(inputFile[l].strip())
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
    SOILparam_fn = inputFile[l].strip()
    l = l+1
    IRR_fn = inputFile[l].strip()
    l = l+1
    inputObs_fn = inputFile[l].strip()
except:
    print '\nType error in the input\n%s' % (messagemanual)
    sys.exit()
fin.close()
print ('\nMARMITES workspace:\n%s\n\nMARMITESsurface workspace:\n%s\n\nMODFLOW workspace:\n%s' % (MARM_ws, inputFOLDER_fn, MF_ws))

# #############################
# ###  MARMITES SURFACES  #####
# #############################

print'\n##############'
print 'MARMITESsurface running...'

if MARMsurf_yn>0:
    MARMsurf_fn = startMMsurf.MMsurf(inputFOLDER_fn, inputFile_TS_fn, inputFile_PAR_fn, outputFILE_fn, MARM_ws)

inputFile = []
inputFile_fn = os.path.join(MARM_ws,MARMsurf_fn)
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
except:
    print 'Type error in file [' + inputFile_fn + ']%s' % (messagemanual)
    sys.exit()
fin.close()


# ##########################
# ###  MODFLOW FILES   #####
# ##########################

if MFpp_yn>0:
    print'\n##############'
    print 'MODFLOW pre-processing'
    nrow, ncol, delr, delc, perlen, nper, top, hnoflo, hdry, ibound, AqType, heads, rch_fn = ppMF.ppMF(MF_ws)

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
print 'MARMITES initialization'

# MARMITES INITIALIZATION
MARM_PROCESS = MARMunsat.process(MARM_ws                  = MARM_ws,
                        MF_ws                    = MF_ws,
                        nrow                     = nrow,
                        ncol                     = ncol,
                        xllcorner                = xllcorner,
                        yllcorner                = yllcorner,
                        cellsizeMF               = delr[0],
                        perlen                   = perlen,
                        hnoflo                   = hnoflo
                        )
MARM_UNSAT = MARMunsat.UNSAT(hnoflo = hnoflo)
MARM_LINRES = MARMunsat.LINRES()
MARM_SATFLOW = MARMunsat.SATFLOW()

# READ input ESRI ASCII rasters # missing gridIRR_fn
gridMETEO = MARM_PROCESS.inputEsriAscii(grid_fn = gridMETEO_fn, datatype = int)

gridSOIL = MARM_PROCESS.inputEsriAscii(grid_fn = gridSOIL_fn, datatype = int)

gridSOILthick = MARM_PROCESS.inputEsriAscii(grid_fn = gridSOILthick_fn,
 datatype = float)


##gridIRR = MARM_PROCESS.inputEsriAscii(grid_fn                  = gridIRR_fn)

# READ input time series and parameters   # missing IRR_fn
_Sm, _Sr, _Sfc, _Si, _Ks, _SUSTm, _n, _f, gridVEGarea, RFzonesTS, E0zonesTS, PETvegzonesTS, RFevegzonesTS, PEsoilzonesTS, inputDate = MARM_PROCESS.inputTS(
            outputFILE_fn            = outputFILE_fn,
            SOILparam_fn             = SOILparam_fn,
            NMETEO                   = NMETEO,
            NVEG                     = NVEG,
            NSOIL                    = NSOIL,
            inputDate_fn             = inputDate_fn,
            inputZON_TS_RF_fn        = inputZON_TS_RF_fn,
            inputZON_TS_PET_fn       = inputZON_TS_PET_fn,
            inputZON_TS_RFe_fn       = inputZON_TS_RFe_fn,
            inputZON_TS_PE_fn        = inputZON_TS_PE_fn,
            inputZON_TS_E0_fn        = inputZON_TS_E0_fn
 ) # IRR_fn

# READ observations time series (heads and soil moisture)
obsCHECK = 1  # 0: no obs, 1: obs
if obsCHECK==1:
    print "\nReading observations time series (hydraulic heads and soil moisture), please wait..."
    obs, outpathname, obs_h, obs_sm = MARM_PROCESS.inputObs(
                            inputObs_fn = inputObs_fn,
                            inputDate   = inputDate
                            )
    # Write first output in a txt file
    outFileExport = []
    for o in range(len(obs.keys())):
        outFileExport.append(open(outpathname[o], 'w'))
        outFileExport[o].write('Date,RF,E0,PET,RFe,Inter,ETu,Es,S,SUST,SUSTcount,Qs, Qscount,Rp,R,hSATFLOW,hMF,hmeas,Smeas,SSATpart, FLOODcount,MB\n')
    outPESTheads_fn      = 'PESTheads.dat'
    outPESTsm_fn         = 'PESTsm.dat'
    outPESTheads=open(os.path.join(MARM_ws,outPESTheads_fn), 'w')
    outPESTsm=open(os.path.join(MARM_ws,outPESTsm_fn), 'w')
else:
    print "\nNo reading of observations time series (hydraulic heads and soil moisture) required."

# ######################
#   ### main loop: calculation of soil water balance in each cell-grid for each time step inside each stress period
# ######################
outFileRF_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(outFile_fn         = 'outRF.asc')
outFilePET_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outPET.asc')
outFileRFe_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outRFe.asc')
outFileSUST_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outSUST.asc')
outFileQs_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outQs.asc')
outFileETu_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outETu.asc')
outFileEs_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outEs.asc')
outFileS_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outS.asc')
outFileRp_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outRp.asc')
outFileR_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outR.asc')
outFileFLOOD_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outFLOOD.asc')
outFileSATpart_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outSATpart.asc')
outFileMB_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outMB.asc')
outFileINTER_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outINTER.asc')
outFilePOND_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outPOND.asc')
outFileRunoff_PERall, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn         = 'outRunoff.asc')
del gridoutDUMMY

# order of output is 0-RF, 1-PET, 2-RFe, 3-SUST, 4-Qs, 5-ETu, 6-S, 7-Rp, 8-countFLOOD, 9-Es, 10-countSATpart, 11-MB, 12-R, 13-INTER, 14-countPOND, 15-countRunoff, 16-E0
results_PERall=np.zeros([nrow,ncol,17], dtype=float)
rech=np.zeros([nper,nrow,ncol], dtype=float)
t0=0
print'\n##############'
print 'MARMITES computing...'
# computing fluxes fopr each stress period in the whole grid
Si_tmp_array = np.zeros([nrow,ncol], dtype=float)
Rp_tmp_array = np.zeros([nrow,ncol], dtype=float)
for n in range(nper):
    # ######################
    # ###  create output files #####
    # ######################
    outFileRF, gridoutRF = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outRF_PER' + str(n+1) + '.asc')
    outFilePET, gridoutPET = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outPET_PER' + str(n+1) + '.asc')
    outFileRFe, gridoutRFe = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outRFe_PER' + str(n+1) + '.asc')
    outFileSUST, gridoutSUST = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outSUST_PER' + str(n+1) + '.asc')
    outFileQs, gridoutQs = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outQs_PER' + str(n+1) + '.asc')
    outFileETu, gridoutETu = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outETu_PER' + str(n+1) + '.asc')
    outFileEs, gridoutEs = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outEs_PER' + str(n+1) + '.asc')
    outFileS, gridoutS = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outS_PER' + str(n+1) + '.asc')
    outFileRp, gridoutRp = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outRp_PER' + str(n+1) + '.asc')
    outFileR, gridoutR = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outR_PER' + str(n+1) + '.asc')
    outFileRPER, gridoutDUMMY = MARM_PROCESS.outputEAgrd(
                        outFile_fn               = rch_fn[n],
                        outFolder                = MARM_PROCESS.MF_ws)
    del gridoutDUMMY
    outFileFLOOD, gridoutFLOOD = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outFLOOD_PER' + str(n+1) + '.asc')
    outFileSATpart, gridoutSATpart = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outSATpart_PER' + str(n+1) + '.asc')
    outFileMB, gridoutMB = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outMB_PER' + str(n+1) + '.asc')
    outFileINTER, gridoutINTER = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outINTER_PER' + str(n+1) + '.asc')
    outFilePOND, gridoutPOND = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outPOND_PER' + str(n+1) + '.asc')
    outFileRunoff, gridoutRunoff = MARM_PROCESS.outputEAgrd(
                            outFile_fn         = 'outRunoff_PER' + str(n+1) + '.asc')
    tstart = 0
    for t in range(n):
        tstart = tstart + perlen[t]
    tend = tstart + perlen[n]
    # order of output is 0-RF, 1-PET, 2-RFe, 3-SUST, 4-Qs, 5-ETu, 6-S, 7-Rp, 8-countFLOOD, 9-Es, 10-countSATpart, 11-MB, 12-R, 13-INTER, 14-countPOND, 15-countRunoff, 16-E0
    results=np.zeros([nrow,ncol,17,perlen[n]], dtype=float)
    for i in range(nrow):
        for j in range(ncol):
            SOILzone_tmp = gridSOIL[i,j]-1
            METEOzone_tmp = gridMETEO[i,j]-1
            if ibound[i,j,0]<>0:
                Sm_tmp    = _Sm[SOILzone_tmp]
                Sfc_tmp   = _Sfc[SOILzone_tmp]
                Sr_tmp    = _Sr[SOILzone_tmp]
                if n==0:
                    Si_tmp    = _Si[SOILzone_tmp]
                else:
                    Si_tmp = Si_tmp_array[i,j]
                n_tmp=_n[SOILzone_tmp]
                f_tmp=_f[SOILzone_tmp]
                Ks_tmp    = _Ks[SOILzone_tmp]
                SUSTm_tmp = _SUSTm[SOILzone_tmp]
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
                # cal functions for reservoirs calculations
                # order of output is 0-RF, 1-PET, 2-RFe, 3-SUST, 4-Qs, 5-ETu, 6-S, 7-Rp, 8-countFLOOD, 9-Es, 10-countSATpart, 11-MB, 12-R
                results_temp = MARM_UNSAT.run(i, j,
                                             Sm    = Sm_tmp,
                                             Sfc   = Sfc_tmp,
                                             Sr    = Sr_tmp,
                                             Si    = Si_tmp,
                                             D     = D_tmp,
                                             Ks    = Ks_tmp,
                                             SUSTm = SUSTm_tmp,
                                             ELEV    = top[i,j],
                                             HEADS   = heads[tstart:tend,i,j,0],
                                             RF      = RFzonesTS[METEOzone_tmp][tstart:tend],
                                             E0      = E0zonesTS[METEOzone_tmp][tstart:tend],
                                             PETveg  = PETvegzonesTS_tmp,
                                             RFeveg  = RFevegzonesTS_tmp,
                                             PEsoil  = PEsoilzonesTS_tmp,
                                             VEGarea = VEGarea_tmp,
                                             perlen  = perlen[n],
                                             AqType  = AqType,
                                             hdry    = hdry)
                #RF, PET_tot, RFe_tot, SUST, Qs, ETu, S, Rp, countFLOOD, Es, countSATpart, MB, INTER, countPOND, countRunoff, E0
                iRF=0
                iPET=1
                iRFe=2
                iSUST=3
                iQs=4
                iETu=5
                iS=6
                iRp=7
                iFLOOD=8
                iEs=9
                iSATpart=10
                iMB=11
                iINTER=12
                iPOND=13
                iRunoff=14
                iE0=15
                # WARNING: R has to be the last one
                iR=16
                for k in range(iR):
                    results[i,j,k,:]=results_temp[k]
                Rp_tmp = results[i,j,iRp,:]*1.0
                Rp_tmp[0] = Rp_tmp[0] + Rp_tmp_array[i,j]
                results[i,j,iR,:] = MARM_LINRES.run(Rp_tmp, n_tmp, f_tmp)

                 # output cumulative (1) or averaged by time (sum(perlen))
                divfact=float(perlen[n])
                gridoutRF[i,j]      = results[i,j,iRF,:].sum()/divfact
                results_PERall[i,j,iRF] = results_PERall[i,j,iRF] + gridoutRF[i,j]
                gridoutPET[i,j]     = results[i,j,iPET,:].sum()/divfact
                results_PERall[i,j,iPET] = results_PERall[i,j,iPET] + gridoutPET[i,j]
                gridoutRFe[i,j]     = results[i,j,iRFe,:].sum()/divfact
                results_PERall[i,j,iRFe] = results_PERall[i,j,iRFe] + gridoutRFe[i,j]
                gridoutSUST[i,j]    = results[i,j,iSUST,:].sum()/divfact
                results_PERall[i,j,iSUST] = results_PERall[i,j,iSUST] + gridoutSUST[i,j]
                gridoutQs[i,j]      = results[i,j,iQs,:].sum()/divfact
                results_PERall[i,j,iQs] = results_PERall[i,j,iQs] + gridoutQs[i,j]
                gridoutETu[i,j]     = results[i,j,iETu,:].sum()/divfact
                results_PERall[i,j,iETu] = results_PERall[i,j,iETu] + gridoutETu[i,j]
                gridoutEs[i,j]     = results[i,j,iEs,:].sum()/divfact
                results_PERall[i,j,iEs] = results_PERall[i,j,iEs] + gridoutEs[i,j]
                gridoutS[i,j]       = results[i,j,iS,:].sum()/divfact
                results_PERall[i,j,iS] = results_PERall[i,j,iS] + gridoutS[i,j]
                gridoutRp[i,j]      = results[i,j,iRp,:].sum()/divfact
                results_PERall[i,j,iRp] = results_PERall[i,j,iRp] + gridoutRp[i,j]
                gridoutR[i,j]       = results[i,j,iR,:].sum()/divfact
                results_PERall[i,j,iR] = results_PERall[i,j,iR] + gridoutR[i,j]
                gridoutFLOOD[i,j]   = results[i,j,iFLOOD,:].sum()
                results_PERall[i,j,iFLOOD] = results_PERall[i,j,iFLOOD] + gridoutFLOOD[i,j]
                gridoutSATpart[i,j] = results[i,j,iSATpart,:].sum()
                results_PERall[i,j,iSATpart] = results_PERall[i,j,iSATpart] + gridoutSATpart[i,j]
                gridoutMB[i,j]      = results[i,j,iMB,:].sum()/divfact
                results_PERall[i,j,iMB] = results_PERall[i,j,iMB] + gridoutMB[i,j]
                gridoutINTER[i,j]   = results[i,j,iINTER,:].sum()/divfact
                results_PERall[i,j,iINTER] = results_PERall[i,j,iINTER] + gridoutINTER[i,j]
                gridoutPOND[i,j]   = results[i,j,iPOND,:].sum()/divfact
                results_PERall[i,j,iPOND] = results_PERall[i,j,iPOND] + gridoutPOND[i,j]
                gridoutRunoff[i,j]   = results[i,j,iRunoff,:].sum()/divfact
                results_PERall[i,j,iRunoff] = results_PERall[i,j,iRunoff] + gridoutRunoff[i,j]
                rech[n,i,j] = results[i,j,iR,:].sum()/(divfact*1000)

                # setting initial conditions for the next SP
                Si_tmp_array[i,j] = results[i,j,iS,perlen[n]-1]
                Rp_tmp_array[i,j] = Rp_tmp.sum() - results[i,j,iR,:].sum()

                outFileRF.write('%.6f'%(gridoutRF[i,j]) + ' ')
                outFilePET.write('%.6f'%(gridoutPET[i,j]) + ' ')
                outFileRFe.write('%.6f'%(gridoutRFe[i,j]) + ' ')
                outFileSUST.write('%.6f'%(gridoutSUST[i,j]) + ' ')
                outFileQs.write('%.6f'%(gridoutQs[i,j]) + ' ')
                outFileETu.write('%.6f'%(gridoutETu[i,j]) + ' ')
                outFileEs.write('%.6f'%(gridoutEs[i,j]) + ' ')
                outFileS.write('%.6f'%(gridoutS[i,j]) + ' ')
                outFileRp.write('%.6f'%(gridoutRp[i,j]) + ' ')
                outFileR.write('%.6f'%(gridoutR[i,j]) + ' ')
                outFileFLOOD.write('%.6f'%(gridoutFLOOD[i,j]) + ' ')
                outFileSATpart.write('%.6f'%(gridoutSATpart[i,j]) + ' ')
                outFileMB.write('%.6f'%(gridoutMB[i,j]) + ' ')
                outFileINTER.write('%.6f'%(gridoutINTER[i,j]) + ' ')
                outFilePOND.write('%.6f'%(gridoutPOND[i,j]) + ' ')
                outFileRunoff.write('%.6f'%(gridoutRunoff[i,j]) + ' ')
                if rech[n,i,j]>1e-10:
                    outFileRPER.write(str(rech[n,i,j])+' ')
                else:
                    outFileRPER.write(str(0.0)+' ')
            else:
                results[i,j,:,:]    = hnoflo
                results_PERall[i,j,:] = hnoflo
                gridoutRF[i,j]      = hnoflo
                gridoutPET[i,j]     = hnoflo
                gridoutRFe[i,j]     = hnoflo
                gridoutSUST[i,j]    = hnoflo
                gridoutQs[i,j]      = hnoflo
                gridoutETu[i,j]     = hnoflo
                gridoutEs[i,j]     = hnoflo
                gridoutS[i,j]       = hnoflo
                gridoutRp[i,j]      = hnoflo
                gridoutR[i,j]       = hnoflo
                gridoutFLOOD[i,j]   = hnoflo
                gridoutSATpart[i,j] = hnoflo
                gridoutMB[i,j]      = hnoflo
                gridoutINTER[i,j]   = hnoflo
                gridoutFLOOD[i,j]   = hnoflo
                gridoutRunoff[i,j]   = hnoflo
                rech[n,i,j]         = hnoflo
                outFileRF.write('%.6f'%(gridoutRF[i,j]) + ' ')
                outFilePET.write('%.6f'%(gridoutPET[i,j]) + ' ')
                outFileRFe.write('%.6f'%(gridoutRFe[i,j]) + ' ')
                outFileSUST.write('%.6f'%(gridoutSUST[i,j]) + ' ')
                outFileQs.write('%.6f'%(gridoutQs[i,j]) + ' ')
                outFileETu.write('%.6f'%(gridoutETu[i,j]) + ' ')
                outFileEs.write('%.6f'%(gridoutEs[i,j]) + ' ')
                outFileS.write('%.6f'%(gridoutS[i,j]) + ' ')
                outFileRp.write('%.6f'%(gridoutRp[i,j]) + ' ')
                outFileR.write('%.6f'%(gridoutR[i,j]) + ' ')
                outFileFLOOD.write('%.6f'%(gridoutFLOOD[i,j]) + ' ')
                outFileSATpart.write('%.6f'%(gridoutSATpart[i,j]) + ' ')
                outFileMB.write('%.6f'%(gridoutMB[i,j]) + ' ')
                outFileINTER.write('%.6f'%(gridoutINTER[i,j]) + ' ')
                outFilePOND.write('%.6f'%(gridoutPOND[i,j]) + ' ')
                outFileRunoff.write('%.6f'%(gridoutRunoff[i,j]) + ' ')
                outFileRPER.write(str(hnoflo)+' ')
            # export data for observation cells and show calib graphs
            if obsCHECK==1:
                for o in range(len(obs.keys())):
                    if obs.get(obs.keys()[o])['i'] == i and obs.get(obs.keys()[o])['j'] == j:
                     # order of output is 0-RF, 1-PET, 2-RFe, 3-SUST, 4-Qs, 5-ETu, 6-S, 7-Rp, 8-countFLOOD, 9-Es, 10-countSATpart, 11-MB, 12-R, 13-INTER, 14-countPOND, 15-countRunoff
                        h_satflow=MARM_SATFLOW.run(results[i,j,iR,:], float(obs.get(obs.keys()[o])['hi']),float(obs.get(obs.keys()[o])['h0']),float(obs.get(obs.keys()[o])['RC']),float(obs.get(obs.keys()[o])['STO']))
                        # reinitialize hi for the next time period
                        obs.get(obs.keys()[o])['hi']=h_satflow[perlen[n]-1]-float(obs.get(obs.keys()[o])['h0'])
                        TS = [results[i,j,iRF,:], results[i,j,iPET,:], results[i,j,iRFe,:],results[i,j,iINTER,:], results[i,j,iETu,:], results[i,j,iEs,:],results[i,j,iS,:], results[i,j,iSUST,:],results[i,j,iPOND,:],results[i,j,iQs,:],results[i,j,iRunoff,:],results[i,j,iFLOOD,:],results[i,j,iRp,:], results[i,j,iR,:], h_satflow, heads[tstart:tend,i,j,0], results[i,j,iSATpart,:],results[i,j,iMB,:], results[i,j,iE0,:]]
                        MARM_PROCESS.ExportResults(inputDate[tstart:tend], TS, obs_h[o][tstart:tend], obs_sm[o][tstart:tend], outFileExport[o], outPESTheads, outPESTsm, obs.keys()[o])
                ##    DLPgraphs.calibGRAPH(inputDate, results[i,j,0,:], results[i,j,1,:],
                ##                          results[i,j,2,:], results[i,j,5,:], results[i,j,6,:],
                ##                          results[i,j,9,:], h, obsh[o], obssm[o], Sm[gridSOIL[i,j]-1], Sr[gridSOIL[i,j]-1])
        outFileRPER.write('\n')
        outFileRF.write('\n')
        outFilePET.write('\n')
        outFileRFe.write('\n')
        outFileSUST.write('\n')
        outFileQs.write('\n')
        outFileETu.write('\n')
        outFileEs.write('\n')
        outFileS.write('\n')
        outFileRp.write('\n')
        outFileR.write('\n')
        outFileFLOOD.write('\n')
        outFileSATpart.write('\n')
        outFileMB.write('\n')
        outFileINTER.write('\n')
        outFilePOND.write('\n')
        outFileRunoff.write('\n')
#           print ('Row '+str(i+1)+'/'+str(nrow)+' done')
    # close export ASCII files
    outFileRF.close()
    outFilePET.close()
    outFileRFe.close()
    outFileSUST.close()
    outFileQs.close()
    outFileETu.close()
    outFileEs.close()
    outFileS.close()
    outFileRp.close()
    outFileR.close()
    outFileFLOOD.close()
    outFileSATpart.close()
    outFileMB.close()
    outFileINTER.close()
    outFilePOND.close()
    outFileRunoff.close()
    outFileRPER.close()
    t0=t0+int(perlen[n])
    print '\nSTRESS PERIOD %i/%i DONE!' % (n+1, nper)
# write fluxes for the total time into grid
for i in range(nrow):
    for j in range(ncol):
        if ibound[i,j,0]<>0:
            outFileRF_PERall.write('%.6f'%(results_PERall[i,j,iRF]/nper) + ' ')
            outFilePET_PERall.write('%.6f'%(results_PERall[i,j,iPET]/nper) + ' ')
            outFileRFe_PERall.write('%.6f'%(results_PERall[i,j,iRFe]/nper) + ' ')
            outFileSUST_PERall.write('%.6f'%(results_PERall[i,j,iSUST]/nper) + ' ')
            outFileQs_PERall.write('%.6f'%(results_PERall[i,j,iQs]/nper) + ' ')
            outFileETu_PERall.write('%.6f'%(results_PERall[i,j,iETu]/nper) + ' ')
            outFileEs_PERall.write('%.6f'%(results_PERall[i,j,iEs]/nper) + ' ')
            outFileS_PERall.write('%.6f'%(results_PERall[i,j,iS]/nper) + ' ')
            outFileRp_PERall.write('%.6f'%(results_PERall[i,j,iRp]/nper) + ' ')
            outFileR_PERall.write('%.6f'%(results_PERall[i,j,iR]/nper) + ' ')
            outFileFLOOD_PERall.write('%.6f'%(results_PERall[i,j,iFLOOD]) + ' ')
            outFileSATpart_PERall.write('%.6f'%(results_PERall[i,j,iSATpart]) + ' ')
            outFileMB_PERall.write('%.6f'%(results_PERall[i,j,iMB]) + ' ')
            outFileINTER_PERall.write('%.6f'%(results_PERall[i,j,iINTER]/nper) + ' ')
            outFilePOND_PERall.write('%.6f'%(results_PERall[i,j,iPOND]) + ' ')
            outFileRunoff_PERall.write('%.6f'%(results_PERall[i,j,iRunoff]) + ' ')
        else:
            outFileRF_PERall.write('%.6f'%(hnoflo) + ' ')
            outFilePET_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileRFe_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileSUST_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileQs_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileETu_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileEs_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileS_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileRp_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileR_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileFLOOD_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileSATpart_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileMB_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileINTER_PERall.write('%.6f'%(hnoflo) + ' ')
            outFilePOND_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileRunoff_PERall.write('%.6f'%(hnoflo) + ' ')
    outFileRF_PERall.write('\n')
    outFilePET_PERall.write('\n')
    outFileRFe_PERall.write('\n')
    outFileSUST_PERall.write('\n')
    outFileQs_PERall.write('\n')
    outFileETu_PERall.write('\n')
    outFileEs_PERall.write('\n')
    outFileS_PERall.write('\n')
    outFileRp_PERall.write('\n')
    outFileR_PERall.write('\n')
    outFileFLOOD_PERall.write('\n')
    outFileSATpart_PERall.write('\n')
    outFileMB_PERall.write('\n')
    outFileINTER_PERall.write('\n')
    outFilePOND_PERall.write('\n')
    outFileRunoff_PERall.write('\n')
outFileRF_PERall.close()
outFilePET_PERall.close()
outFileRFe_PERall.close()
outFileSUST_PERall.close()
outFileQs_PERall.close()
outFileETu_PERall.close()
outFileEs_PERall.close()
outFileS_PERall.close()
outFileRp_PERall.close()
outFileR_PERall.close()
outFileFLOOD_PERall.close()
outFileSATpart_PERall.close()
outFileMB_PERall.close()
outFileINTER_PERall.close()
outFilePOND_PERall.close()
outFileRunoff_PERall.close()

if obsCHECK == 1:
    for o in range(len(obs.keys())):
        outFileExport[o].close()
    outPESTheads.close()
    outPESTsm.close()

# final report of successful run
timeend = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
duration=(timeend-timestart)
print ('\n############### REPORT ####################\nMARMITES executed successfully!')
print ('%s time steps\n%sx%s cells (rows x cols)') % (int(sum(perlen)),str(nrow),str(ncol))
print ('Run time: %s minute(s) and %.1f second(s)') % (str(int(duration*24.0*60.0)), (duration*24.0*60.0-int(duration*24.0*60.0))*60)
print ('Output files written in folder: \n%s\n################# EOF #####################\n') % MARM_ws
# ReSoGraph.calibGRAPH(DateInput, P, PET, Pe, ETu, S, R, h, hmeas, Smeas, Sm, Sr)

del gridoutRF, gridoutRFe, gridoutR, gridoutRp, gridoutETu, gridoutEs, gridoutS, results_temp, results

##except (ValueError, TypeError, KeyboardInterrupt, IOError), e:
##    print e
##    raise e
#    os.system('pause')

