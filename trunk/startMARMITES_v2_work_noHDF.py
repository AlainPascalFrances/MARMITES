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
import ppMODFLOW_flopy_noHDF as ppMF
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESplot')
import MARMITESplot as MMplot
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESsurf')
import CreateColors
import StringIO

#####################################
#try:
# workspace (ws) definition
timestart = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
print '\n##############\nMARMITES started!\n##############'

messagemanual="Please read the manual!\n(that by the way still doesn't exist...)"

# read input file (called _input.ini in the MARMITES workspace
# the first character on the first line has to be the character used to comment
# the file can contain any comments as the user wish, but the sequence of the input has to be respected
inputFile = []
inputFile_fn = r'E:\00code\00ws\zz_TESTS\MARMITESv2_r13c6l2\_inputMM.ini'
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
    print 'Error in the input file %s, check format!\n%s' % (inputFile_fn, messagemanual)
    sys.exit()
except e:
    print "Unexpected error in the input file %s:\n", sys.exc_info()[0] % (inputFile_fn)
    print messagemanual
    sys.exit()
l=0
try:
    # report file (0 - report and no interpreter verbose, 1 - no report and interpreter verbose)
    verbose = int(inputFile[l].strip())
    l = l + 1
    # output plot (1 is YES, 0 is NO)
    plot_out  = int(inputFile[l].strip())
    l = l + 1
    # Define MARMITES ws folders
    MM_ws = inputFile[l]
    l = l+1
    #run MARMITESsurface  (1 is YES, 0 is NO)
    MMsurf_yn = int(inputFile[l].strip())
    l = l+1
    # Define MARMITESsurface folder
    MMsurf_ws = inputFile[l].strip()
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
    MF_ws = inputFile[l].strip()
    l = l+1
    rch_input = inputFile[l].strip()
    l = l+1
    wel_input = inputFile[l].strip()
    l = l+1
    plot_freq =  int(inputFile[l].strip())
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
    gridPONDhmax_fn =  inputFile[l].strip()
    l = l+1
    gridPONDw_fn =  inputFile[l].strip()
    l = l+1
    SOILparam_fn = inputFile[l].strip()
    l = l+1
    IRR_fn = inputFile[l].strip()
    l = l+1
    inputObs_fn = inputFile[l].strip()
    l = l+1
    inputObsHEADS_fn = inputFile[l].strip()
    l = l+1
    inputObsSM_fn = inputFile[l].strip()
    l = l+1
    convcrit = float(inputFile[l].strip())
    l = l+1
    ccnum = int(inputFile[l].strip())
    l = l+1
    plt_ConvLoop_fn = inputFile[l].strip()
except:
    print '\nType error in the input file %s\n%s' % (inputFile_fn, messagemanual)
    sys.exit()
fin.close()

if verbose == 0:
#capture interpreter output to be written in to a report file
    s = StringIO.StringIO()
    sys.stdout = s
    report_fn = os.path.join(MM_ws,'00_MM_report.txt')
    report = open(report_fn, 'w')
    print '\n##############\nMARMITES started!\n##############'
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
# ###  MARMITES SURFACES  #####
# #############################

print'\n##############'
print 'MARMITESsurf RUN'

if MMsurf_yn>0:
    MMsurf_fn = startMMsurf.MMsurf(MMsurf_ws, inputFile_TS_fn, inputFile_PAR_fn, outputFILE_fn, MM_ws)

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
k_Tu_slp = []
k_Tu_inter = []
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
        k_Tu_slp.append(float(line[v]))
    l = l +1
    line = inputFile[l].split()
    for v in range(NVEG):
        k_Tu_inter.append(float(line[v]))
    l = l +1
    line = inputFile[l].split()
    for v in range(NSOIL):
        TRANS_sdw.append(int(line[v]))
except:
    print 'Type error in file [' + inputFile_fn + ']%s' % (messagemanual)
    sys.exit()
fin.close()

# #############################
# ### 1st MODFLOW RUN with initial user-input recharge
# #############################
durationMF = 0.0
timestartMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
print'\n##############'
print 'MODFLOW RUN (initial user-input fluxes)'
SP_d, nrow, ncol, delr, delc, nlay, perlen, nper, top, hnoflo, hdry, ibound, laytyp, h_MF, cbc, cbc_nam_tmp, top_array, inputFileMF_fn, lenuni = ppMF.ppMF(MF_ws, rch_input = rch_input, rch_dft = 0.001, wel_input = wel_input)

cbc_nam = []
for c in cbc_nam_tmp:
    cbc_nam.append(c.strip())
iDRN = cbc_nam.index('DRAINS')
iSTO = cbc_nam.index('STORAGE')
iRCH = cbc_nam.index('RECHARGE')
iWEL = cbc_nam.index('WELLS')
# convert cbc from volume to mm
# cbc format is: (kstp), kper, textprocess, ncol, nrow, nlay
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

timeendMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
durationMF = durationMF + (timeendMF-timestartMF)

# ####   SUMMARY OF MODFLOW READINGS   ####
# active cells in layer 1_____________________ibound[0]
# elevation___________________________________top
# heads in layer n____________________________heads[timestep, row, col, layer]
# aquifer type of layer 1_____________________laytyp[0]
# numero total de time step___________________sum(perlen)
# code for dry cell___________________________hdry
# cell size___________________________________delr[0]

print'\n##############'
print 'MARMITESunsat initialization'

# MARMITES INITIALIZATION
MM_PROCESS = MMunsat.PROCESS(MM_ws               = MM_ws,
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
                                inputZON_TS_E0_fn        = inputZON_TS_E0_fn
 ) # IRR_fn

# SOIL PARAMETERS
_nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _Si, _Ks = MM_PROCESS.inputSoilParam(SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
_nslmax = max(_nsl)

# READ observations time series (heads and soil moisture)
obsCHECK = 1  # 0: no obs, 1: obs
if obsCHECK==1:
    print "\nReading observations time series (hydraulic heads and soil moisture)..."
    obs, outpathname, obs_h, obs_S = MM_PROCESS.inputObs(
                            inputObs_fn = inputObs_fn,
                            inputObsHEADS_fn = inputObsHEADS_fn,
                            inputObsSM_fn = inputObsSM_fn,
                            inputDate   = inputDate,
                            _nslmax    = _nslmax
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
            Smeasout = Smeasout+ 'Smeas_' + str(l+1) + ','
        header='Date,RF,E0,PET,PE,RFe,Inter,'+Eu_str+Tu_str+'Eg,Tg,ETg,WEL_MF,Es,'+S_str+Spc_str+dS_str+'dPOND,POND,Ro,SEEPAGE,DRN_MF,'+Rp_str+'R,Rn,R_MF,hSATFLOW,hMF,hmeas,' + Smeasout + 'MB\n'
        outFileExport[o].write(header)
    outPESTheads_fn      = 'PESTheads.dat'
    outPESTsm_fn         = 'PESTsm.dat'
    outPESTheads=open(os.path.join(MM_ws,outPESTheads_fn), 'w')
    outPESTsm=open(os.path.join(MM_ws,outPESTsm_fn), 'w')
else:
    print "\nNo reading of observations time series (hydraulic heads and soil moisture) required."

h_pSP = 0
LOOP = 0
LOOPlst = [LOOP]
h_diff = [10]
h_diff_log = [1]

plt_ConvLoop_fn = os.path.join(MM_ws, plt_ConvLoop_fn + '.png')
duration = 0.0

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
    # OUTPUT: average flux for the whole simulated period [L/T]
    outFileRF_PERall  = MM_PROCESS.outputEAgrd(outFile_fn = 'outRF.asc')
    outFilePET_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outPET.asc')
    outFilePE_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outPE.asc')
    outFileRFe_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outRFe.asc')
    outFiledPOND_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outdPOND.asc')
    outFilePOND_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outPOND.asc')
    outFileRo_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outRo.asc')
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
    outFileETg_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outETg.asc')
    outFileSEEPAGE_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outSEEPAGE.asc')
    outFileMB_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outMB.asc')
    outFileINTER_PERall = MM_PROCESS.outputEAgrd(outFile_fn = 'outINTER.asc')

    # ###############
    # Create arrays to store output
    # arrays for fluxes independent of the soil layering
    index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iPOND':4, 'iRo':5, 'iSEEPAGE':6, 'iEs':7, 'iMB':8, 'iINTER':9, 'iE0':10, 'iEg':11, 'iTg':12, 'iR':13, 'iRn':14, 'idPOND':15, 'iETg':16}
    # for the whole simulated period
    resavg_PERall = np.zeros([nrow,ncol,len(index)], dtype=float)
    # for each SP
    res_PERall = np.zeros([nrow,ncol,len(index),sum(perlen)], dtype=float)
    # arrays for fluxes in each soil layer
    index_S = {'iEu':0, 'iTu':1,'iSpc':2, 'iRp':3, 'idS':4, 'iS':5}
    # for the whole simulated period
    resavg_PERall_S = np.zeros([nrow,ncol,len(index_S),_nslmax], dtype=float)
    # for each SP
    res_PERall_S = np.zeros([nrow,ncol,len(index_S),_nslmax,sum(perlen)], dtype=float)
    # to compute net recharge to be exported to MF
    R4MF=np.zeros([nper,nrow,ncol], dtype=float)
    ETg4MF=np.zeros([nper,nrow,ncol], dtype=float)

    # ###############
    # # main loop: calculation of soil water balance in each cell-grid for each time step inside each stress period
    t0=0
    print '\nComputing...'

    # initial values of SP
    Si_tmp_array = np.zeros([nrow,ncol,_nslmax], dtype=float)
    Rpi_tmp_array = np.zeros([nrow,ncol,_nslmax], dtype=float)
    PONDi_tmp_array = np.zeros([nrow,ncol], dtype=float)
    DRNi_tmp_array = np.zeros([nrow,ncol], dtype=float)

    for n in range(nper):
        # ###############
        # OUTPUT
        gridoutRF = np.zeros([nrow,ncol], dtype=float)
        gridoutPET = np.zeros([nrow,ncol], dtype=float)
        gridoutPE = np.zeros([nrow,ncol], dtype=float)
        gridoutRFe = np.zeros([nrow,ncol], dtype=float)
        gridoutdPOND = np.zeros([nrow,ncol], dtype=float)
        gridoutPOND = np.zeros([nrow,ncol], dtype=float)
        gridoutRo = np.zeros([nrow,ncol], dtype=float)
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
        gridoutETg = np.zeros([nrow,ncol], dtype=float)
        gridoutR = np.zeros([nrow,ncol], dtype=float)
        gridoutSEEPAGE = np.zeros([nrow,ncol], dtype=float)
        gridoutMB = np.zeros([nrow,ncol], dtype=float)
        gridoutINTER = np.zeros([nrow,ncol], dtype=float)
        outFileRPER  = MM_PROCESS.outputEAgrd(outFile_fn = rch_input + str(n+1) + '.asc', outFolder  = MM_PROCESS.MF_ws)
        outFileETgPER  = MM_PROCESS.outputEAgrd(outFile_fn = wel_input + str(n+1) + '.asc', outFolder  = MM_PROCESS.MF_ws)
        if SP_d == 0:
            # average flux for each SP [L/T]
            outFileRF = MM_PROCESS.outputEAgrd(outFile_fn = 'outRF_PER' + str(n+1) + '.asc')
            outFilePET = MM_PROCESS.outputEAgrd(outFile_fn = 'outPET_PER' + str(n+1) + '.asc')
            outFilePE = MM_PROCESS.outputEAgrd(outFile_fn = 'outPE_PER' + str(n+1) + '.asc')
            outFileRFe = MM_PROCESS.outputEAgrd(outFile_fn = 'outRFe_PER' + str(n+1) + '.asc')
            outFiledPOND = MM_PROCESS.outputEAgrd(outFile_fn = 'outdPOND_PER' + str(n+1) + '.asc')
            outFilePOND = MM_PROCESS.outputEAgrd(outFile_fn = 'outPOND_PER' + str(n+1) + '.asc')
            outFileRo = MM_PROCESS.outputEAgrd(outFile_fn = 'outRo_PER' + str(n+1) + '.asc')
            outFileEu = []
            outFileTu = []
            outFileS = []
            outFileRp = []
            for l in range(_nslmax):
                outFileEu.append(MM_PROCESS.outputEAgrd(outFile_fn = 'outEu_l'+str(l+1)+'_PER' + str(n+1) + '.asc'))
                outFileTu.append(MM_PROCESS.outputEAgrd(outFile_fn = 'outTu_l'+str(l+1)+'_PER' + str(n+1) + '.asc'))
                outFileS.append(MM_PROCESS.outputEAgrd(outFile_fn = 'outSpc_l'+str(l+1)+'_PER' + str(n+1) + '.asc'))
                outFileRp.append(MM_PROCESS.outputEAgrd(outFile_fn = 'outRp_l'+str(l+1)+'_PER' + str(n+1) + '.asc'))
            outFileEs = MM_PROCESS.outputEAgrd(outFile_fn = 'outEs_PER' + str(n+1) + '.asc')
            outFileEg = MM_PROCESS.outputEAgrd(outFile_fn = 'outEg_PER' + str(n+1) + '.asc')
            outFileTg = MM_PROCESS.outputEAgrd(outFile_fn = 'outTg_PER' + str(n+1) + '.asc')
            outFileRn = MM_PROCESS.outputEAgrd(outFile_fn = 'outRn_PER' + str(n+1) + '.asc')
            outFileETg = MM_PROCESS.outputEAgrd(outFile_fn = 'outETg_PER' + str(n+1) + '.asc')
            outFileR = MM_PROCESS.outputEAgrd(outFile_fn = 'outR_PER' + str(n+1) + '.asc')
            outFileSEEPAGE  = MM_PROCESS.outputEAgrd(outFile_fn = 'outSEEPAGE_PER' + str(n+1) + '.asc')
            outFileMB = MM_PROCESS.outputEAgrd(outFile_fn = 'outMB_PER' + str(n+1) + '.asc')
            outFileINTER = MM_PROCESS.outputEAgrd(outFile_fn = 'outINTER_PER' + str(n+1) + '.asc')
            outFileHEADS_MF = MM_PROCESS.outputEAgrd(outFile_fn = 'outHEADS_PER' + str(n+1) + '.asc')
        tstart = 0
        for t in range(n):
            tstart = tstart + perlen[t]
        tend = tstart + perlen[n]
        results=np.zeros([nrow,ncol,len(index),perlen[n]], dtype=float)
        results_S=np.zeros([nrow,ncol,len(index_S),_nslmax,perlen[n]], dtype=float)
        ncell = 0
        # loop into the grid
        for i in range(nrow):
            for j in range(ncol):
                SOILzone_tmp = gridSOIL[i,j]-1
                METEOzone_tmp = gridMETEO[i,j]-1
                if n == 0:
                    cbc[:,:,i,j,:] = conv_fact*cbc[:,:,i,j,:]/(delr[j]*delc[i])
                if ibound[i,j,0]<>0.0:
                    ncell = ncell + 1
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
                    PONDm_tmp = 1.12*gridPONDhmax[i,j]*gridPONDw[i,j]/delr[j]
                    PONDratio = 1.12*gridPONDw[i,j]/delr[j]
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
##                    if n==0:
##                        h_MF_tmp = np.zeros([perlen[n]], dtype = float)
##                        h_MF_tmp[0] = h_MF[0,i,j,0]
##                        h_MF_tmp[1:perlen[n]] = h_MF[tstart:tend-1,i,j,0]
##                        cbc_tmp = np.zeros([perlen[n]], dtype = float)
##                        cbc_tmp[0] = -cbc[0,iDRN,i,j,0]
##                        cbc_tmp[1:perlen[n]] = -cbc[tstart:tend-1,iDRN,i,j,0]
##                    else:
##                        h_MF_tmp = h_MF[tstart-1:tend-1,i,j,0]
##                        cbc_tmp = -cbc[tstart-1:tend-1,iDRN,i,j,0]
                    h_MF_tmp = h_MF[tstart:tend,i,j,0]
                    cbc_tmp = -cbc[tstart:tend,iDRN,i,j,0]
                    # for the while loop
                    h_MF_L0_m = np.ma.masked_values(h_MF_tmp, hdry, atol = 1E+25)
                    h_MF_L0_m_avg = np.nansum(h_MF_L0_m)/perlen[n]
                    if isinstance(h_MF_L0_m_avg, float):
                        h_MFsum = h_MFsum + h_MF_L0_m_avg
                    # cal functions for reservoirs calculations
                    results1_temp, results2_temp = MM_UNSAT.run(
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
                                                 ELEV    = top[i,j],
                                                 HEADS   = h_MF_tmp,
                                                 DRN     = cbc_tmp,
                                                 DRNi    = DRNi,
                                                 RF      = RFzonesTS[METEOzone_tmp][tstart:tend],
                                                 E0      = E0zonesTS[METEOzone_tmp][tstart:tend],
                                                 PETveg  = PETvegzonesTS_tmp,
                                                 RFeveg  = RFevegzonesTS_tmp,
                                                 PEsoil  = PEsoilzonesTS_tmp,
                                                 VEGarea = VEGarea_tmp,
                                                 Zr      = Zr,
                                                 perlen  = perlen[n],
                                                 hdry    = hdry,
                                                 k_Tu_slp = k_Tu_slp,
                                                 k_Tu_inter = k_Tu_inter)
                    for k in range(len(index)):
                        results[i,j,k,:] = results1_temp[k,:]
                        res_PERall[i,j,k,tstart:tend] = results1_temp[k,:]
                    for k in range(len(index_S)):
                        for l in range(nsl_tmp):
                           results_S[i,j,k,l,:] = results2_temp[k,l,:]
                           res_PERall_S[i,j,k,l,tstart:tend] = results2_temp[k,l,:]
                    # output cumulative (1) or averaged by time (sum(perlen))
                    divfact = float(perlen[n])
                    gridoutRF[i,j]      = results[i,j,index.get('iRF'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iRF')] = resavg_PERall[i,j,index.get('iRF')] + gridoutRF[i,j]
                    gridoutPET[i,j]     = results[i,j,index.get('iPET'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iPET')] = resavg_PERall[i,j,index.get('iPET')] + gridoutPET[i,j]
                    gridoutPE[i,j]     = results[i,j,index.get('iPE'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iPE')] = resavg_PERall[i,j,index.get('iPE')] + gridoutPE[i,j]
                    gridoutRFe[i,j]     = results[i,j,index.get('iRFe'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iRFe')] = resavg_PERall[i,j,index.get('iRFe')] + gridoutRFe[i,j]
                    gridoutdPOND[i,j]    = results[i,j,index.get('idPOND'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('idPOND')] = resavg_PERall[i,j,index.get('idPOND')] + gridoutdPOND[i,j]
                    gridoutPOND[i,j]    = results[i,j,index.get('iPOND'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iPOND')] = resavg_PERall[i,j,index.get('iPOND')] + gridoutPOND[i,j]
                    gridoutRo[i,j]      = results[i,j,index.get('iRo'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iRo')] = resavg_PERall[i,j,index.get('iRo')] + gridoutRo[i,j]
                    gridoutEs[i,j]     = results[i,j,index.get('iEs'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iEs')] = resavg_PERall[i,j,index.get('iEs')] + gridoutEs[i,j]
                    gridoutEg[i,j]       = results[i,j,index.get('iEg'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iEg')] = resavg_PERall[i,j,index.get('iEg')] + gridoutEg[i,j]
                    gridoutTg[i,j]       = results[i,j,index.get('iTg'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iTg')] = resavg_PERall[i,j,index.get('iTg')] + gridoutTg[i,j]
                    gridoutRn[i,j]       = results[i,j,index.get('iRn'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iRn')] = resavg_PERall[i,j,index.get('iRn')] + gridoutRn[i,j]
                    gridoutETg[i,j]       = results[i,j,index.get('iETg'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iETg')] = resavg_PERall[i,j,index.get('iETg')] + gridoutETg[i,j]
                    gridoutR[i,j]       = results[i,j,index.get('iR'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iR')] = resavg_PERall[i,j,index.get('iR')] + gridoutR[i,j]
                    gridoutSEEPAGE[i,j]   = results[i,j,index.get('iSEEPAGE'),:].sum()
                    resavg_PERall[i,j,index.get('iSEEPAGE')] = resavg_PERall[i,j,index.get('iSEEPAGE')] + gridoutSEEPAGE[i,j]
                    gridoutMB[i,j]      = results[i,j,index.get('iMB'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iMB')] = resavg_PERall[i,j,index.get('iMB')] + gridoutMB[i,j]
                    gridoutINTER[i,j]   = results[i,j,index.get('iINTER'),:].sum()/divfact
                    resavg_PERall[i,j,index.get('iINTER')] = resavg_PERall[i,j,index.get('iINTER')] + gridoutINTER[i,j]

                    for l in range(nsl_tmp):
                        gridoutEu[l][i,j] = results_S[i,j,index_S.get('iEu'),l,:].sum()/divfact
                        resavg_PERall_S[i,j,index_S.get('iEu'),l] = resavg_PERall_S[i,j,index_S.get('iEu'),l] + gridoutEu[l][i,j]
                        gridoutTu[l][i,j] = results_S[i,j,index_S.get('iTu'),l,:].sum()/divfact
                        resavg_PERall_S[i,j,index_S.get('iTu'),l] = resavg_PERall_S[i,j,index_S.get('iTu'),l] + gridoutTu[l][i,j]
                        gridoutS[l][i,j] = results_S[i,j,index_S.get('iS'),l,:].sum()/divfact
                        resavg_PERall_S[i,j,index_S.get('iSpc'),l] = resavg_PERall_S[i,j,index_S.get('iSpc'),l] + gridoutS[l][i,j]
                        gridoutRp[l][i,j] = results_S[i,j,index_S.get('iRp'),l,:].sum()/divfact
                        resavg_PERall_S[i,j,index_S.get('iRp'),l] = resavg_PERall_S[i,j,index_S.get('iRp'),l] + gridoutRp[l][i,j]

                    R4MF[n,i,j] = results[i,j,index.get('iR'),:].sum()/(divfact*1000)
                    ETg4MF[n,i,j] = (results[i,j,index.get('iETg'),:].sum()/(divfact*1000))*(delr[j]*delc[i])
                    # setting initial conditions for the next SP
                    for l in range(nsl_tmp):
                        Si_tmp_array[i,j,l] = results_S[i,j,index_S.get('iSpc'),l,perlen[n]-1]
                        Rpi_tmp_array[i,j,l] = results_S[i,j,index_S.get('iRp'),l,perlen[n]-1]
                        PONDi_tmp_array[i,j] = results[i,j,index.get('iPOND'),perlen[n]-1]
                        DRNi_tmp_array[i,j]  = cbc_tmp[-1]
                    if SP_d==0:
                        outFileRF.write('%.6f'%(gridoutRF[i,j]) + ' ')
                        outFilePET.write('%.6f'%(gridoutPET[i,j]) + ' ')
                        outFilePE.write('%.6f'%(gridoutPE[i,j]) + ' ')
                        outFileRFe.write('%.6f'%(gridoutRFe[i,j]) + ' ')
                        outFiledPOND.write('%.6f'%(gridoutdPOND[i,j]) + ' ')
                        outFilePOND.write('%.6f'%(gridoutPOND[i,j]) + ' ')
                        outFileRo.write('%.6f'%(gridoutRo[i,j]) + ' ')
                        for l in range(_nslmax):
                            outFileEu[l].write('%.6f'%(gridoutEu[l][i,j]) + ' ')
                            outFileTu[l].write('%.6f'%(gridoutTu[l][i,j]) + ' ')
                            outFileS[l].write('%.6f'%(gridoutS[l][i,j]) + ' ')
                            outFileRp[l].write('%.6f'%(gridoutRp[l][i,j]) + ' ')
                        outFileEs.write('%.6f'%(gridoutEs[i,j]) + ' ')
                        outFileEg.write('%.6f'%(gridoutEg[i,j]) + ' ')
                        outFileTg.write('%.6f'%(gridoutTg[i,j]) + ' ')
                        outFileRn.write('%.6f'%(gridoutRn[i,j]) + ' ')
                        outFileETg.write('%.6f'%(gridoutETg[i,j]) + ' ')
                        outFileR.write('%.6f'%(gridoutR[i,j]) + ' ')
                        outFileSEEPAGE.write('%.6f'%(gridoutSEEPAGE[i,j]) + ' ')
                        outFileMB.write('%.6f'%(gridoutMB[i,j]) + ' ')
                        outFileINTER.write('%.6f'%(gridoutINTER[i,j]) + ' ')
                        outFileHEADS_MF.write('%.6f'%(h_MF[tend-1,i,j,0]) + ' ')
                    if R4MF[n,i,j]>0.0:
                        outFileRPER.write(str(R4MF[n,i,j])+' ')
                    else:
                        outFileRPER.write(str(0.0)+' ')
                    if ETg4MF[n,i,j]<0.0:
                        outFileETgPER.write(str(ETg4MF[n,i,j])+' ')
                    else:
                        outFileETgPER.write(str(0.0)+' ')
                else:
                    results[i,j,:,:]    = hnoflo
                    results_S[i,j,:,:,:] = hnoflo
                    resavg_PERall[i,j,:] = hnoflo
                    resavg_PERall_S[i,j,:,:] = hnoflo
                    gridoutRF[i,j]      = hnoflo
                    gridoutPET[i,j]     = hnoflo
                    gridoutPE[i,j]      = hnoflo
                    gridoutRFe[i,j]     = hnoflo
                    gridoutdPOND[i,j]    = hnoflo
                    gridoutPOND[i,j]    = hnoflo
                    gridoutRo[i,j]      = hnoflo
                    for l in range(nsl_tmp):
                        gridoutEu[l][i,j]     = hnoflo
                        gridoutTu[l][i,j]     = hnoflo
                        gridoutS[l][i,j]       = hnoflo
                        gridoutRp[l][i,j]      = hnoflo
                    gridoutEs[i,j]      = hnoflo
                    gridoutEg[i,j]      = hnoflo
                    gridoutTg[i,j]      = hnoflo
                    gridoutRn[i,j]      = hnoflo
                    gridoutETg[i,j]      = hnoflo
                    gridoutR[i,j]       = hnoflo
                    gridoutSEEPAGE[i,j]   = hnoflo
                    gridoutMB[i,j]      = hnoflo
                    gridoutINTER[i,j]   = hnoflo
                    R4MF[n,i,j]        = hnoflo
                    ETg4MF[n,i,j]        = hnoflo
                    outFileRPER.write(str(hnoflo)+' ')
                    outFileETgPER.write(str(hnoflo)+' ')
                    if SP_d==0:
                        outFileRF.write('%.6f'%(gridoutRF[i,j]) + ' ')
                        outFilePET.write('%.6f'%(gridoutPET[i,j]) + ' ')
                        outFilePE.write('%.6f'%(gridoutPE[i,j]) + ' ')
                        outFileRFe.write('%.6f'%(gridoutRFe[i,j]) + ' ')
                        outFiledPOND.write('%.6f'%(gridoutdPOND[i,j]) + ' ')
                        outFilePOND.write('%.6f'%(gridoutPOND[i,j]) + ' ')
                        outFileRo.write('%.6f'%(gridoutRo[i,j]) + ' ')
                        for l in range(_nslmax):
                            outFileEu[l].write('%.6f'%(gridoutEu[l][i,j]) + ' ')
                            outFileTu[l].write('%.6f'%(gridoutTu[l][i,j]) + ' ')
                            outFileS[l].write('%.6f'%(gridoutS[l][i,j]) + ' ')
                            outFileRp[l].write('%.6f'%(gridoutRp[l][i,j]) + ' ')
                        outFileEs.write('%.6f'%(gridoutEs[i,j]) + ' ')
                        outFileEg.write('%.6f'%(gridoutEg[i,j]) + ' ')
                        outFileTg.write('%.6f'%(gridoutTg[i,j]) + ' ')
                        outFileRn.write('%.6f'%(gridoutRn[i,j]) + ' ')
                        outFileETg.write('%.6f'%(gridoutETg[i,j]) + ' ')
                        outFileR.write('%.6f'%(gridoutR[i,j]) + ' ')
                        outFileSEEPAGE.write('%.6f'%(gridoutSEEPAGE[i,j]) + ' ')
                        outFileMB.write('%.6f'%(gridoutMB[i,j]) + ' ')
                        outFileINTER.write('%.6f'%(gridoutINTER[i,j]) + ' ')
                        outFileHEADS_MF.write('%.6f'%(h_MF[tend-1,i,j,0]) + ' ')

            outFileRPER.write('\n')
            outFileETgPER.write('\n')
            if SP_d==0:
                outFileRF.write('\n')
                outFilePET.write('\n')
                outFilePE.write('\n')
                outFileRFe.write('\n')
                outFiledPOND.write('\n')
                outFilePOND.write('\n')
                outFileRo.write('\n')
                for l in range(_nslmax):
                    outFileEu[l].write('\n')
                    outFileTu[l].write('\n')
                    outFileS[l].write('\n')
                    outFileRp[l].write('\n')
                outFileEs.write('\n')
                outFileEg.write('\n')
                outFileTg.write('\n')
                outFileRn.write('\n')
                outFileETg.write('\n')
                outFileR.write('\n')
                outFileSEEPAGE.write('\n')
                outFileMB.write('\n')
                outFileINTER.write('\n')
                outFileHEADS_MF.write('\n')
        # close export ASCII files
        outFileRPER.close()
        outFileETgPER.close()
        if SP_d==0:
            outFileRF.close()
            outFilePET.close()
            outFilePE.close()
            outFileRFe.close()
            outFiledPOND.close()
            outFilePOND.close()
            outFileRo.close()
            for l in range(_nslmax):
                outFileEu[l].close()
                outFileTu[l].close()
                outFileS[l].close()
                outFileRp[l].close()
            outFileEs.close()
            outFileEg.close()
            outFileTg.close()
            outFileRn.close()
            outFileETg.close()
            outFileR.close()
            outFileSEEPAGE.close()
            outFileMB.close()
            outFileINTER.close()
            outFileHEADS_MF.close()

        t0=t0+int(perlen[n])
     #   print '\nSTRESS PERIOD %i/%i DONE!' % (n+1, nper)

    # #############################
    # ###  PRINT CONVERG. PLOT ###
    # #############################
    h_MFsum = h_MFsum/ncell/nper
    h_diff.append(h_MFsum - h_pSP)
    LOOP = LOOP + 1
    LOOPlst.append(LOOP)
    h_pSP = h_MFsum
    if pylab.absolute(h_diff[LOOP])>0.0:
        h_diff_log.append(pylab.log10(pylab.absolute(h_diff[LOOP])))
    else:
        h_diff_log.append(pylab.log10(convcrit))
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
    plt.close()
    del fig
    if verbose == 0:
        report.write(s.getvalue())
        s.close()
        s = StringIO.StringIO()
        sys.stdout = s
    if LOOP <2:
        print "\nInitial average heads:\n%.3f m" % h_diff[LOOP]
    else:
        print "\nHeads diff. from previous conv. loop:\n%.3f m" % h_diff[LOOP]
    if h_MFsum == 0.0:
        print '\nFirst layer of the model DRY!'
    elif abs(h_diff[LOOP]) < convcrit:
        print '\nSuccessfull convergence between MARMITES and MODFLOW!\n(Conv. criterion = %.3G)' % convcrit
        break
    elif LOOP>ccnum:
        print'\nNo convergence between MARMITES and MODFLOW!\n(Conv. criterion = %.3G)' % convcrit
        break

    # #############################
    # ### MODFLOW RUN with MM-computed recharge
    # #############################
    durationMF = 0.0
    timestartMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
    print'\n##############'
    print 'MODFLOW RUN (MARMITES fluxes)'
    SP_d, nrow, ncol, delr, delc, nlay, perlen, nper, top, hnoflo, hdry, ibound, laytyp, h_MF, cbc, cbc_nam_tmp, top_array, inputFileMF_fn, lenuni = ppMF.ppMF(MF_ws, rch_input = rch_input, rch_dft = 0.0001, wel_input = wel_input)

    timeendMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
    durationMF = durationMF + (timeendMF-timestartMF)

plt.close()
timeend = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
duration = duration + (timeend-timestart) -durationMF

# #############################
# ###  END CONVERGENCE LOOP ###
# #############################

# #############################
# ### MODFLOW RUN with MM-computed recharge
# #############################
durationMF = 0.0
timestartMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
print'\n##############'
print 'MODFLOW RUN (MARMITES fluxes after conv. loop)'
SP_d, nrow, ncol, delr, delc, nlay, perlen, nper, top, hnoflo, hdry, ibound, laytyp, h_MF, cbc, cbc_nam_tmp, top_array, inputFileMF_fn, lenuni = ppMF.ppMF(MF_ws, rch_input = rch_input, rch_dft = 0.0001, wel_input = wel_input)

h_MF_m = np.ma.masked_values(h_MF, hnoflo, atol = 0.09)
top_array_m = np.ma.masked_values(top_array, hnoflo, atol = 0.09)

for i in range(nrow):
    for j in range(ncol):
        cbc[:,:,i,j,:] = conv_fact*cbc[:,:,i,j,:]/(delr[j]*delc[i])

timeendMF = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
durationMF = durationMF + (timeendMF-timestartMF)

# #############################
# ###  OUTPUT EXPORT   ########
# #############################

timestartExport = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
print '\n##############\nMARMITES exporting...'

# write fluxes for the total time into ESRI ASCII grid
for i in range(nrow):
    for j in range(ncol):
        if ibound[i,j,0]<>0:
            outFileRF_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iRF')]/nper) + ' ')
            outFilePET_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iPET')]/nper) + ' ')
            outFilePE_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iPE')]/nper) + ' ')
            outFileRFe_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iRFe')]/nper) + ' ')
            outFiledPOND_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('idPOND')]/nper) + ' ')
            outFilePOND_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iPOND')]/nper) + ' ')
            outFileRo_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iRo')]/nper) + ' ')
            for l in range(_nslmax):
                outFileEu_PERall[l].write('%.6f'%(resavg_PERall_S[i,j,index_S.get('iEu'),l]/nper) + ' ')
                outFileTu_PERall[l].write('%.6f'%(resavg_PERall_S[i,j,index_S.get('iTu'),l]/nper) + ' ')
                outFileS_PERall[l].write('%.6f'%(resavg_PERall_S[i,j,index_S.get('iSpc'),l]/nper) + ' ')
                outFileRp_PERall[l].write('%.6f'%(resavg_PERall_S[i,j,index_S.get('iRp'),l]/nper) + ' ')
            outFileEs_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iEs')]/nper) + ' ')
            outFileEg_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iEg')]/nper) + ' ')
            outFileTg_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iTg')]/nper) + ' ')
            outFileRn_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iRn')]/nper) + ' ')
            outFileETg_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iETg')]/nper) + ' ')
            outFileR_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iR')]/nper) + ' ')
            outFileSEEPAGE_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iSEEPAGE')]) + ' ')
            outFileMB_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iMB')]) + ' ')
            outFileINTER_PERall.write('%.6f'%(resavg_PERall[i,j,index.get('iINTER')]/nper) + ' ')
        else:
            outFileRF_PERall.write('%.6f'%(hnoflo) + ' ')
            outFilePET_PERall.write('%.6f'%(hnoflo) + ' ')
            outFilePE_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileRFe_PERall.write('%.6f'%(hnoflo) + ' ')
            outFiledPOND_PERall.write('%.6f'%(hnoflo) + ' ')
            outFilePOND_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileRo_PERall.write('%.6f'%(hnoflo) + ' ')
            for l in range(_nslmax):
                outFileEu_PERall[l].write('%.6f'%(hnoflo) + ' ')
                outFileTu_PERall[l].write('%.6f'%(hnoflo) + ' ')
                outFileS_PERall[l].write('%.6f'%(hnoflo) + ' ')
                outFileRp_PERall[l].write('%.6f'%(hnoflo) + ' ')
            outFileEs_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileEg_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileTg_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileRn_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileETg_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileR_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileSEEPAGE_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileMB_PERall.write('%.6f'%(hnoflo) + ' ')
            outFileINTER_PERall.write('%.6f'%(hnoflo) + ' ')
    outFileRF_PERall.write('\n')
    outFilePET_PERall.write('\n')
    outFilePE_PERall.write('\n')
    outFileRFe_PERall.write('\n')
    outFiledPOND_PERall.write('\n')
    outFilePOND_PERall.write('\n')
    outFileRo_PERall.write('\n')
    for l in range(_nslmax):
        outFileEu_PERall[l].write('\n')
        outFileTu_PERall[l].write('\n')
        outFileS_PERall[l].write('\n')
        outFileRp_PERall[l].write('\n')
    outFileEs_PERall.write('\n')
    outFileEg_PERall.write('\n')
    outFileTg_PERall.write('\n')
    outFileRn_PERall.write('\n')
    outFileETg_PERall.write('\n')
    outFileR_PERall.write('\n')
    outFileSEEPAGE_PERall.write('\n')
    outFileMB_PERall.write('\n')
    outFileINTER_PERall.write('\n')
outFileRF_PERall.close()
outFilePET_PERall.close()
outFilePE_PERall.close()
outFileRFe_PERall.close()
outFiledPOND_PERall.close()
outFilePOND_PERall.close()
outFileRo_PERall.close()
for l in range(_nslmax):
    outFileEu_PERall[l].close()
    outFileTu_PERall[l].close()
    outFileS_PERall[l].close()
    outFileRp_PERall[l].close()
outFileEs_PERall.close()
outFileEg_PERall.close()
outFileTg_PERall.close()
outFileRn_PERall.close()
outFileETg_PERall.close()
outFileR_PERall.close()
outFileSEEPAGE_PERall.close()
outFileMB_PERall.close()
outFileINTER_PERall.close()

# exporting ASCII files at observations and plots

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
    DRNmax.append(np.nanmax(-cbc[:,iDRN,:,:,:]).flatten())
    DRNmin.append(np.nanmin(-cbc[:,iDRN,:,:,:]).flatten())
    cbcmax.append(np.nanmax(-cbc[:,:,:,:,:]).flatten())
    cbcmin.append(np.nanmin(-cbc[:,:,:,:,:]).flatten())
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
plt_export_fn = os.path.join(MM_ws, '00_UNSATandGWbudgets.png')
index_cbc = [iRCH, iSTO, iDRN, iWEL]
flxlbl = ['RF', 'INTER', 'SEEPAGE', 'dS', 'dPOND', 'Ro', 'Es', 'Eu', 'Tu', 'R', 'Rn', 'Eg', 'Tg', 'ETg']
flxlst =[res_PERall[:,:,index.get('iRF'),:].sum()/sum(perlen)/ncell,
        -res_PERall[:,:,index.get('iINTER'),:].sum()/sum(perlen)/ncell,
        res_PERall[:,:,index.get('iSEEPAGE'),:].sum()/sum(perlen)/ncell,
        res_PERall_S[:,:,index_S.get('idS'),:,:].sum()/sum(perlen)/ncell,
        res_PERall[:,:,index.get('idPOND'),:].sum()/sum(perlen)/ncell,
        -res_PERall[:,:,index.get('iRo'),:].sum()/sum(perlen)/ncell,
        -res_PERall[:,:,index.get('iEs'),:].sum()/sum(perlen)/ncell,
        -res_PERall_S[:,:,index_S.get('iEu'),:,:].sum()/sum(perlen)/ncell,
        -res_PERall_S[:,:,index_S.get('iTu'),:,:].sum()/sum(perlen)/ncell,
        -res_PERall[:,:,index.get('iR'),:].sum()/sum(perlen)/ncell,
        -res_PERall[:,:,index.get('iRn'),:].sum()/sum(perlen)/ncell,
        -res_PERall[:,:,index.get('iEg'),:].sum()/sum(perlen)/ncell,
        -res_PERall[:,:,index.get('iTg'),:].sum()/sum(perlen)/ncell,
        res_PERall[:,:,index.get('iETg'),:].sum()/sum(perlen)/ncell]
for l in range(nlay):
    for x in range(len(index_cbc)):
        flxlst.append(cbc[:,index_cbc[x],:,:,l].sum()/sum(perlen)/ncell)
        flxlbl.append('GW_' + cbc_nam[index_cbc[x]] + '_L' + str(l+1))
colors_flx = CreateColors.main(hi=0, hf=180, numbcolors = len(flxlbl))
MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = 'MARMITES and MODFLOW water flux budget for the whole catchment', fluxmax = cbcmax, fluxmin = cbcmin)

# export data for observation cells and show calib graphs
if obsCHECK == 1:
    colors_nsl = CreateColors.main(hi=00, hf=180, numbcolors = (_nslmax+1))
    for o in range(len(obs.keys())):
        i = obs.get(obs.keys()[o])['i']
        j = obs.get(obs.keys()[o])['j']
        # SATFLOW
        h_satflow = MM_SATFLOW.run(res_PERall[i,j,index.get('iR'),:], float(obs.get(obs.keys()[o])['hi']),float(obs.get(obs.keys()[o])['h0']),float(obs.get(obs.keys()[o])['RC']),float(obs.get(obs.keys()[o])['STO']))
##        # correct heads from MF (1 day delay)
##        h_MF_tmp = np.zeros([sum(perlen)], dtype = float)
##        h_MF_tmp[0] = h_MF[0,i,j,0]
##        h_MF_tmp[1:sum(perlen)] = h_MF[0:sum(perlen)-1,i,j,0]
        # export ASCII file at piezometers location
        #TODO extract heads at piezo location and not center of cell
        MM_PROCESS.ExportResults(i, j, inputDate, _nslmax, res_PERall, index, res_PERall_S, index_S, -cbc[:,iDRN,i,j,0], cbc[:,iRCH,i,j,0], -cbc[:,iWEL,i,j,0], h_satflow, h_MF[:,i,j,0], obs_h[o][0,:], obs_S[o], outFileExport[o], outPESTheads, outPESTsm, obs.keys()[o])
        outFileExport[o].close()
        # plot
        if plot_out == 1:
            # DateInput, P, PET, Pe, POND, dPOND, Ro, ETa, S, Rp, R, h, hmeas, Smeas, Sm, Sr):
            # index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iPOND':4, 'iRo':5, 'iSEEPAGE':6, 'iEs':7, 'iMB':8, 'iINTER':9, 'iE0':10, 'iEg':11, 'iTg':12, 'iR':13, 'iRn':14, 'idPOND':15, 'iETg':16}
            plt_export_fn = os.path.join(MM_ws, '00_'+ obs.keys()[o] + '.png')
            plt_title = obs.keys()[o]
            # def allPLOT(DateInput, P, PET, PE, Pe, dPOND, POND, Ro, Eu, Tu, Eg, Tg, S, dS, Spc, Rp, SEEPAGE, R, Rn, Es, MB, h_MF, h_SF, hmeas, Smeas, Sm, Sr, hnoflo, plt_export_fn, plt_title, colors_nsl, hmax, hmin):
            MMplot.allPLOT(
            inputDate,
            res_PERall[i,j,index.get('iRF'),:],
            res_PERall[i,j,index.get('iPET'),:],
            res_PERall[i,j,index.get('iPE'),:],
            res_PERall[i,j,index.get('iRFe'),:],
            res_PERall[i,j,index.get('idPOND'),:],
            res_PERall[i,j,index.get('iPOND'),:],
            res_PERall[i,j,index.get('iRo'),:],
            res_PERall_S[i,j,index_S.get('iEu'),0:_nsl[gridSOIL[i,j]-1],:],
            res_PERall_S[i,j,index_S.get('iTu'),0:_nsl[gridSOIL[i,j]-1],:],
            res_PERall[i,j,index.get('iEg'),:],
            res_PERall[i,j,index.get('iTg'),:],
            res_PERall_S[i,j,index_S.get('iS'),0:_nsl[gridSOIL[i,j]-1],:],
            res_PERall_S[i,j,index_S.get('idS'),0:_nsl[gridSOIL[i,j]-1],:],
            res_PERall_S[i,j,index_S.get('iSpc'),0:_nsl[gridSOIL[i,j]-1],:],
            res_PERall_S[i,j,index_S.get('iRp'),0:_nsl[gridSOIL[i,j]-1]-1,:],
            res_PERall[i,j,index.get('iSEEPAGE'),:],
            res_PERall[i,j,index.get('iR'),:],
            res_PERall[i,j,index.get('iRn'),:],
            res_PERall[i,j,index.get('iEs'),:],
            res_PERall[i,j,index.get('iMB'),:],
            h_MF[:,i,j,0], h_satflow, obs_h[o][0,:], obs_S[o],
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
            plt_export_fn = os.path.join(MM_ws, '00_'+ obs.keys()[o] + 'UNSATandGWbudgets.png')
            flxlst =[res_PERall[i,j,index.get('iRF'),:].sum()/sum(perlen),
                    -res_PERall[i,j,index.get('iINTER'),:].sum()/sum(perlen),
                    res_PERall[i,j,index.get('iSEEPAGE'),:].sum()/sum(perlen),
                    res_PERall_S[i,j,index_S.get('idS'),:,:].sum()/sum(perlen),
                    res_PERall[i,j,index.get('idPOND'),:].sum()/sum(perlen),
                    -res_PERall[i,j,index.get('iRo'),:].sum()/sum(perlen),
                    -res_PERall[i,j,index.get('iEs'),:].sum()/sum(perlen),
                    -res_PERall_S[i,j,index_S.get('iEu'),:,:].sum()/sum(perlen),
                    -res_PERall_S[i,j,index_S.get('iTu'),:,:].sum()/sum(perlen),
                    -res_PERall[i,j,index.get('iR'),:].sum()/sum(perlen),
                    -res_PERall[i,j,index.get('iRn'),:].sum()/sum(perlen),
                    -res_PERall[i,j,index.get('iEg'),:].sum()/sum(perlen),
                    -res_PERall[i,j,index.get('iTg'),:].sum()/sum(perlen),
                    res_PERall[i,j,index.get('iETg'),:].sum()/sum(perlen)]
            for l in range(nlay):
                for x in range(len(index_cbc)):
                    flxlst.append(cbc[:,index_cbc[x],i,j,l].sum()/sum(perlen))
            MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title, fluxmax = cbcmax, fluxmin = cbcmin)
    # output for PEST
    outPESTheads.close()
    outPESTsm.close()

if plot_out == 1:
    # plot heads (grid + contours), DRN, etc... at specified TS
    TSlst = []
    TS = 0
    while TS < len(h_MF):
        TSlst.append(TS)
        TS = TS + plot_freq
    TSlst.append(len(h_MF)-1)
    for TS in TSlst:
        # plot heads [m]
        V=[]
        for L in range(nlay):
            V.append(h_MF_m[TS,:,:,L])
        MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation (m)', msg = 'DRY', plt_title = 'HEADS', MM_ws = MM_ws, interval_type = 'arange', interval_diff = 0.5, Vmax = hmax, Vmin = hmin)
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
            V.append(-cbc[TS,iDRN,:,:,L])
        MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage (mm/time step)', msg = '- no drainage', plt_title = 'DRN', MM_ws = MM_ws, interval_type = 'linspace', interval_num = 5, Vmin = DRNmin, Vmax = DRNmax, fmt='%.3G')

timeendExport = pylab.datestr2num(pylab.datetime.datetime.today().isoformat())
durationExport=(timeendExport-timestartExport)

# final report of successful run

print ('\n##############\nMARMITES executed successfully!')
print ('%s time steps\n%sx%s cells (rows x cols)') % (int(sum(perlen)),str(nrow),str(ncol))
print ('\nMARMITES run time: %s minute(s) and %.1f second(s)') % (str(int(duration*24.0*60.0)), (duration*24.0*60.0-int(duration*24.0*60.0))*60)
print ('MODFLOW run time: %s minute(s) and %.1f second(s)') % (str(int(durationMF*24.0*60.0)), (durationMF*24.0*60.0-int(durationMF*24.0*60.0))*60)
print ('Export run time: %s minute(s) and %.1f second(s)') % (str(int(durationExport*24.0*60.0)), (durationExport*24.0*60.0-int(durationExport*24.0*60.0))*60)
print ('\nOutput written in folder: \n%s\n##############\n') % MM_ws

if verbose == 0:
    report.write(s.getvalue())
    s.close()
    report.close()
    # print '\nMARMITES succesfully run!\Check report %s' % (report_fn)
##except (ValueError, TypeError, KeyboardInterrupt, IOError), e:
##    print e
##    raise e

#os.system('pause')