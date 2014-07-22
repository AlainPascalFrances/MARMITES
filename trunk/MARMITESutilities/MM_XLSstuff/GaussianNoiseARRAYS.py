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
import startMARMITESsurface as startMMsurf
import MARMITESunsat_v2 as MMunsat
import MARMITESprocess as MMproc
import ppMODFLOW_flopy as ppMF
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
MM_ws = r'E:\00code_ws\00_TESTS\MARMITESv2_r13c6l2_REF'
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

# #############################
# ###  READ MODFLOW CONFIG ####
# #############################

print'\n##############'
print 'Importing MODFLOW configuration file'
nrow, ncol, delr, delc, reggrid, nlay, nper, perlen, nstp, hnoflo, hdry, laytyp, lenuni, itmuni = ppMF.ppMFini(MF_ws, MF_ini_fn, out = 'MM')

modelname, namefile_ext, exe_name, dum_sssp1, ext_dis, nlay, ncol, nrow, nper, itmuni, lenuni,laycbd, delr, delc, top_fn, botm_fn, perlen, nstp, tsmult, Ss_tr, ext_bas, ibound_fn, strt_fn, hnoflo,ext_lpf, ilpfcb, hdry, nplpf, laytyp, layavg, chani, layvka, laywet, hk_fn, vka_fn, ss_fn, sy_fn,ext_oc, ihedfm, iddnfm, ext_cbc, ext_heads, ext_ddn, ext_rch, rch_input_user, rch_dft, nrchop, ext_wel, wel_input_user, wel_dft, ext_drn, drn_elev_fn, drn_cond_fn = ppMF.ppMFini(MF_ws, MF_ini_fn, out = 'MF')

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

# 1 - reaf asc file and convert in np.array

print "\nImporting ESRI ASCII files to initialize the MODFLOW packages..."

for e, (ws, fn, treshold) in enumerate(zip((MF_ws,MF_ws,MF_ws,MM_ws),(hk_fn,ss_fn,sy_fn,[gridSOILthick_fn]), (0.0,0.0,0.0,0.0))):
    for l in fn:
        path = os.path.join(ws, l)
        array_in = np.zeros((nrow,ncol))
        array_out = np.zeros((nrow,ncol))
        array_in[:,:] = MM_PROCESS.convASCIIraster2array(path, array_in[:,:])
        path = path.partition('.')[0] + '_RAND.asc'
        outFile=open(path, 'w')
        outFile=MM_PROCESS.writeHeaderESRIraster(outFile)
        for i in range(nrow):
            for j in range(ncol):
                 if array_in[i,j]<>treshold:
                    noise=np.random.normal(0,np.std(array_in[:,:])/4)  #mean, std dev, num pts
                    array_out[i,j] = array_in[i,j] + noise
                 outFile.write(str(array_out[i,j])+ '  ')
            outFile.write('\n')
        outFile.close()

print '\nDone!'