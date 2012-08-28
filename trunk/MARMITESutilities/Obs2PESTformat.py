# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:        obs2PESTformat.py
# Purpose:
#
# Author:      frances08512
#
# Created:     29-06-2012
# Copyright:   (c) frances08512 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.3"
__date__ = "2012"

import sys, os, traceback, h5py
import numpy as np
import MARMITESprocess_v3 as MMproc
import ppMODFLOW_flopy_v3 as ppMF

#####################################

# workspace (ws) definition
# read input file (called _input.ini in the MARMITES workspace
# the first character on the first line has to be the character used to comment
# the file can contain any comments as the user wish, but the sequence of the input has to be respected
# 00_TESTS\MARMITESv3_r13c6l2'  00_TESTS\r40c20'  00_TESTS\r20c40'
# SARDON'  CARRIZAL' LAMATA'  LaMata_new'
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
    #run MARMITESsoil  (1 is YES, 0 is NO)
    MMsoil_yn = int(inputFile[l].strip())
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
    gridSshmax_fn =  inputFile[l].strip()
    l += 1
    gridSsw_fn =  inputFile[l].strip()
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
    if MMsoil_yn == 1 and MF_yn != 1:
        MF_yn == 1
    if MMsurf_plot == 1:
        plt_out = 0
        plt_out_obs = 0
        print "\nYou required the MMsurf plots to appear on the screen. Due to backends limitation, MM and MF plots were disabled. Run again MM with MMsurf_plot = 0 to obtain the MM and MF plots."
except:
    raise SystemExit('\nFATAL ERROR!\nType error in the input file %s' % (MM_fn))
del inputFile

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

try:
# #############################
# ###  MARMITES surface  ######
# #############################
    inputFile = MMproc.readFile(MM_ws,outMMsurf_fn)
    l=0
    Zr = []
    kTu_min = []
    kTu_n = []
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
    raise SystemExit('\nFATAL ERROR!\Error in reading file [' + outMMsurf_fn + '], run MMsurf first!')
del inputFile
numDays = len(MMproc.readFile(MM_ws, inputDate_fn))

try:
    # #############################
    # ###  READ MODFLOW CONFIG ####
    # #############################
    print'\n##############'
    print 'MODFLOW initialization'
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
        mask_tmp += (np.asarray(cMF.ibound)[:,:,l] <> 0)
    maskAllL = (mask_tmp == 0)
    del iboundBOL, mask_tmp

    # SOIL PARAMETERS
    _nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _Su_ini, _Ks = cMF.MM_PROCESS.inputSoilParam(MM_ws = MM_ws, SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
    _nslmax = max(_nsl)

    # #############################
    # ### MF time processing
    # #############################
    # if required by user, compute nper, perlen,etc based on RF analysis in the METEO zones
    if cMF.timedef >= 0:
        if isinstance(cMF.nper, str):
            try:
                perlenmax = int(cMF.nper.split()[1].strip())
            except:
                raise SystemExit('\nFATAL ERROR!\nError in nper format of the MODFLOW ini file!\n')
        if irr_yn == 0:
            cMF.ppMFtime(MM_ws, MF_ws, inputDate_fn, inputZON_dSP_RF_veg_fn, inputZON_dSP_RFe_veg_fn, inputZON_dSP_PT_fn,input_dSP_LAI_veg_fn, inputZON_dSP_PE_fn, inputZON_dSP_E0_fn, NMETEO, NVEG, NSOIL)
        else:
            cMF.ppMFtime(MM_ws, MF_ws, inputDate_fn, inputZON_dSP_RF_veg_fn, inputZON_dSP_RFe_veg_fn, inputZON_dSP_PT_fn, input_dSP_LAI_veg_fn, inputZON_dSP_PE_fn, inputZON_dSP_E0_fn, NMETEO, NVEG, NSOIL, inputZON_dSP_RF_irr_fn, inputZON_dSP_RFe_irr_fn, inputZON_dSP_PT_irr_fn, input_dSP_crop_irr_fn, NFIELD)

    # READ observations time series (heads and soil moisture)
    print "\nReading observations time series (hydraulic heads and soil moisture)..."
    obs, obs_list = cMF.MM_PROCESS.inputObs(MM_ws            = MM_ws,
                                  inputObs_fn      = inputObs_fn,
                                  inputObsHEADS_fn = inputObsHEADS_fn,
                                  inputObsSM_fn    = inputObsSM_fn,
                                  inputDate        = cMF.inputDate,
                                  _nslmax          = _nslmax,
                                  nlay             = cMF.nlay
                                  )
    outPESTheads_fn      = 'h_obs4PEST.smp'
    outPESTsm_fn         = 'sm_obs4PEST.smp'
    outPESTheads = open(os.path.join(MF_ws,outPESTheads_fn), 'w')
    outPESTsm = open(os.path.join(MF_ws,outPESTsm_fn), 'w')

    # exporting MM time series results to ASCII files and plots at observations cells
    print "\nExporting observations time series (hydraulic heads and soil moisture) in PEST format (smp file)..."
    for o_ref in obs_list:
        for o in obs.keys():
            if o == o_ref:
                i = obs.get(o)['i']
                j = obs.get(o)['j']
                cMF.MM_PROCESS.smMMname.append(o)
                obs_h = obs.get(o)['obs_h']
                obs_S = obs.get(o)['obs_S']
                # export ASCII file at piezometers location
                #TODO extract heads at piezo location and not center of cell
                if obs_h != []:
                    obs_h_tmp = obs_h[0,:]
                else:
                    obs_h_tmp = []
                # Export time series observations as ASCII file for PEST
                cMF.MM_PROCESS.ExportResultsPEST(i, j, cMF.inputDate, _nslmax, obs_h_tmp, obs_S, outPESTheads, outPESTsm, o)
                del obs_h, obs_S
    del obs
    outPESTheads.close()
    # write PEST obs smp file
    inputFile = MMproc.readFile(MM_ws,inputObs_fn)
    ind = 0
    for i in range(len(inputFile)):
        line = inputFile[i].split()
        name = line[0]
        for j in cMF.MM_PROCESS.smMMname:
            if j == name:
                for l in cMF.MM_PROCESS.smMM[ind]:
                    outPESTsm.write(l)
        ind += 1
    outPESTsm.close()
    print "\nDONE!\nFiles exported:\n%s\n%s" % (os.path.join(MF_ws,outPESTheads_fn),os.path.join(MF_ws,outPESTsm_fn))
    print ''

except StandardError, e:  #Exception
    raise SystemExit('\nFATAL ERROR!\nMM run interruption in the export phase!\nError description:\n%s' % traceback.print_exc(file=sys.stdout))
#    traceback.print_exc(limit=1, file=sys.stdout)