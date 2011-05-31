# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:        startMARMITES
# Purpose:
#
# Author:      frances08512
#
# Created:     01-03-2011
# Copyright:   (c) frances08512 2010
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.1"
__date__ = "March 2011"

import pylab
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import mf
import mfreadbinaries as mfrdbin
import h5py
import MARMITESprocess as MMproc

#####################################

def ppMFini(MF_ws, MF_ini_fn, out = 'MF'):

    ''' read input file (called _input.ini in the MARMITES workspace
    the first character on the first line has to be the character used to comment
    the file can contain any comments as the user wish, but the sequence of the input has to be respected'''

    inputFile = MMproc.readFile(MF_ws, MF_ini_fn)
    l=0
    try:
        modelname =  str(inputFile[l].strip())
        l += 1
        namefile_ext = str(inputFile[l].strip())
        l += 1
        exe_name = str(inputFile[l].strip())
        l += 1
        dum_sssp1 = int(inputFile[l].strip())
        l += 1
        ext_dis = str(inputFile[l].strip())
        l += 1
        nlay = int(inputFile[l].strip())
        l += 1
        ncol = int(inputFile[l].strip())
        l += 1
        nrow = int(inputFile[l].strip())
        l += 1
        nper = inputFile[l].strip()
        try:
            nper = int(nper)
        except:
            pass
        if isinstance(nper, int):
            if nper <0:
            # daily nper, perlen, nstp, tsmult
                nper = abs(nper)
                timedef = -1
            else:
                # read nper, perlen, nstp, tsmult from user-input ini file
                nper = nper
                timedef = 1
        elif isinstance(nper, str):
            # compute nper, perlen, nstp, tsmult based on daily RF analysis
            nper = nper[0]
            timedef = 0
            nper = inputFile[l].strip()
        l += 1
        itmuni = int(inputFile[l].strip())
        l += 1
        lenuni = int(inputFile[l].strip())
        l += 1
        laycbd = []
        laycbd_tmp =  inputFile[l].split()
        for i in range(nlay):
            laycbd.append(int(laycbd_tmp[i]))
        l += 1
        delr = []
        delr_tmp =  inputFile[l].split()
        if len(delr_tmp)>1:
            for i in range(ncol):
                delr.append(int(delr_tmp[i]))
        else:
            for i in range(ncol):
                delr.append(int(delr_tmp[0]))
        l += 1
        delc = []
        delc_tmp =  inputFile[l].split()
        if len(delc_tmp)>1:
            for i in range(nrow):
                delc.append(int(delc_tmp[i]))
        else:
            for i in range(nrow):
                delc.append(int(delc_tmp[0]))
        l += 1
        reggrid = int(inputFile[l].strip())
        l += 1
        top_fn = str(inputFile[l].strip())
        l += 1
        botm_fn = []
        botm_tmp =  inputFile[l].split()
        for i in range(nlay):
            botm_fn.append(str(botm_tmp[i]))
        l += 1
        perlen = []
        nstp = []
        tsmult = []
        Ss_tr = []
        if timedef > 0:
            perlen_tmp =  inputFile[l].split()
            for i in range(nper):
                perlen.append(int(perlen_tmp[i]))
            l += 1
            nstp_tmp =  inputFile[l].split()
            for i in range(nper):
                nstp.append(int(nstp_tmp[i]))
            l += 1
            tsmult_tmp =  inputFile[l].split()
            for i in range(nper):
                tsmult.append(int(tsmult_tmp[i]))
            l += 1
            Ss_tr_tmp =  inputFile[l].split()
            for i in range(nper):
                Ss_tr.append(str(Ss_tr_tmp[i]))
            l += 1
        elif timedef < 0:
            for t in range(nper):
                perlen.append(1)
                nstp.append(1)
                tsmult.append(1)
                Ss_tr.append('TR')
            l += 4
        elif timedef == 0:
            perlen = 'ToBeDefined'
            nstp = 'ToBeDefined'
            tsmult = 'ToBeDefined'
            l += 4
        ext_bas = str(inputFile[l].strip())
        l += 1
        ibound_fn = []
        ibound_tmp =  inputFile[l].split()
        for i in range(nlay):
            ibound_fn.append(str(ibound_tmp[i]))
        l += 1
        strt_fn = []
        strt_tmp =  inputFile[l].split()
        for i in range(nlay):
            strt_fn.append(str(strt_tmp[i]))
        l += 1
        hnoflo = float(inputFile[l].strip())
        l += 1
        ext_lpf = str(inputFile[l].strip())
        l += 1
        ilpfcb = int(inputFile[l].strip())
        l += 1
        hdry = eval(inputFile[l].strip())
        l += 1
        nplpf = int(inputFile[l].strip())
        l += 1
        laytyp = []
        laytyp_tmp =  inputFile[l].split()
        for i in range(nlay):
            laytyp.append(int(laytyp_tmp[i]))
        l += 1
        layavg = []
        layavg_tmp =  inputFile[l].split()
        for i in range(nlay):
            layavg.append(int(layavg_tmp[i]))
        l += 1
        chani = []
        chani_tmp =  inputFile[l].split()
        for i in range(nlay):
            chani.append(int(chani_tmp[i]))
        l += 1
        layvka = []
        layvka_tmp =  inputFile[l].split()
        for i in range(nlay):
            layvka.append(int(layvka_tmp[i]))
        l += 1
        laywet = []
        laywet_tmp =  inputFile[l].split()
        for i in range(nlay):
            laywet.append(int(laywet_tmp[i]))
        l += 1
        hk_fn = []
        hk_tmp =  inputFile[l].split()
        for i in range(nlay):
            hk_fn.append(str(hk_tmp[i]))
        l += 1
        vka_fn = []
        vka_tmp =  inputFile[l].split()
        for i in range(nlay):
            vka_fn.append(str(vka_tmp[i]))
        l += 1
        ss_fn = []
        ss_tmp =  inputFile[l].split()
        for i in range(nlay):
            ss_fn.append(str(ss_tmp[i]))
        l += 1
        sy_fn = []
        sy_tmp =  inputFile[l].split()
        for i in range(nlay):
            sy_fn.append(str(sy_tmp[i]))
        l += 1
        ext_oc = str(inputFile[l].strip())
        l += 1
        ihedfm = int(inputFile[l].strip())
        l += 1
        iddnfm = int(inputFile[l].strip())
        l += 1
        ext_cbc = str(inputFile[l].strip())
        l += 1
        ext_heads = str(inputFile[l].strip())
        l += 1
        ext_ddn = str(inputFile[l].strip())
        l += 1
        ext_rch = str(inputFile[l].strip())
        l += 1
        nrchop = int(inputFile[l].strip())
        l += 1
        rch_user = float(inputFile[l].strip())
        l += 1
        ext_wel = str(inputFile[l].strip())
        l += 1
        wel_user = float(inputFile[l].strip())
        l += 1
        ext_drn = str(inputFile[l].strip())
        l += 1
        drn_elev_fn = []
        drn_elev_tmp =  inputFile[l].split()
        for i in range(nlay-1):
            drn_elev_fn.append(str(drn_elev_tmp[i]))
        l += 1
        drn_cond_fn = []
        drn_cond_tmp =  inputFile[l].split()
        for i in range(nlay-1):
            drn_cond_fn.append(str(drn_cond_tmp[i]))
    except:
        print "Unexpected error in the input file:\n", sys.exc_info()[0]
        sys.exit()
    del inputFile

    if out == 'MF':
        return modelname, namefile_ext, exe_name, dum_sssp1, ext_dis, nlay, ncol, nrow, nper, itmuni, lenuni,laycbd, delr, delc, top_fn, botm_fn, perlen, nstp, tsmult, Ss_tr, ext_bas, ibound_fn, strt_fn, hnoflo,ext_lpf, ilpfcb, hdry, nplpf, laytyp, layavg, chani, layvka, laywet, hk_fn, vka_fn, ss_fn, sy_fn,ext_oc, ihedfm, iddnfm, ext_cbc, ext_heads, ext_ddn, ext_rch, rch_user, nrchop, ext_wel, wel_user, ext_drn, drn_elev_fn, drn_cond_fn
    elif out == 'MM':
        return nrow, ncol, delr, delc, reggrid, nlay, nper, perlen, nstp, hnoflo, hdry, laytyp, lenuni, itmuni, ibound_fn

#####################################

def ppMFtime(MM_ws, MF_ws, MFtime_fn, perlenmax, inputDate_fn, inputZON_TS_RF_fn, inputZON_TS_PET_fn, inputZON_TS_RFe_fn, inputZON_TS_PE_fn, inputZON_TS_E0_fn, NMETEO, NVEG, NSOIL):

    ''' RF analysis to define automatically nper/perlen/nstp
    Daily RF>0 creates a nper
    Succesive days with RF=0 are averaged (PET, PE, E0) to a user-input maximum # of days'''

    def ExportResults1(TS, outFileExport):
        """
        Write the processed data in a open txt file readable by  MARMITES
        INPUT:      output flux time series and open file
        """
        for i in range(len(TS)):
            out_line =  '%14.9G' %TS[i],'\n'
            for l in out_line:
                outFileExport.write(l)

    #####################################

    print'\nComputing MODFLOW time discretization based on rainfall analysis in the METEO zones.'

    inputFile = MMproc.readFile(MM_ws, inputDate_fn)
    d = []
    for l in inputFile:
        d.append(l)
    inputFileRF = MMproc.readFile(MM_ws, inputZON_TS_RF_fn)
    RF_d = np.zeros([NMETEO, len(d)])
    inputFilePET = MMproc.readFile(MM_ws, inputZON_TS_PET_fn)
    inputFileRFe = MMproc.readFile(MM_ws, inputZON_TS_RFe_fn)
    PET_d = np.zeros([NMETEO, NVEG, len(d)])
    RFe_d = np.zeros([NMETEO, NVEG, len(d)])
    inputFilePE = MMproc.readFile(MM_ws, inputZON_TS_PE_fn)
    PE_d = np.zeros([NMETEO, NSOIL, len(d)])
    inputFileE0 = MMproc.readFile(MM_ws, inputZON_TS_E0_fn)
    E0_d = np.zeros([NMETEO, len(d)])
    for n in range(NMETEO):
        for t in range(len(d)):
            RF_d[n,t] = float(inputFileRF[t+len(d)*n].strip())
            E0_d[n,t] = float(inputFileE0[t+len(d)*n].strip())
        for v in range(NVEG):
            for t in range(len(d)):
                    PET_d[n,v,t] = float(inputFilePET[t+(n*NVEG+v)*len(d)].strip())
                    RFe_d[n,v,t] = float(inputFileRFe[t+(n*NVEG+v)*len(d)].strip())
        for s in range(NSOIL):
            for t in range(len(d)):
                PE_d[n,s,t] = float(inputFilePE[t+(n*NSOIL+s)*len(d)].strip())
    del inputFileRF, inputFilePET, inputFilePE, inputFileRFe, inputFileE0

    RF_stp = []
    PET_stp = []
    RFe_stp = []
    PE_stp = []
    E0_stp = []
    PET_d_tmp = []
    PE_d_tmp=[]
    E0_d_tmp = []
    nper = 1
    perlen = [1]
    c = 0
    for n in range(NMETEO):
        RF_stp.append([])
        PET_stp.append([])
        RFe_stp.append([])
        PET_d_tmp.append([])
        for v in range(NVEG):
            PET_stp[n].append([])
            RFe_stp[n].append([])
            PET_d_tmp[n].append(0.0)
        PE_stp.append([])
        PE_d_tmp.append([])
        for v in range(NSOIL):
            PE_stp[n].append([])
            PE_d_tmp[n].append(0.0)
        E0_stp.append([])
        E0_d_tmp.append(0.0)
    if perlenmax < 2:
        print '\nperlenmax must be higher than 1!\nCorrect perlenmax in the MODFLOW ini file.'
        sys.exit()
    for j in range(len(d)):
            if RF_d[:,j].sum()>0.0:
                if c == 1:
                    for n in range(NMETEO):
                        RF_stp[n].append(0.0)
                        for v in range(NVEG):
                            PET_stp[n][v].append(PET_d_tmp[n][v]/perlen[nper-1])
                            RFe_stp[n][v].append(0.0)
                            PET_d_tmp[n][v] = 0.0
                        for s in range(NSOIL):
                            PE_stp[n][s].append(PE_d_tmp[n][s]/perlen[nper-1])
                            PE_d_tmp[n][s] = 0.0
                        E0_stp[n].append(E0_d_tmp[n]/perlen[nper-1])
                        E0_d_tmp[n] = 0.0
                for n in range(NMETEO):
                    RF_stp[n].append(RF_d[n,j])
                    for v in range(NVEG):
                        PET_stp[n][v].append(PET_d[n,v,j])
                        RFe_stp[n][v].append(RFe_d[n,v,j])
                    for s in range(NSOIL):
                        PE_stp[n][s].append(PE_d[n,s,j])
                    E0_stp[n].append(E0_d[n,j])
                if j>0:
                    nper += 1
                    if j < len(d):
                        perlen.append(1)
                c = 0
            else:
                if perlen[nper-1] < perlenmax:
                    for n in range(NMETEO):
                        for v in range(NVEG):
                            PET_d_tmp[n][v] += PET_d[n,v,j]
                        for s in range(NSOIL):
                            PE_d_tmp[n][s]  += PE_d[n,s,j]
                        E0_d_tmp[n]  += E0_d[n,j]
                    if j>0:
                        if c == 0:
                            nper += 1
                            if j < len(d):
                                perlen.append(1)
                        elif c == 1:
                            perlen[nper-1] += 1
                    c = 1
                else:
                    for n in range(NMETEO):
                        RF_stp[n].append(0.0)
                        for v in range(NVEG):
                            PET_stp[n][v].append(PET_d_tmp[n][v]/perlen[nper-1])
                            RFe_stp[n][v].append(0.0)
                            PET_d_tmp[n][v] = 0.0
                        for s in range(NSOIL):
                            PE_stp[n][s].append(PE_d_tmp[n][s]/perlen[nper-1])
                            PE_d_tmp[n][s] = 0.0
                        E0_stp[n].append(E0_d_tmp[n]/perlen[nper-1])
                        E0_d_tmp[n] = 0.0
                    if j>0:
                        nper += 1
                        if j < len(d):
                            perlen.append(1)
                    c = 1
    del RF_d, RFe_d, PET_d, PE_d, E0_d, PET_d_tmp, PE_d_tmp, E0_d_tmp, c
    perlen = np.asarray(perlen, dtype = np.int)
    nstp = np.ones(nper, dtype = np.int)
    tsmult = nstp
    Ss_tr = []
    for n in range(nper):
        Ss_tr.append('TR')

    outpathname = os.path.join(MF_ws, MFtime_fn)
    outFileExport = open(outpathname, 'w')
    outFileExport.write('#')
    outFileExport.write('\n# nper\n')
    outFileExport.write(str(nper))
    outFileExport.write('\n# perlen nstp tsmult Ss_tr')
    for n in range(nper):
        outFileExport.write('\n' + str(perlen[n]) + ' ')
        outFileExport.write(str(nstp[n]) + ' ')
        outFileExport.write(str(tsmult[n]) + ' ')
        outFileExport.write(str(Ss_tr[n]))
    outFileExport.close()

    inputZON_TS_RF_fn = "inputZONRF_stp.txt"
    inputZONRF_fn = os.path.join(MM_ws, inputZON_TS_RF_fn)
    inputZONRF = open(inputZONRF_fn, 'w')
    inputZONRF.write('#\n')

    inputZON_TS_PET_fn = "inputZONPET_stp.txt"
    inputZONPET_fn = os.path.join(MM_ws, inputZON_TS_PET_fn)
    inputZONPET = open(inputZONPET_fn, 'w')
    inputZONPET.write('#\n')

    inputZON_TS_RFe_fn = "inputZONRFe_stp.txt"
    inputZONRFe_fn = os.path.join(MM_ws, inputZON_TS_RFe_fn)
    inputZONRFe = open(inputZONRFe_fn, 'w')
    inputZONRFe.write('#\n')

    inputZON_TS_PE_fn = "inputZONPE_stp.txt"
    inputZONPE_fn = os.path.join(MM_ws, inputZON_TS_PE_fn)
    inputZONPE = open(inputZONPE_fn, 'w')
    inputZONPE.write('#\n')

    inputZON_TS_E0_fn = "inputZONE0_stp.txt"
    inputZONE0_fn = os.path.join(MM_ws, inputZON_TS_E0_fn)
    inputZONE0 = open(inputZONE0_fn, 'w')
    inputZONE0.write('#\n')

    for n in range(NMETEO):
        try:
            ExportResults1(RF_stp[n], inputZONRF)
            ExportResults1(E0_stp[n], inputZONE0)
            #VEG
            if NVEG>0:
                for v in range(NVEG):
                    ExportResults1(PET_stp[n][v], inputZONPET)
                    ExportResults1(RFe_stp[n][v], inputZONRFe)
            # SOIL
            if NSOIL>0:
                for s in range(NSOIL):
                    ExportResults1(PE_stp[n][s], inputZONPE)
        except:
            print "\nError in output file, some output files were not exported."

    inputZONRF.close()
    inputZONRFe.close()
    inputZONPET.close()
    inputZONPE.close()
    inputZONE0.close()

    return nper, perlen, nstp, inputZON_TS_RF_fn, inputZON_TS_PET_fn, inputZON_TS_RFe_fn, inputZON_TS_PE_fn, inputZON_TS_E0_fn

#####################################

def ppMF(MM_ws, xllcorner, yllcorner, MF_ws, MF_ini_fn, rch_MM = "", rch_user = None, wel_MM = "", wel_user = None, report = None, verbose = 1, chunks = 0):

    messagemanual="Please read the manual!\n(that by the way still doesn't exist...)"

    if verbose == 0:
        print '--------------'

    modelname, namefile_ext, exe_name, dum_sssp1, ext_dis, nlay, ncol, nrow, nper, itmuni, lenuni,laycbd, delr, delc, top_fn, botm_fn, perlen, nstp, tsmult, Ss_tr, ext_bas, ibound_fn, strt_fn, hnoflo,ext_lpf, ilpfcb, hdry, nplpf, laytyp, layavg, chani, layvka, laywet, hk_fn, vka_fn, ss_fn, sy_fn,ext_oc, ihedfm, iddnfm, ext_cbc, ext_heads, ext_ddn, ext_rch, rch_user, nrchop, ext_wel, wel_user, ext_drn, drn_elev_fn, drn_cond_fn = ppMFini(MF_ws, MF_ini_fn, out = 'MF')

    if os.path.exists(rch_MM[0]):
        rch_input = rch_MM
        wel_input = wel_MM
    else:
        rch_input = rch_user
        wel_input = wel_user

    if isinstance(nper, str):
        inputFile = MMproc.readFile(MF_ws, nper.split()[0])
        l=0
        try:
            nper =  int(inputFile[l].strip())
            perlen = []
            nstp = []
            tsmult = []
            Ss_tr = []
            for l in inputFile[1:]:
                perlen.append(int(l.split()[0]))
                nstp.append(int(l.split()[1]))
                tsmult.append(int(l.split()[2]))
                Ss_tr.append(l.split()[3])
        except:
            print "Unexpected error in the time file:\n", sys.exc_info()[0]
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

    top_path = os.path.join(MF_ws, top_fn)
    top_array = np.zeros((nrow,ncol))
    top_array = MM_PROCESS.convASCIIraster2array(top_path, top_array)
    top_array = list(top_array)

    botm_array = np.zeros((nrow,ncol, nlay))
    for l in range(nlay):
        botm_path = os.path.join(MF_ws, botm_fn[l])
        botm_array[:,:,l] = MM_PROCESS.convASCIIraster2array(botm_path, botm_array[:,:,l])
    botm_array = list(botm_array)

    ibound_array = np.zeros((nrow,ncol, nlay), dtype = int)
    for l in range(nlay):
        ibound_path = os.path.join(MF_ws, ibound_fn[l])
        ibound_array[:,:,l] = MM_PROCESS.convASCIIraster2array(ibound_path, ibound_array[:,:,l])
    ibound_array = list(ibound_array)

    strt_array = np.zeros((nrow,ncol, nlay))
    for l in range(nlay):
        strt_path = os.path.join(MF_ws, strt_fn[l])
        strt_array[:,:,l] = MM_PROCESS.convASCIIraster2array(strt_path, strt_array[:,:,l])
    strt_array = list(strt_array)

    hk_array = np.zeros((nrow,ncol, nlay))
    for l in range(nlay):
        hk_path = os.path.join(MF_ws, hk_fn[l])
        hk_array[:,:,l] = MM_PROCESS.convASCIIraster2array(hk_path, hk_array[:,:,l])
    hk_array = list(hk_array)

    vka_array = np.zeros((nrow,ncol, nlay))
    for l in range(nlay):
        vka_path = os.path.join(MF_ws, vka_fn[l])
        vka_array[:,:,l] = MM_PROCESS.convASCIIraster2array(vka_path, vka_array[:,:,l])
    vka_array = list(vka_array)

    ss_array = np.zeros((nrow,ncol, nlay))
    for l in range(nlay):
        ss_path = os.path.join(MF_ws, ss_fn[l])
        ss_array[:,:,l] = MM_PROCESS.convASCIIraster2array(ss_path, ss_array[:,:,l])
    ss_array = list(ss_array)

    sy_array = np.zeros((nrow,ncol, nlay))
    for l in range(nlay):
        sy_path = os.path.join(MF_ws, sy_fn[l])
        sy_array[:,:,l] = MM_PROCESS.convASCIIraster2array(sy_path, sy_array[:,:,l])
    sy_array = list(sy_array)

# RECHARGE
    if isinstance(rch_input,float):
        rch_array = rch_input
        print '\nRecharge input: %s' % str(rch_input)
    else:
        rch_array = []
        print '\nRecharge input: %s' % rch_input[0]
        try:
            h5_rch = h5py.File(rch_input[0], 'r')
            for n in range(nper):
                rch_array.append(h5_rch[rch_input[1]][n])
            h5_rch.close
        except:
            rch_array = rch_dft
            print 'WARNING!\nNo valid recharge package file(s) provided, running MODFLOW using default recharge value: %.3G' % rch_dft
            rch_input = rch_dft

# WELL
    if isinstance(wel_input,float):
        wel_array = wel_input
        print '\nWell input: %s' % str(wel_input)
    else:
        wel_array = []
        print '\nWell input: %s' % wel_input[0]
        try:
            h5_wel = h5py.File(wel_input[0], 'r')
            for n in range(nper):
                wel_array.append(h5_wel[wel_input[1]][n])
            h5_wel.close
        except:
            wel_array = wel_dft
            print 'WARNING!\nNo valid well package file(s) provided, running MODFLOW using default well value: %.3G' % wel_dft
            wel_input = wel_dft

# DRAIN
    # in layer 1, DRN are attributed to topgraphy level
    layer_row_column_elevation_cond = [[]]
    ri=0
    for r in top_array:
        ri += 1
        ci=0
        for v in r:
            ci += 1
            layer_row_column_elevation_cond[0].append([1,ri,ci,v,1E5])
    # DRN in other layers depend on user inputESRI ASCII files
    if nlay>1 and isinstance(drn_cond_fn[0], str):
        drn_elev_array = np.zeros((nrow,ncol, nlay-1))
        for l in range(nlay-1):
            drn_elev_path = os.path.join(MF_ws, drn_elev_fn[l])
            drn_elev_array[:,:,l] = MM_PROCESS.convASCIIraster2array(drn_elev_path, drn_elev_array[:,:,l])
        drn_cond_array = np.zeros((nrow,ncol, nlay-1))
        for l in range(nlay-1):
            drn_cond_path = os.path.join(MF_ws, drn_cond_fn[l])
            drn_cond_array[:,:,l] = MM_PROCESS.convASCIIraster2array(drn_cond_path, drn_cond_array[:,:,l])
        for l in range(nlay-1):
            for i in range(nrow):
                for j in range(ncol):
                    if drn_elev_array[i,j,l]<>0:
                        if drn_elev_array[i,j,l]<0:
                            drn_elev = botm_array[i][j][l] - (botm_array[i][j][l]-botm_array[i][j][l-1])/10
                        else:
                            drn_elev = drn_elev_array[i,j,l]
                        layer_row_column_elevation_cond[0].append([l+2, i+1, j+1, drn_elev, drn_cond_array[i,j,l]])   #+np.random.normal(0,0.05,1)])

# average for 1st SS stress period
    if dum_sssp1 == 1:
        if isinstance(rch_input,tuple):
            rch_array = np.asarray(rch_array)
            rch_SS = np.zeros((nrow,ncol))
            for n in range(nper):
                rch_SS += rch_array[n,:,:]
            rch_SS = rch_SS/nper
            rch_array = list(rch_array)
            rch_array.insert(0, rch_SS)
        if isinstance(wel_input,tuple):
            wel_array = np.asarray(wel_array)
            wel_SS = np.zeros((nrow,ncol))
            for n in range(nper):
                wel_SS += wel_array[n,:,:]
            wel_SS = wel_SS/nper
            wel_array = list(wel_array)
            wel_array.insert(0, wel_SS)
        nper +=  1
        perlen.insert(0,1)
        nstp.insert(0,1)
        tsmult.insert(0,1)
        Ss_tr.insert(0, 'SS')

    layer_row_column_Q = []
    for n in range(nper):
        layer_row_column_Q.append([])
        for r in range(nrow):
            for c in range(ncol):
                if isinstance(wel_array, float):
                    layer_row_column_Q[n].append([1,r+1,c+1,wel_array])
                else:
                    layer_row_column_Q[n].append([1,r+1,c+1,wel_array[n][r][c]])

    # 2 - create the modflow packages files
    # MFfile initialization
    mfmain = mf.modflow(modelname = modelname, exe_name = exe_name, namefile_ext = namefile_ext, version = 'mf2000', model_ws= MF_ws)
    # dis initialization
    for i in range(nper):
        if Ss_tr[i] == 'TR':
            Ss_tr[i] = False
        elif Ss_tr[i] == 'SS':
            Ss_tr[i] = True
        else:
            print '\nVariable Ss_tr from the DIS package is not correct, check the MODFLOW manual'
            sys.exit()
    # dis package
    dis = mf.mfdis(model = mfmain, nrow = nrow, ncol = ncol, nlay = nlay, nper = nper, delr = delr, delc = delc, laycbd = laycbd, top = top_array, botm = botm_array, perlen = perlen, nstp = nstp, tsmult = tsmult, itmuni = itmuni, lenuni = lenuni, steady = Ss_tr, extension = ext_dis)
    del botm_array
    dis.write_file()
    # bas package
    bas = mf.mfbas(model = mfmain, ibound = ibound_array, strt = strt_array, hnoflo = hnoflo, extension = ext_bas)
    del strt_array
    bas.write_file()
    # lpf initialization
    lpf = mf.mflpf(model = mfmain, hdry = hdry, laytyp = laytyp, layavg = layavg, chani = chani, layvka = layvka, laywet = laywet, hk = hk_array, vka = vka_array, ss = ss_array, sy = sy_array, extension=ext_lpf)
    del hk_array, vka_array, ss_array, sy_array
    lpf.write_file()
    # rch initialization
    rch = mf.mfrch(mfmain, irchcb=lpf.ilpfcb, nrchop=nrchop, rech=rch_array, extension = ext_rch)
    del rch_array
    rch.write_file
    # wel initialization
    if isinstance(wel_input,float):
        if wel_input<>0:
            wel = mf.mfwel(mfmain, iwelcb = lpf.ilpfcb, layer_row_column_Q = layer_row_column_Q, extension = ext_wel)
            wel.write_file()
    else:
        wel = mf.mfwel(mfmain, iwelcb = lpf.ilpfcb, layer_row_column_Q = layer_row_column_Q, extension = ext_wel)
        wel.write_file()
    del layer_row_column_Q
    # drn package initialization
    drn = mf.mfdrn(model = mfmain, idrncb=lpf.ilpfcb, layer_row_column_elevation_cond = layer_row_column_elevation_cond, extension = ext_drn)
    del layer_row_column_elevation_cond
    drn.write_file()
    # ghb package initialization
#    ghb = mf.mfghb(model=mfmain, igbhcb = lpf.ilpfcb, layer_row_column_head_cond = [[2,13,3,90,50]])
#   ghb.write_file()
    # output control initialization
    oc = mf.mfoc(mfmain, ihedfm=ihedfm, iddnfm=iddnfm, item2=[[0,1,1,1]], item3=[[0,0,1,0]], extension=[ext_oc,ext_cbc,ext_heads,ext_ddn])
    oc.write_file()
    # select one of the 3 below (i.e. pcg or sip or sor)
    # preconditionned conjugate-gradient initialization
    pcg = mf.mfpcg(mfmain, mxiter = 150, iter1=75, hclose=1e-1, rclose=1e-1, npcond = 1, relax = 1)
    pcg.write_file()
    # sip
#    sip = mf.mfsip(mfmain, hclose=1e-3)
#    sip.write_file()
    # sor
#    sor = mf.mfsor(mfmain, hclose=1e-3)
#    sor.write_file()
    h_MF_fn = os.path.join(MF_ws, modelname + "." + ext_heads)
    if os.path.exists(h_MF_fn):
        os.remove(h_MF_fn)
    cbc_MF_fn = os.path.join(MF_ws, modelname + "." + ext_cbc)
    if os.path.exists(cbc_MF_fn):
        os.remove(cbc_MF_fn)

    # run MODFLOW and read the heads back into Python
    mfmain.write_name_file()
    mfmain.run_model(pause = False, report = report)

    # extract heads
    print''
    h5_MF_fn = os.path.join(MF_ws, '_h5_MF.h5')
    h5_MF = h5py.File(h5_MF_fn, 'w')
    try:
        h = mfrdbin.mfhdsread(mfmain, 'LF95').read_all(h_MF_fn)
    except:
        raise ValueError, '\nMODFLOW error!\nCheck the MODFLOW list file in folder:\n%s' % MF_ws
    if len(h[1])<sum(nstp):
        raise ValueError, '\nMODFLOW error!\nCheck the MODFLOW list file in folder:\n%s' % MF_ws
    print ''

    # extract cell-by-cell budget
    cbc = mfrdbin.mfcbcread(mfmain, 'LF95').read_all(cbc_MF_fn)
    h5_MF.create_dataset('cbc_nam', data = np.asarray(cbc[2]))

    print'\nStoring heads and cbc terms into HDF5 file %s' % (h5_MF_fn)

    if dum_sssp1 == 1:
        if chunks == 1:
            h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1][1:]), chunks = (1,nrow,ncol,nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf
            h5_MF.create_dataset(name = 'cbc', data = np.asarray(cbc[1][1:]), chunks = (1,len(h5_MF['cbc_nam']),nrow,ncol,nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
        else:
            h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1][1:]))
            h5_MF.create_dataset(name = 'cbc', data = np.asarray(cbc[1][1:]))
        del h, cbc
        nper = nper - 1
        perlen = perlen[1:]
        nstp = nstp[1:]
    else:
        # Not tested
        if chunks == 1:
            h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1]), chunks = (1,nrow,ncol,nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
            h5_MF.create_dataset(name = 'cbc', data = np.asarray(cbc[1]), chunks = (1,len(h5_MF['cbc_nam']),nrow,ncol,nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
        else:
            h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1]))
            h5_MF.create_dataset(name = 'cbc', data = np.asarray(cbc[1]))
        del h, cbc
    h5_MF.close()
    del mfmain, dis, bas, lpf, rch, wel, drn, oc, pcg
    # to delete MF bin files and save disk space
    os.remove(h_MF_fn)
    os.remove(cbc_MF_fn)

    return np.asarray(top_array), h5_MF_fn
    del nrow, ncol, delr, delc, nlay, perlen, nper, hnoflo, hdry, ibound_array, laytyp, h5_MF_fn, top_array, inputFileMF_fn, lenuni