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

import sys
import os
import numpy as np
import mf
import mfreadbinaries as mfrdbin
import h5py
import MARMITESprocess as MMproc

#####################################
# TODO modify the reading of the MF ini file
def ppMFini(MF_ws, MF_ini_fn, out = 'MF', numDays = -1):

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
        version =  str(inputFile[l].strip())
        l += 1
        dum_sssp1 = int(inputFile[l].strip())
        l += 1
        # dis
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
                nper = numDays
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
        top = [inputFile[l].strip()]
        l += 1
        botm = []
        for i in range(nlay):
            botm.append(inputFile[l].split()[i])
        l += 1
        perlen = []
        nstp = []
        tsmult = []
        Ss_tr = []
        if timedef > 0:
            perlen_tmp =  inputFile[l].split()
            l += 1
            nstp_tmp =  inputFile[l].split()
            l += 1
            tsmult_tmp =  inputFile[l].split()
            l += 1
            Ss_tr_tmp =  inputFile[l].split()
            for i in range(nper):
                perlen.append(int(perlen_tmp[i]))
                nstp.append(int(nstp_tmp[i]))
                tsmult.append(int(tsmult_tmp[i]))
                Ss_tr.append(str(Ss_tr_tmp[i]))
                if nstp[i]>perlen[i]:
                    print "ERROR!\nMM doesn't accept nstp < perlen!"
                    sys.exit()
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
        # bas
        ext_bas = str(inputFile[l].strip())
        l += 1
        ibound = []
        for i in range(nlay):
            ibound.append(inputFile[l].split()[i])
        l += 1
        strt = []
        for i in range(nlay):
            strt.append(inputFile[l].split()[i])
        l += 1
        hnoflo = float(inputFile[l].strip())
        # lpf
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
        for i in range(nlay):
            layavg.append(int(inputFile[l].split()[i]))
        l += 1
        chani = []
        for i in range(nlay):
            chani.append(int(inputFile[l].split()[i]))
        l += 1
        layvka = []
        for i in range(nlay):
            layvka.append(int(inputFile[l].split()[i]))
        l += 1
        laywet = []
        for i in range(nlay):
            laywet.append(int(inputFile[l].split()[i]))
        l += 1
        hk = []
        for i in range(nlay):
            hk.append(inputFile[l].split()[i])
        l += 1
        vka = []
        for i in range(nlay):
            vka.append(inputFile[l].split()[i])
        l += 1
        ss = []
        for i in range(nlay):
            ss.append(inputFile[l].split()[i])
        l += 1
        sy = []
        for i in range(nlay):
            sy.append(inputFile[l].split()[i])
        # uzf
        l += 1
        ext_uzf = str(inputFile[l].strip())
        l += 1
        nuztop = int(inputFile[l].strip())
        l += 1
        iuzfopt = int(inputFile[l].strip())
        l += 1
        irunflg = int(inputFile[l].strip())
        l += 1
        ietflg = int(inputFile[l].strip())
        l += 1
        iuzfcb1 = int(inputFile[l].strip())
        l += 1
        iuzfcb2 = int(inputFile[l].strip())
        l += 1
        ntrail2 = int(inputFile[l].strip())
        l += 1
        nsets = int(inputFile[l].strip())
        l += 1
        nuzgag = int(inputFile[l].strip())
        l += 1
        surfdep = float(inputFile[l].strip())
        l += 1
        eps = float(inputFile[l].strip())
        l += 1
        thts = float(inputFile[l].strip())
        l += 1
        thti = float(inputFile[l].strip())
        row_col_iftunit_iuzfopt = []
        if nuzgag > 0:
            l += 1
            iuzrow =  inputFile[l].split()
            l += 1
            iuzcol = inputFile[l].split()
            l += 1
            iftunit = inputFile[l].split()
            l += 1
            iuzopt = inputFile[l].split()
            l += 1
            for g in range(nuzgag):
                row_col_iftunit_iuzfopt.append([[int(iuzrow[g]),int(iuzcol[g]),int(iftunit[g]),int(iuzopt[g])]])
        finf_user = float(inputFile[l].strip())
        # wel
        l += 1
        ext_wel = str(inputFile[l].strip())
        l += 1
        wel_user = float(inputFile[l].strip())
        l += 1
        # drn
        ext_drn = str(inputFile[l].strip())
        l += 1
        drn_elev = []
        for i in range(nlay):
            drn_elev.append(inputFile[l].split()[i])
        l += 1
        drn_cond = []
        for i in range(nlay):
            drn_cond.append(inputFile[l].split()[i])
        # oc
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
        # pcg
        hclose = float(inputFile[l].strip())
        l += 1
        rclose = float(inputFile[l].strip())
        l += 1
    except:
        print "Unexpected error in the input file:\n", sys.exc_info()[0]
        sys.exit()
    del inputFile

    if out == 'MF':
        return modelname, namefile_ext, exe_name, version, dum_sssp1, ext_dis, nlay, ncol, nrow, nper, itmuni, lenuni,laycbd, delr, delc, top, botm, perlen, nstp, tsmult, Ss_tr, ext_bas, ibound, strt, hnoflo,ext_lpf, ilpfcb, hdry, nplpf, laytyp, layavg, chani, layvka, laywet, hk, vka, ss, sy,ext_oc, ihedfm, iddnfm, ext_cbc, ext_heads, ext_ddn, hclose, rclose, ext_wel, wel_user, ext_drn, drn_elev, drn_cond, ext_uzf, nuztop, iuzfopt, irunflg, ietflg, iuzfcb1, iuzfcb2, ntrail2, nsets, nuzgag, surfdep, eps, thts, thti, row_col_iftunit_iuzfopt, finf_user
    elif out == 'MM':
        return nrow, ncol, delr, delc, top, reggrid, nlay, nper, perlen, nstp, timedef, hnoflo, hdry, laytyp, lenuni, itmuni, ibound

#####################################

def ppMFtime(MM_ws, MF_ws, outpathname, nper, perlen, nstp, timedef, inputDate_fn, inputZON_TS_RF_fn, inputZON_TS_PET_fn, inputZON_TS_RFe_fn, inputZON_TS_PE_fn, inputZON_TS_E0_fn, NMETEO, NVEG, NSOIL):

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
    del inputFile, inputFileRF, inputFilePET, inputFilePE, inputFileRFe, inputFileE0

    if timedef == 0:
        try:
            if isinstance(nper, str):
                perlenmax = int(nper.split()[1].strip())
        except:
            print '\nError in nper format of the MODFLOW ini file!\n'
            sys.exit()
        RF_stp = []
        PET_stp = []
        RFe_stp = []
        PE_stp = []
        E0_stp = []
        RF_stp_tmp = []
        PET_stp_tmp = []
        PE_stp_tmp=[]
        E0_stp_tmp = []
        nper = 0
        perlen = []
        perlen_tmp = 0
        c = 0
        for n in range(NMETEO):
            RF_stp.append([])
            RF_stp_tmp.append(0.0)
            PET_stp.append([])
            RFe_stp.append([])
            PET_stp_tmp.append([])
            for v in range(NVEG):
                PET_stp[n].append([])
                RFe_stp[n].append([])
                PET_stp_tmp[n].append(0.0)
            PE_stp.append([])
            PE_stp_tmp.append([])
            for v in range(NSOIL):
                PE_stp[n].append([])
                PE_stp_tmp[n].append(0.0)
            E0_stp.append([])
            E0_stp_tmp.append(0.0)
        if perlenmax < 2:
            print '\nperlenmax must be higher than 1!\nCorrect perlenmax in the MODFLOW ini file or select the daily option.'
            sys.exit()
        for j in range(len(d)):
                if RFe_d[:,:,j].sum()>0.0:
                    if c == 1:
                        for n in range(NMETEO):
                            RF_stp[n].append(RF_stp_tmp[n]/perlen_tmp)
                            RF_stp_tmp[n] = 0.0
                            for v in range(NVEG):
                                PET_stp[n][v].append(PET_stp_tmp[n][v]/perlen_tmp)
                                RFe_stp[n][v].append(0.0)
                                PET_stp_tmp[n][v] = 0.0
                            for s in range(NSOIL):
                                PE_stp[n][s].append(PE_stp_tmp[n][s]/perlen_tmp)
                                PE_stp_tmp[n][s] = 0.0
                            E0_stp[n].append(E0_stp_tmp[n]/perlen_tmp)
                            E0_stp_tmp[n] = 0.0
                        perlen.append(perlen_tmp)
                        perlen_tmp = 0
                    for n in range(NMETEO):
                        RF_stp[n].append(RF_d[n,j])
                        for v in range(NVEG):
                            PET_stp[n][v].append(PET_d[n,v,j])
                            RFe_stp[n][v].append(RFe_d[n,v,j])
                        for s in range(NSOIL):
                            PE_stp[n][s].append(PE_d[n,s,j])
                        E0_stp[n].append(E0_d[n,j])
                    nper += 1
                    perlen.append(1)
                    c = 0
                else:
                    if perlen_tmp < perlenmax:
                        for n in range(NMETEO):
                            RF_stp_tmp[n]  += RF_d[n,j]
                            for v in range(NVEG):
                                PET_stp_tmp[n][v] += PET_d[n,v,j]
                            for s in range(NSOIL):
                                PE_stp_tmp[n][s]  += PE_d[n,s,j]
                            E0_stp_tmp[n]  += E0_d[n,j]
                        if c == 0:
                            nper += 1
                        perlen_tmp += 1
                        c = 1
                    else:
                        for n in range(NMETEO):
                            RF_stp[n].append(RF_stp_tmp[n]/perlen_tmp)
                            RF_stp_tmp[n] = 0.0
                            for v in range(NVEG):
                                PET_stp[n][v].append(PET_stp_tmp[n][v]/perlen_tmp)
                                RFe_stp[n][v].append(0.0)
                                PET_stp_tmp[n][v] = 0.0
                            for s in range(NSOIL):
                                PE_stp[n][s].append(PE_stp_tmp[n][s]/perlen_tmp)
                                PE_stp_tmp[n][s] = 0.0
                            E0_stp[n].append(E0_stp_tmp[n]/perlen_tmp)
                            E0_stp_tmp[n] = 0.0
                        perlen.append(perlen_tmp)
                        for n in range(NMETEO):
                            RF_stp_tmp[n]  += RF_d[n,j]
                            for v in range(NVEG):
                                PET_stp_tmp[n][v] += PET_d[n,v,j]
                            for s in range(NSOIL):
                                PE_stp_tmp[n][s]  += PE_d[n,s,j]
                            E0_stp_tmp[n]  += E0_d[n,j]
                        nper += 1
                        perlen_tmp = 1
                        c = 1
        if c == 1:
            for n in range(NMETEO):
                RF_stp[n].append(RF_stp_tmp[n]/perlen_tmp)
                RF_stp_tmp[n] = 0.0
                for v in range(NVEG):
                    PET_stp[n][v].append(PET_stp_tmp[n][v]/perlen_tmp)
                    RFe_stp[n][v].append(0.0)
                    PET_stp_tmp[n][v] = 0.0
                for s in range(NSOIL):
                    PE_stp[n][s].append(PE_stp_tmp[n][s]/perlen_tmp)
                    PE_stp_tmp[n][s] = 0.0
                E0_stp[n].append(E0_stp_tmp[n]/perlen_tmp)
                E0_stp_tmp[n] = 0.0
            perlen.append(perlen_tmp)
        del perlen_tmp
        del RF_d, RFe_d, PET_d, PE_d, E0_d, RF_stp_tmp, PET_stp_tmp, PE_stp_tmp, E0_stp_tmp, c
        perlen = np.asarray(perlen, dtype = np.int)
        nstp = np.ones(nper, dtype = np.int)
        tsmult = nstp
        Ss_tr = []
        for n in range(nper):
            Ss_tr.append('TR')
    elif timedef>0:
        RF_stp = np.zeros([NMETEO,sum(nstp)])
        PET_stp = np.zeros([NMETEO,NVEG,sum(nstp)])
        RFe_stp = np.zeros([NMETEO,NVEG,sum(nstp)])
        PE_stp = np.zeros([NMETEO,NSOIL,sum(nstp)])
        E0_stp = np.zeros([NMETEO,sum(nstp)])
        stp = 0
        for per in range(nper):
            for stp in range(nstp[per]):
                tstart = sum(perlen[0:per])+stp*(perlen[per]/nstp[per])
                tend = sum(perlen[0:per])+(perlen[per]/nstp[per])*(1+stp)
                for n in range(NMETEO):
                    RF_stp[n,sum(nstp[0:per])+stp] += RF_d[n,tstart:tend].sum()/(perlen[per]/nstp[per])
                    for v in range(NVEG):
                        PET_stp[n,v,sum(nstp[0:per])+stp] += PET_d[n,v,tstart:tend].sum()/(perlen[per]/nstp[per])
                        RFe_stp[n,v,sum(nstp[0:per])+stp] += RFe_d[n,v,tstart:tend].sum()/(perlen[per]/nstp[per])
                    for s in range(NSOIL):
                        PE_stp[n,s,sum(nstp[0:per])+stp] += PE_d[n,s,tstart:tend].sum()/(perlen[per]/nstp[per])
                    E0_stp[n,sum(nstp[0:per])+stp] += E0_d[n,tstart:tend].sum()/(perlen[per]/nstp[per])
        tsmult = np.ones(nper, dtype = np.int)
        Ss_tr = []
        for n in range(nper):
            Ss_tr.append('TR')

    outpathname = os.path.join(MF_ws, outpathname)
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

def ppMF(MM_ws, xllcorner, yllcorner, MF_ws, MF_ini_fn, finf_MM = "", finf_user = None, wel_MM = "", wel_user = None, report = None, verbose = 1, chunks = 0, MFtime_fn = None, numDays = -1):

    ##################

    def checkarray(nrow, ncol, nlay, MF_ws, var):
        try:
            if len(var)>1:
                lst_out = []
                for v in var:
                    lst_out.append(float(v))
            else:
                lst_out = float(var[0])
        except:
            array = np.zeros((nrow,ncol,len(var)))
            l = 0
            for v in var:
                if isinstance(v, str):
                    array_path = os.path.join(MF_ws, v)
                    array[:,:,l] = MM_PROCESS.convASCIIraster2array(array_path, array[:,:,l])
                else:
                    print'\nERROR!\nMODFLOW ini file incorrect, check files or values %s' % var
                l += 1
            if len(var)>1:
                lst_out = list(array)
            else:
                lst_out = list(array[:,:,0])

        return lst_out

    ##################

    if verbose == 0:
        print '--------------'

    modelname, namefile_ext, exe_name, version, dum_sssp1, ext_dis, nlay, ncol, nrow, nper, itmuni, lenuni,laycbd, delr, delc, top, botm, perlen, nstp, tsmult, Ss_tr, ext_bas, ibound, strt, hnoflo,ext_lpf, ilpfcb, hdry, nplpf, laytyp, layavg, chani, layvka, laywet, hk, vka, ss, sy,ext_oc, ihedfm, iddnfm, ext_cbc, ext_heads, ext_ddn, hclose, rclose, ext_wel, wel_user, ext_drn, drn_elev, drn_cond, ext_uzf, nuztop, iuzfopt, irunflg, ietflg, iuzfcb1, iuzfcb2, ntrail2, nsets, nuzgag, surfdep, eps, thts, thti, row_col_iftunit_iuzfopt, finf_user = ppMFini(MF_ws, MF_ini_fn, out = 'MF', numDays = numDays)

    if os.path.exists(finf_MM[0]):
        finf_input = finf_MM
        wel_input = wel_MM
    else:
        finf_input = finf_user
        wel_input = wel_user

    if MFtime_fn <> None:
        inputFile = MMproc.readFile(MF_ws, MFtime_fn)
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

    top = checkarray(nrow, ncol, nlay, MF_ws, top)
    botm = checkarray(nrow, ncol, nlay, MF_ws, botm)
    ibound = checkarray(nrow, ncol, nlay, MF_ws, ibound)
    iuzfbnd = list(np.asarray(np.asarray(ibound)[:,:,0] > 0.0, dtype = np.int))
    strt = checkarray(nrow, ncol, nlay, MF_ws, strt)
    hk = checkarray(nrow, ncol, nlay, MF_ws, hk)
    vka = checkarray(nrow, ncol, nlay, MF_ws, vka)
    ss = checkarray(nrow, ncol, nlay, MF_ws, ss)
    sy = checkarray(nrow, ncol, nlay, MF_ws, sy)

# FINF
    if isinstance(finf_input,float):
        finf_array = finf_input
        print '\nInfiltration (UZF1 package) input: %s' % str(finf_input)
    else:
        finf_array = []
        print '\nInfiltration (UZF1 package) input: %s' % finf_input[0]
        try:
            h5_finf = h5py.File(finf_input[0], 'r')
            for n in range(nper):
                finf_array.append(h5_finf[finf_input[1]][n])
            h5_finf.close()
        except:
            finf_array = finf_user
            print 'WARNING!\nNo valid UZF1 package file(s) provided, running MODFLOW using user-input UZF1infiltration value: %.3G' % finf_user
            finf_input = finf_user

# WELL
    if isinstance(wel_input,float):
        wel_array = wel_input
        print '\nWell (WEL package) input: %s' % str(wel_input)
    else:
        wel_array = []
        print '\nWell (WEL package) input: %s' % wel_input[0]
        try:
            h5_wel = h5py.File(wel_input[0], 'r')
            for n in range(nper):
                wel_array.append(h5_wel[wel_input[1]][n])
            h5_wel.close()
        except:
            wel_array = wel_user
            print 'WARNING!\nNo valid WEL package file(s) provided, running MODFLOW using user-input well value: %.3G' % wel_user
            wel_input = wel_user

# DRAIN
    l = 0
    layer_row_column_elevation_cond = [[]]
    for d in drn_cond:
        drn_check = 1
        drn_elev_array = np.zeros((nrow,ncol))
        drn_cond_array = np.zeros((nrow,ncol))
        if isinstance(d, str):
            drn_elev_path = os.path.join(MF_ws, drn_elev[l])
            drn_elev_array[:,:] = MM_PROCESS.convASCIIraster2array(drn_elev_path, drn_elev_array[:,:])
            drn_cond_path = os.path.join(MF_ws, drn_cond[l])
            drn_cond_array[:,:] = MM_PROCESS.convASCIIraster2array(drn_cond_path, drn_cond_array[:,:])
        else:
            drn_elev_array[:,:] = drn_elev[l]
            drn_cond_array[:,:] = drn_cond[l]
        for i in range(nrow):
            for j in range(ncol):
                if drn_elev_array[i,j]<>0:
                    if drn_elev_array[i,j]<0:
                        if isinstance(botm[l], float):
                            drn_elev_tmp = botm[l]
                        else:
                            drn_elev_tmp = botm[i][j][l] #- (botm_array[i][j][l]-botm_array[i][j][l-1])/10
                    else:
                        drn_elev_tmp = drn_elev_array[i,j]
                    layer_row_column_elevation_cond[0].append([l+1, i+1, j+1, drn_elev_tmp, drn_cond_array[i,j]])
        l += 1

# average for 1st SS stress period
    if dum_sssp1 == 1:
        if isinstance(finf_input,tuple):
            finf_array = np.asarray(finf_array)
            finf_SS = np.zeros((nrow,ncol))
            for n in range(nper):
                finf_SS += finf_array[n,:,:]
            finf_SS = finf_SS/nper
            finf_array = list(finf_array)
            finf_array.insert(0, finf_SS)
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

##    if isinstance(wel_input,float) and wel_input == 0.0:
##        layer_row_column_Q = None
##    else:
    layer_row_column_Q = []
    for n in range(nper):
        layer_row_column_Q.append([])
        for r in range(nrow):
            for c in range(ncol):
                if isinstance(wel_array, float):
                    layer_row_column_Q[n].append([1,r+1,c+1,-wel_array*delr[c]*delc[r]])
                else:
                    layer_row_column_Q[n].append([1,r+1,c+1,-(wel_array[n][r][c])*delr[c]*delc[r]])

    # 2 - create the modflow packages files
    # MFfile initialization
    mfmain = mf.modflow(modelname = modelname, exe_name = exe_name, namefile_ext = namefile_ext, version = version, model_ws= MF_ws)
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
    dis = mf.mfdis(model = mfmain, nrow = nrow, ncol = ncol, nlay = nlay, nper = nper, delr = delr, delc = delc, laycbd = laycbd, top = top, botm = botm, perlen = perlen, nstp = nstp, tsmult = tsmult, itmuni = itmuni, lenuni = lenuni, steady = Ss_tr, extension = ext_dis)
    del botm
    dis.write_file()
    # bas package
    bas = mf.mfbas(model = mfmain, ibound = ibound, strt = strt, hnoflo = hnoflo, extension = ext_bas)
    del strt
    bas.write_file()
    # lpf initialization
    lpf = mf.mflpf(model = mfmain, hdry = hdry, laytyp = laytyp, layavg = layavg, chani = chani, layvka = layvka, laywet = laywet, hk = hk, vka = vka, ss = ss, sy = sy, extension=ext_lpf)
    del hk, vka, ss, sy
    lpf.write_file()
    # wel initialization
    if layer_row_column_Q <> None:
        wel = mf.mfwel(mfmain, iwelcb = lpf.ilpfcb, layer_row_column_Q = layer_row_column_Q, extension = ext_wel)
        wel.write_file()
    del layer_row_column_Q
    # drn package initialization
    if drn_check == 1:
        drn = mf.mfdrn(model = mfmain, idrncb=lpf.ilpfcb, layer_row_column_elevation_cond = layer_row_column_elevation_cond, extension = ext_drn)
        del layer_row_column_elevation_cond
        drn.write_file()
    # uzf initialization
    uzf = mf.mfuzf1(mfmain, nuztop = nuztop, iuzfopt = iuzfopt, irunflg = irunflg, ietflg = ietflg, iuzfcb1 = iuzfcb1, iuzfcb2 = iuzfcb2, ntrail2 = ntrail2, nsets = nsets, nuzgag = nuzgag, surfdep = surfdep, iuzfbnd = iuzfbnd, eps = eps, thts = thts, thti = thti, row_col_iftunit_iuzfopt = row_col_iftunit_iuzfopt, finf = finf_array, extension = ext_uzf)
    del finf_array
    uzf.write_file()
    # output control initialization
    oc = mf.mfoc(mfmain, ihedfm=ihedfm, iddnfm=iddnfm, item2=[[0,1,1,1]], item3=[[0,0,1,0]], extension=[ext_oc,ext_heads,ext_ddn,ext_cbc])
    oc.write_file()
    # select one of the 3 below (i.e. pcg or sip or sor)
    # preconditionned conjugate-gradient initialization
    pcg = mf.mfpcg(mfmain, mxiter = 150, iter1=75, hclose=hclose, rclose=rclose, npcond = 1, relax = 1)
    pcg.write_file()
    # sip
#    sip = mf.mfsip(mfmain, hclose=1e-3)
#    sip.write_file()
    # sor
#    sor = mf.mfsor(mfmain, hclose=1e-3)
#    sor.write_file()

    h_MF_fn = os.path.join(MF_ws, modelname + "." + ext_heads)
    cbc_MF_fn = os.path.join(MF_ws, modelname + "." + ext_cbc)
    cbc_MFuzf_fn = os.path.join(MF_ws, modelname + ".uzfbt1")

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
        h5_MF.close()
        raise ValueError, '\nMODFLOW error!\nCheck the MODFLOW list file in folder:\n%s' % MF_ws
    if len(h[1])<sum(nstp):
        h5_MF.close()
        raise ValueError, '\nMODFLOW error!\nCheck the MODFLOW list file in folder:\n%s' % MF_ws
    print ''

    # extract cell-by-cell budget
    cbc = mfrdbin.mfcbcread(mfmain, 'LF95').read_all(cbc_MF_fn)
    h5_MF.create_dataset('cbc_nam', data = np.asarray(cbc[2][1]))
    print ''

    # extract cell-by-cell uzf budget
    cbc_uzf = mfrdbin.mfcbcread(mfmain, 'LF95').read_all(cbc_MFuzf_fn)
    h5_MF.create_dataset('cbc_uzf_nam', data = np.asarray(cbc_uzf[2][1]))

    print'\nStoring heads and cbc terms into HDF5 file\n%s' % (h5_MF_fn)

    if dum_sssp1 == 1:
        if chunks == 1:
            h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1][1:]), chunks = (1,nrow,ncol,nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf
            h5_MF.create_dataset(name = 'cbc', data = np.asarray(cbc[1][1:]), chunks = (1,len(h5_MF['cbc_nam']),nrow,ncol,nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
            h5_MF.create_dataset(name = 'cbc_uzf', data = np.asarray(cbc_uzf[1][1:]), chunks = (1,len(h5_MF['cbc_uz_nam']),nrow,ncol,nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
        else:
            h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1][1:]))
            h5_MF.create_dataset(name = 'cbc', data = np.asarray(cbc[1][1:]))
            h5_MF.create_dataset(name = 'cbc_uzf', data = np.asarray(cbc_uzf[1][1:]))
        del h, cbc
        nper = nper - 1
        perlen = perlen[1:]
        nstp = nstp[1:]
    else:
        if chunks == 1:
            h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1]), chunks = (1,nrow,ncol,nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
            h5_MF.create_dataset(name = 'cbc', data = np.asarray(cbc[1]), chunks = (1,len(h5_MF['cbc_nam']),nrow,ncol,nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
            h5_MF.create_dataset(name = 'cbc_uzf', data = np.asarray(cbc_uzf[1]), chunks = (1,len(h5_MF['cbc_uzf_nam']),nrow,ncol,nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
        else:
            h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1]))
            h5_MF.create_dataset(name = 'cbc', data = np.asarray(cbc[1]))
            h5_MF.create_dataset(name = 'cbc_uzf', data = np.asarray(cbc_uzf[1]))
        del h, cbc, cbc_uzf
    h5_MF.close()
    # to delete MF bin files and save disk space
#    os.remove(h_MF_fn)
#    os.remove(cbc_MF_fn)
#    os.remove(cbc_MFuzf_fn)

    return h5_MF_fn
    del nrow, ncol, delr, delc, nlay, perlen, nper, hnoflo, hdry, ibound_array, laytyp, h5_MF_fn, top_array, inputFileMF_fn, lenuni