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

import sys, os, traceback
import numpy as np
import mf
import mfreadbinaries as mfrdbin
import h5py
import MARMITESprocess_v3 as MMproc

#####################################
class MF():
    def __init__(self, MM_ws, MF_ws, MF_ini_fn, xllcorner, yllcorner, numDays = -1):

        ''' read input file (called _input.ini in the MARMITES workspace
        the first character on the first line has to be the character used to comment
        the file can contain any comments as the user wish, but the sequence of the input has to be respected'''

        inputFile = MMproc.readFile(MF_ws, MF_ini_fn)
        l=0
        try:
            self.h5_MF_fn = os.path.join(MF_ws, '_h5_MF.h5')
            self.xllcorner = xllcorner
            self.yllcorner = yllcorner
            self.modelname =  str(inputFile[l].strip())
            l += 1
            self.namefile_ext = str(inputFile[l].strip())
            l += 1
            self.exe_name = str(inputFile[l].strip())
            l += 1
            self.version =  str(inputFile[l].strip())
            l += 1
            self.dum_sssp1 = int(inputFile[l].strip())
            l += 1
            # COMPULSORY PACKAGES
            # dis
            self.ext_dis = str(inputFile[l].strip())
            l += 1
            self.nlay = int(inputFile[l].strip())
            l += 1
            self.ncol = int(inputFile[l].strip())
            l += 1
            self.nrow = int(inputFile[l].strip())
            l += 1
            nper = inputFile[l].strip()
            try:
                nper = int(nper)
            except:
                pass
            if isinstance(nper, int):
                if nper <0:
                # daily nper, perlen, nstp, tsmult
                    self.nper = numDays
                    self.timedef = -1
                else:
                    # read nper, perlen, nstp, tsmult from user-input ini file
                    self.nper = nper
                    self.timedef = 1
            elif isinstance(nper, str):
                # compute nper, perlen, nstp, tsmult based on daily RF analysis
                nper = nper[0]
                self.timedef = 0
                self.nper = inputFile[l].strip()
            l += 1
            self.perlen = []
            self.nstp = []
            self.tsmult = []
            self.Ss_tr = []
            if self.timedef > 0:
                #read user-input time discretization
                perlen_tmp =  inputFile[l].split()
                l += 1
                nstp_tmp =  inputFile[l].split()
                l += 1
                tsmult_tmp =  inputFile[l].split()
                l += 1
                Ss_tr_tmp =  inputFile[l].split()
                for i in range(self.nper):
                    self.perlen.append(int(perlen_tmp[i]))
                    self.nstp.append(int(nstp_tmp[i]))
                    self.tsmult.append(int(tsmult_tmp[i]))
                    self.Ss_tr.append(str(Ss_tr_tmp[i]))
                    if self.nstp[i]>self.perlen[i]:
                        raise SystemExit("FATAL ERROR!\nMM doesn't accept nstp < perlen!")
                    if self.Ss_tr[i] == 'TR':
                        self.Ss_tr[i] = False
                    elif self.Ss_tr[i] == 'SS':
                        self.Ss_tr[i] = True
                    else:
                        raise SystemExit('\nVariable Ss_tr from the DIS package is not correct, check the MODFLOW manual')
                l += 1
            elif self.timedef < 0:
                # daily data
                for t in range(self.nper):
                    self.perlen.append(1)
                    self.nstp.append(1)
                    self.tsmult.append(1)
                    self.Ss_tr.append(False)
                l += 4
            elif self.timedef == 0:
                # will be defined later in def MF.ppMFtime
                l += 4
            self.itmuni = int(inputFile[l].strip())
            l += 1
            self.lenuni = int(inputFile[l].strip())
            l += 1
            self.laycbd = []
            laycbd_tmp =  inputFile[l].split()
            for i in range(self.nlay):
                self.laycbd.append(int(laycbd_tmp[i]))
            l += 1
            self.delr = []
            delr_tmp =  inputFile[l].split()
            if len(delr_tmp)>1:
                for i in range(self.ncol):
                    self.delr.append(int(delr_tmp[i]))
            else:
                for i in range(self.ncol):
                    self.delr.append(int(delr_tmp[0]))
            l += 1
            self.delc = []
            delc_tmp =  inputFile[l].split()
            if len(delc_tmp)>1:
                for i in range(self.nrow):
                    self.delc.append(int(delc_tmp[i]))
            else:
                for i in range(self.nrow):
                    self.delc.append(int(delc_tmp[0]))
            l += 1
            self.reggrid = int(inputFile[l].strip())
            l += 1
            self.top = [inputFile[l].strip()]
            l += 1
            self.botm = []
            for i in range(self.nlay):
                self.botm.append(inputFile[l].split()[i])
            l += 1
            # bas
            self.ext_bas = str(inputFile[l].strip())
            l += 1
            self.ibound = []
            for i in range(self.nlay):
                self.ibound.append(inputFile[l].split()[i])
            l += 1
            self.strt = []
            for i in range(self.nlay):
                self.strt.append(inputFile[l].split()[i])
            l += 1
            self.hnoflo = float(inputFile[l].strip())
            # layer packages (lpf or upw)
            if self.version != 'mfnwt':
                # lpf
                l += 1
                self.ext_lpf = str(inputFile[l].strip())
                l += 1
                self.ilpfcb = int(inputFile[l].strip())
                l += 1
                self.hdry = eval(inputFile[l].strip())
                l += 1
                self.nplpf = int(inputFile[l].strip())
                l += 1
                self.laytyp = []
                laytyp_tmp =  inputFile[l].split()
                for i in range(self.nlay):
                    self.laytyp.append(int(laytyp_tmp[i]))
                l += 1
                self.layavg = []
                for i in range(self.nlay):
                    self.layavg.append(int(inputFile[l].split()[i]))
                l += 1
                self.chani = []
                for i in range(self.nlay):
                    self.chani.append(int(inputFile[l].split()[i]))
                l += 1
                self.layvka = []
                for i in range(self.nlay):
                    self.layvka.append(int(inputFile[l].split()[i]))
                l += 1
                self.laywet = []
                for i in range(self.nlay):
                    self.laywet.append(int(inputFile[l].split()[i]))
                l += 1
                self.hk = []
                for i in range(self.nlay):
                    self.hk.append(inputFile[l].split()[i])
                l += 1
                self.vka = []
                for i in range(self.nlay):
                    self.vka.append(inputFile[l].split()[i])
                l += 1
                self.ss = []
                for i in range(self.nlay):
                    self.ss.append(inputFile[l].split()[i])
                l += 1
                self.sy = []
                for i in range(self.nlay):
                    self.sy.append(inputFile[l].split()[i])
                l += 1
                if int(inputFile[l].strip()) == 1:
                    self.storagecoefficient = True
                else:
                    self.storagecoefficient = False
                l += 1
                if int(inputFile[l].strip()) == 1:
                    self.constantcv = True
                    self.nocvcorrection = True
                else:
                    self.constantcv = False
                l += 1
                if int(inputFile[l].strip()) == 1:
                    self.thickstrt = True
                else:
                    self.thickstrt = False
                l += 1
                if int(inputFile[l].strip()) == 1:
                    self.nocvcorrection = True
                else:
                    if self.constantcv != True:
                        self.nocvcorrection = False
                l += 1
                if int(inputFile[l].strip()) == 1:
                    self.novfc = True
                else:
                    self.novfc = False
                l += 14
            elif self.version == 'mfnwt':
                l += 18
                # upw
                l += 1
                self.ext_upw = str(inputFile[l].strip())
                l += 1
                self.iupwcb = int(inputFile[l].strip())
                l += 1
                self.hdry = eval(inputFile[l].strip())
                l += 1
                self.npupw = int(inputFile[l].strip())
                l += 1
                self.iphdry = int(inputFile[l].strip())
                l += 1
                self.laytyp = []
                laytyp_tmp =  inputFile[l].split()
                for i in range(self.nlay):
                    self.laytyp.append(int(laytyp_tmp[i]))
                l += 1
                self.layavg = []
                for i in range(self.nlay):
                    self.layavg.append(int(inputFile[l].split()[i]))
                l += 1
                self.chani = []
                for i in range(self.nlay):
                    self.chani.append(int(inputFile[l].split()[i]))
                l += 1
                self.layvka = []
                for i in range(self.nlay):
                    self.layvka.append(int(inputFile[l].split()[i]))
                l += 1
                self.laywet = []
                for i in range(self.nlay):
                    self.laywet.append(int(inputFile[l].split()[i]))
                l += 1
                self.hk = []
                for i in range(self.nlay):
                    self.hk.append(inputFile[l].split()[i])
                l += 1
                self.vka = []
                for i in range(self.nlay):
                    self.vka.append(inputFile[l].split()[i])
                l += 1
                self.ss = []
                for i in range(self.nlay):
                    self.ss.append(inputFile[l].split()[i])
                l += 1
                self.sy = []
                for i in range(self.nlay):
                    self.sy.append(inputFile[l].split()[i])
            else:
                print 'FATAL ERROR!\nMODFLOW version should be mf2k, mf2005 or mfnwt!'
                print 'Value %s provided in the MF ini file.' % self.versionsys.exit()
            # oc
            l += 1
            self.ext_oc = str(inputFile[l].strip())
            l += 1
            self.ihedfm = int(inputFile[l].strip())
            l += 1
            self.iddnfm = int(inputFile[l].strip())
            l += 1
            self.ext_cbc = str(inputFile[l].strip())
            l += 1
            self.ext_heads = str(inputFile[l].strip())
            l += 1
            self.ext_ddn = str(inputFile[l].strip())
            # solver
            if self.version != 'mfnwt':
                # pcg
                l += 1
                self.ext_pcg = str(inputFile[l].strip())
                l += 1
                self.hclose = float(inputFile[l].strip())
                l += 1
                self.rclose = float(inputFile[l].strip())
                l += 9
            elif self.version == 'mfnwt':
                l += 4
                self.ext_nwt = str(inputFile[l].strip())
                l += 1
                self.headtol = float(inputFile[l].strip())
                l += 1
                self.fluxtol = float(inputFile[l].strip())
                l += 1
                self.maxiterout = int(inputFile[l].strip())
                l += 1
                self.thickfact = float(inputFile[l].strip())
                l += 1
                self.linmeth = int(inputFile[l].strip())
                l += 1
                self.iprnwt = int(inputFile[l].strip())
                l += 1
                self.ibotav = int(inputFile[l].strip())
                l += 1
                self.options = str(inputFile[l].strip())
            else:
                print 'FATAL ERROR!\nMODFLOW version should be mf2k, mf2005 or mfnwt!'
                print 'Value %s provided in the MF ini file.' % self.versionsys.exit()
            # OPTIONNAL PACKAGES
            # uzf
            l += 1
            self.uzf_yn = int(inputFile[l].strip())
            if self.uzf_yn == 1 :
                l += 1
                self.ext_uzf = str(inputFile[l].strip())
                l += 1
                self.nuztop = int(inputFile[l].strip())
                l += 1
                self.iuzfopt = int(inputFile[l].strip())
                l += 1
                self.irunflg = int(inputFile[l].strip())
                l += 1
                self.ietflg = int(inputFile[l].strip())
                l += 1
                self.iuzfcb1 = int(inputFile[l].strip())
                l += 1
                self.iuzfcb2 = int(inputFile[l].strip())
                l += 1
                self.ntrail2 = int(inputFile[l].strip())
                l += 1
                self.nsets = int(inputFile[l].strip())
                l += 1
                self.nuzgag = int(inputFile[l].strip())
                l += 1
                self.surfdep = float(inputFile[l].strip())
                l += 1
                self.iuzfbnd = [inputFile[l].strip()]
                l += 1
                self.vks = [inputFile[l].strip()]
                l += 1
                self.eps = [inputFile[l].strip()]
                l += 1
                self.thts = [inputFile[l].strip()]
                l += 1
                self.thti = [inputFile[l].strip()]
                self.row_col_iftunit_iuzopt = []
                self.uzfbud_ext = []
                if self.nuzgag > 0:
                    l += 1
                    iuzrow =  inputFile[l].split()
                    l += 1
                    iuzcol = inputFile[l].split()
                    l += 1
                    iftunit = inputFile[l].split()
                    l += 1
                    iuzopt = inputFile[l].split()
                    l += 1
                    for g in range(self.nuzgag):
                        self.row_col_iftunit_iuzopt.append([[int(iuzrow[g]),int(iuzcol[g]),int(iftunit[g]),int(iuzopt[g])]])
                        if iftunit[g]<0.0:
                            self.uzfbud_ext.append(self.ext_uzf + '_tot' + str(abs(int(iftunit[g]))))
                        else:
                            self.uzfbud_ext.append(str(abs(int(iftunit[g]))) + '.' + self.ext_uzf)
                self.finf_user = float(inputFile[l].strip())
            else:
                l += 21
            # wel
            l += 1
            self.wel_yn = int(inputFile[l].strip())
            if self.wel_yn == 1 :
                l += 1
                self.ext_wel = str(inputFile[l].strip())
                l += 1
                self.wel_user = float(inputFile[l].strip())
            else:
                l += 2
            # ghb
            l += 1
            self.ghb_yn = int(inputFile[l].strip())
            if self.ghb_yn == 1 :
                l += 1
                self.ext_ghb = str(inputFile[l].strip())
                l += 1
                self.ghb_head = []
                for i in range(self.nlay):
                    self.ghb_head.append(inputFile[l].split()[i])
                l += 1
                self.ghb_cond = []
                for i in range(self.nlay):
                    self.ghb_cond.append(inputFile[l].split()[i])
            else:
                l += 3
            # drn
            l += 1
            self.drn_yn = int(inputFile[l].strip())
            if self.drn_yn ==1 :
                l += 1
                self.ext_drn = str(inputFile[l].strip())
                l += 1
                self.drn_elev = []
                for i in range(self.nlay):
                    self.drn_elev.append(inputFile[l].split()[i])
                l += 1
                self.drn_cond = []
                for i in range(self.nlay):
                        self.drn_cond.append(inputFile[l].split()[i])
            else:
                l += 3
            # rch
            l += 1
            self.rch_yn = int(inputFile[l].strip())
            if self.rch_yn ==1 :
                l += 1
                self.ext_rch = str(inputFile[l].strip())
                l += 1
                self.nrchop = int(inputFile[l].strip())
                l += 1
                self.rch_user = float(inputFile[l].strip())
            else:
                l += 3
            # MF output
            l += 1
            self.MFout_yn = int(inputFile[l].strip())
        except:
            raise SystemExit("Unexpected error in the MODFLOW input file:\n", sys.exc_info()[0], '\n%s' % traceback.print_exc(file=sys.stdout))
        del inputFile

        self.MM_PROCESS = MMproc.PROCESS(MM_ws                = MM_ws,
                                MF_ws                    = MF_ws,
                                nrow                     = self.nrow,
                                ncol                     = self.ncol,
                                xllcorner                = self.xllcorner,
                                yllcorner                = self.yllcorner,
                                cellsizeMF               = self.delr[0],
                                hnoflo                   = self.hnoflo
                                )

        # 1 - reaf asc file and convert in np.array
        print "\nImporting ESRI ASCII files to initialize the MODFLOW packages"

        self.top     = self.MM_PROCESS.checkarray(self.top)
        self.botm    = self.MM_PROCESS.checkarray(self.botm)
        self.ibound  = self.MM_PROCESS.checkarray(self.ibound, dtype = np.int)
        if self.nlay < 2:
            if isinstance(self.ibound, list):
                self.ibound = (np.asarray(self.ibound)).reshape((self.nrow, self.ncol, 1))
            if isinstance(self.botm, list):
                self.botm = (np.asarray(self.botm)).reshape((self.nrow, self.ncol, 1))
        if self.uzf_yn == 1:
            self.iuzfbnd = self.MM_PROCESS.checkarray(self.iuzfbnd, dtype = np.int)
        if self.ghb_yn == 1:
            ghb = np.asarray(self.MM_PROCESS.checkarray(self.ghb_cond, dtype = np.float))
            if self.nlay > 1:
                self.ghbcells = np.zeros((self.nlay), dtype = np.int)
                for l in range(self.nlay):
                    self.ghbcells[l] = (np.asarray(ghb[:,:,l]) != self.hnoflo).sum()
            else:
                self.ghbcells = (np.asarray(ghb[:,:]) != self.hnoflo).sum()
            del ghb
        if self.drn_yn == 1:
            drn = np.asarray(self.MM_PROCESS.checkarray(self.drn_cond, dtype = np.float))
            if self.nlay > 1:
                self.drncells = np.zeros((self.nlay), dtype = np.int)
                for l in range(self.nlay):
                    self.drncells[l] = (np.asarray(drn[:,:,l]) != self.hnoflo).sum()
            else:
                self.drncells = [(np.asarray(drn[:,:]) != self.hnoflo).sum()]
            del drn
        # TODO confirm wel and add well by user to simulate extraction by borehole
        if self.wel_yn == 1:
            self.welcells = np.zeros((self.nlay), dtype = np.int)
            self.welcells[0] = (np.asarray(self.ibound)[:,:,0] != 0).sum()
            for l in range(1,self.nlay):
                self.welcells[l] = 0

#####################################

    def uzf_obs(self, obs):
        self.nuzgag += len(obs)
        n = 0
        for o in range(len(obs)):
                self.row_col_iftunit_iuzopt.append([[obs.get( obs.keys()[o])['i']+1, obs.get(obs.keys()[o])['j']+1, 200+n, 2]])
                self.uzfbud_ext.append(obs.keys()[o] + '.' + self.ext_uzf)
                n += 1

####################################

    def ppMFtime(self,MM_ws, MF_ws, inputDate_fn, inputZON_dTS_RF_fn, inputZON_dTS_PT_fn, inputZON_dTS_RFe_fn, inputZON_dTS_PE_fn, inputZON_dTS_E0_fn, NMETEO, NVEG, NSOIL):

        ''' RF analysis to define automatically nper/perlen/nstp
        Daily RF>0 creates a nper
        Succesive days with RF=0 are averaged (PT, PE, E0) to a user-input maximum # of days'''

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

        print'\nComputing MODFLOW time discretization based on rainfall analysis in each METEO zone.'
        inputFile = MMproc.readFile(MM_ws, inputDate_fn)
        d = []
        for l in inputFile:
            d.append(l)
        inputFileRF = MMproc.readFile(MM_ws, inputZON_dTS_RF_fn)
        RF_d = np.zeros([NMETEO, len(d)])
        inputFilePT = MMproc.readFile(MM_ws, inputZON_dTS_PT_fn)
        inputFileRFe = MMproc.readFile(MM_ws, inputZON_dTS_RFe_fn)
        PT_d = np.zeros([NMETEO, NVEG, len(d)])
        RFe_d = np.zeros([NMETEO, NVEG, len(d)])
        inputFilePE = MMproc.readFile(MM_ws, inputZON_dTS_PE_fn)
        PE_d = np.zeros([NMETEO, NSOIL, len(d)])
        inputFileE0 = MMproc.readFile(MM_ws, inputZON_dTS_E0_fn)
        E0_d = np.zeros([NMETEO, len(d)])
        for n in range(NMETEO):
            for t in range(len(d)):
                RF_d[n,t] = float(inputFileRF[t+len(d)*n].strip())
                E0_d[n,t] = float(inputFileE0[t+len(d)*n].strip())
            for v in range(NVEG):
                for t in range(len(d)):
                        PT_d[n,v,t] = float(inputFilePT[t+(n*NVEG+v)*len(d)].strip())
                        RFe_d[n,v,t] = float(inputFileRFe[t+(n*NVEG+v)*len(d)].strip())
            for s in range(NSOIL):
                for t in range(len(d)):
                    PE_d[n,s,t] = float(inputFilePE[t+(n*NSOIL+s)*len(d)].strip())
        del inputFile, inputFileRF, inputFilePT, inputFilePE, inputFileRFe, inputFileE0

        if self.timedef == 0:
            # summing RF and avering other flux when RF = 0
            try:
                if isinstance(self.nper, str):
                    perlenmax = int(self.nper.split()[1].strip())
            except:
                raise SystemExit('\nError in nper format of the MODFLOW ini file!\n')
            RF_stp = []
            PT_stp = []
            RFe_stp = []
            PE_stp = []
            E0_stp = []
            RF_stp_tmp = []
            PT_stp_tmp = []
            PE_stp_tmp=[]
            E0_stp_tmp = []
            self.nper = 0
            self.perlen = []
            perlen_tmp = 0
            c = 0
            for n in range(NMETEO):
                RF_stp.append([])
                RF_stp_tmp.append(0.0)
                PT_stp.append([])
                RFe_stp.append([])
                PT_stp_tmp.append([])
                for v in range(NVEG):
                    PT_stp[n].append([])
                    RFe_stp[n].append([])
                    PT_stp_tmp[n].append(0.0)
                PE_stp.append([])
                PE_stp_tmp.append([])
                for v in range(NSOIL):
                    PE_stp[n].append([])
                    PE_stp_tmp[n].append(0.0)
                E0_stp.append([])
                E0_stp_tmp.append(0.0)
            if perlenmax < 2:
                raise SystemExit('\nperlenmax must be higher than 1!\nCorrect perlenmax in the MODFLOW ini file or select the daily option.')
            for j in range(len(d)):
                    if RFe_d[:,:,j].sum()>0.0:
                        if c == 1:
                            for n in range(NMETEO):
                                RF_stp[n].append(RF_stp_tmp[n]/perlen_tmp)
                                RF_stp_tmp[n] = 0.0
                                for v in range(NVEG):
                                    PT_stp[n][v].append(PT_stp_tmp[n][v]/perlen_tmp)
                                    RFe_stp[n][v].append(0.0)
                                    PT_stp_tmp[n][v] = 0.0
                                for s in range(NSOIL):
                                    PE_stp[n][s].append(PE_stp_tmp[n][s]/perlen_tmp)
                                    PE_stp_tmp[n][s] = 0.0
                                E0_stp[n].append(E0_stp_tmp[n]/perlen_tmp)
                                E0_stp_tmp[n] = 0.0
                            self.perlen.append(perlen_tmp)
                            perlen_tmp = 0
                        for n in range(NMETEO):
                            RF_stp[n].append(RF_d[n,j])
                            for v in range(NVEG):
                                PT_stp[n][v].append(PT_d[n,v,j])
                                RFe_stp[n][v].append(RFe_d[n,v,j])
                            for s in range(NSOIL):
                                PE_stp[n][s].append(PE_d[n,s,j])
                            E0_stp[n].append(E0_d[n,j])
                        self.nper += 1
                        self.perlen.append(1)
                        c = 0
                    else:
                        if perlen_tmp < perlenmax:
                            for n in range(NMETEO):
                                RF_stp_tmp[n]  += RF_d[n,j]
                                for v in range(NVEG):
                                    PT_stp_tmp[n][v] += PT_d[n,v,j]
                                for s in range(NSOIL):
                                    PE_stp_tmp[n][s]  += PE_d[n,s,j]
                                E0_stp_tmp[n]  += E0_d[n,j]
                            if c == 0:
                                self.nper += 1
                            perlen_tmp += 1
                            c = 1
                        else:
                            for n in range(NMETEO):
                                RF_stp[n].append(RF_stp_tmp[n]/perlen_tmp)
                                RF_stp_tmp[n] = 0.0
                                for v in range(NVEG):
                                    PT_stp[n][v].append(PT_stp_tmp[n][v]/perlen_tmp)
                                    RFe_stp[n][v].append(0.0)
                                    PT_stp_tmp[n][v] = 0.0
                                for s in range(NSOIL):
                                    PE_stp[n][s].append(PE_stp_tmp[n][s]/perlen_tmp)
                                    PE_stp_tmp[n][s] = 0.0
                                E0_stp[n].append(E0_stp_tmp[n]/perlen_tmp)
                                E0_stp_tmp[n] = 0.0
                            self.perlen.append(perlen_tmp)
                            for n in range(NMETEO):
                                RF_stp_tmp[n]  += RF_d[n,j]
                                for v in range(NVEG):
                                    PT_stp_tmp[n][v] += PT_d[n,v,j]
                                for s in range(NSOIL):
                                    PE_stp_tmp[n][s]  += PE_d[n,s,j]
                                E0_stp_tmp[n]  += E0_d[n,j]
                            self.nper += 1
                            perlen_tmp = 1
                            c = 1
            if c == 1:
                for n in range(NMETEO):
                    RF_stp[n].append(RF_stp_tmp[n]/perlen_tmp)
                    RF_stp_tmp[n] = 0.0
                    for v in range(NVEG):
                        PT_stp[n][v].append(PT_stp_tmp[n][v]/perlen_tmp)
                        RFe_stp[n][v].append(0.0)
                        PT_stp_tmp[n][v] = 0.0
                    for s in range(NSOIL):
                        PE_stp[n][s].append(PE_stp_tmp[n][s]/perlen_tmp)
                        PE_stp_tmp[n][s] = 0.0
                    E0_stp[n].append(E0_stp_tmp[n]/perlen_tmp)
                    E0_stp_tmp[n] = 0.0
                self.perlen.append(perlen_tmp)
            del perlen_tmp
            del RF_d, RFe_d, PT_d, PE_d, E0_d, RF_stp_tmp, PT_stp_tmp, PE_stp_tmp, E0_stp_tmp, c
            self.perlen = np.asarray(self.perlen, dtype = np.int)
            self.nstp = np.ones(self.nper, dtype = np.int)
            self.tsmult = self.nstp
            self.Ss_tr = []
            for n in range(self.nper):
                self.Ss_tr.append(False)
        elif self.timedef>0:
            RF_stp = np.zeros([NMETEO,sum(self.nstp)])
            PT_stp = np.zeros([NMETEO,NVEG,sum(self.nstp)])
            RFe_stp = np.zeros([NMETEO,NVEG,sum(self.nstp)])
            PE_stp = np.zeros([NMETEO,NSOIL,sum(self.nstp)])
            E0_stp = np.zeros([NMETEO,sum(self.nstp)])
            stp = 0
            for per in range(self.nper):
                for stp in range(self.nstp[per]):
                    tstart = sum(self.perlen[0:per])+stp*(self.perlen[per]/self.nstp[per])
                    tend = sum(self.perlen[0:per])+(self.perlen[per]/self.nstp[per])*(1+stp)
                    for n in range(NMETEO):
                        RF_stp[n,sum(self.nstp[0:per])+stp] += RF_d[n,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                        for v in range(NVEG):
                            PT_stp[n,v,sum(nstp[0:per])+stp] += PT_d[n,v,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                            RFe_stp[n,v,sum(nstp[0:per])+stp] += RFe_d[n,v,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                        for s in range(NSOIL):
                            PE_stp[n,s,sum(nstp[0:per])+stp] += PE_d[n,s,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                        E0_stp[n,sum(nstp[0:per])+stp] += E0_d[n,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])

        self.inputZON_TS_RF_fn = "inputZONRF_stp.txt"
        self.inputZONRF_fn = os.path.join(MM_ws, self.inputZON_TS_RF_fn)
        inputZONRF = open(self.inputZONRF_fn, 'w')
        inputZONRF.write('#\n')

        self.inputZON_TS_PT_fn = "inputZONPT_stp.txt"
        self.inputZONPT_fn = os.path.join(MM_ws, self.inputZON_TS_PT_fn)
        inputZONPT = open(self.inputZONPT_fn, 'w')
        inputZONPT.write('#\n')

        self.inputZON_TS_RFe_fn = "inputZONRFe_stp.txt"
        self.inputZONRFe_fn = os.path.join(MM_ws, self.inputZON_TS_RFe_fn)
        inputZONRFe = open(self.inputZONRFe_fn, 'w')
        inputZONRFe.write('#\n')

        self.inputZON_TS_PE_fn = "inputZONPE_stp.txt"
        self.inputZONPE_fn = os.path.join(MM_ws, self.inputZON_TS_PE_fn)
        inputZONPE = open(self.inputZONPE_fn, 'w')
        inputZONPE.write('#\n')

        self.inputZON_TS_E0_fn = "inputZONE0_stp.txt"
        self.inputZONE0_fn = os.path.join(MM_ws, self.inputZON_TS_E0_fn)
        inputZONE0 = open(self.inputZONE0_fn, 'w')
        inputZONE0.write('#\n')

        for n in range(NMETEO):
            try:
                ExportResults1(RF_stp[n], inputZONRF)
                ExportResults1(E0_stp[n], inputZONE0)
                #VEG
                if NVEG>0:
                    for v in range(NVEG):
                        ExportResults1(PT_stp[n][v], inputZONPT)
                        ExportResults1(RFe_stp[n][v], inputZONRFe)
                # SOIL
                if NSOIL>0:
                    for s in range(NSOIL):
                        ExportResults1(PE_stp[n][s], inputZONPE)
            except:
                print "\nError in output file, some output files were not exported."

        print 'Found %d days converted into %d stress periods.' % (sum(self.perlen), self.nper)

        inputZONRF.close()
        inputZONRFe.close()
        inputZONPT.close()
        inputZONPE.close()
        inputZONE0.close()

#####################################

    def ppMF(self, MM_ws, MF_ws, MF_ini_fn, finf_MM = "", finf_user = None, wel_MM = "", wel_user = None, report = None, verbose = 1, chunks = 0, numDays = -1):

        if verbose == 0:
            print '--------------'

        if os.path.exists(finf_MM[0]):
            if self.uzf_yn == 1:
                finf_input = finf_MM
            if self.wel_yn == 1:
                wel_input = wel_MM
            if self.rch_yn == 1:
                rch_input = rch_MM
        else:
            if self.uzf_yn == 1:
                finf_input = self.finf_user
            if self.wel_yn == 1:
                wel_input = self.wel_user
            if self.rch_yn == 1:
                rch_input = self.rch_user

        # 1 - reaf asc file and convert in np.array
        print "\nImporting ESRI ASCII files to initialize the MODFLOW packages"

        strt    = self.MM_PROCESS.checkarray(self.strt)
        hk      = self.MM_PROCESS.checkarray(self.hk)
        vka     = self.MM_PROCESS.checkarray(self.vka)
        ss      = self.MM_PROCESS.checkarray(self.ss)
        sy      = self.MM_PROCESS.checkarray(self.sy)
        if self.uzf_yn == 1:
            if self.iuzfopt == 1:
                vks     = self.MM_PROCESS.checkarray(self.vks)
            else:
                vks = 0.0
            eps     = self.MM_PROCESS.checkarray(self.eps)
            thts    = self.MM_PROCESS.checkarray(self.thts)
            thti    = self.MM_PROCESS.checkarray(self.thti)

        # FINF
        if self.uzf_yn == 1:
            print '\nUZF1 package initialization'
            if isinstance(finf_input,float):
                finf_array = finf_input
                print 'infiltration input: %s' % str(finf_input)
            else:
                finf_array = []
                print 'infiltration input: %s' % finf_input[0]
                try:
                    h5_finf = h5py.File(finf_input[0], 'r')
                    for n in range(self.nper):
                        finf_array.append(h5_finf[finf_input[1]][n])
                    h5_finf.close()
                except:
                    finf_array = self.finf_user
                    print 'WARNING!\nNo valid UZF1 package file(s) provided, running MODFLOW using user-input UZF1 infiltration value: %.3G' % self.finf_user
                    finf_input = self.finf_user

        # WELL
        # TODO confirm wel and add well by user to simulate extraction by borehole
        if self.wel_yn == 1:
            print '\nWEL package initialization'
            if isinstance(wel_input,float):
                wel_array = wel_input
                print 'discharge input: %s' % str(wel_input)
            else:
                wel_array = []
                print 'discharge input: %s' % wel_input[0]
                try:
                    h5_wel = h5py.File(wel_input[0], 'r')
                    for n in range(self.nper):
                        wel_array.append(h5_wel[wel_input[1]][n])
                    h5_wel.close()
                except:
                    wel_array = self.wel_user
                    print 'WARNING!\nNo valid WEL package file(s) provided, running MODFLOW using user-input well value: %.3G' % self.wel_user
                    wel_input = self.wel_user

        # RCH
        if self.rch_yn == 1:
            print '\nRCH package initialization'
            if isinstance(rch_input,float):
                rch_array = rch_input
                print 'recharge input: %s' % str(rch_input)
            else:
                rch_array = []
                print 'recharge input: %s' % rch_input[0]
                try:
                    h5_rch = h5py.File(rch_input[0], 'r')
                    for n in range(self.nper):
                        rch_array.append(h5_rch[rch_input[1]][n])
                    h5_rch.close()
                except:
                    rch_array = rch_user
                    print 'WARNING!\nNo valid RCH package file(s) provided, running MODFLOW using user-input recharge value: %.3G' % rch_user
                    rch_input = rch_user

        # DRAIN
        if self.drn_yn == 1:
            print '\nDRN package initialization'
            l = 0
            layer_row_column_elevation_cond = [[]]
            for d in self.drn_cond:
                drn_check = 1
                drn_elev_array = np.zeros((self.nrow,self.ncol))
                drn_cond_array = np.zeros((self.nrow,self.ncol))
                if isinstance(d, str):
                    drn_elev_path = os.path.join(MF_ws, self.drn_elev[l])
                    drn_elev_array[:,:] = self.MM_PROCESS.convASCIIraster2array(drn_elev_path, drn_elev_array[:,:])
                    drn_cond_path = os.path.join(MF_ws, self.drn_cond[l])
                    drn_cond_array[:,:] = self.MM_PROCESS.convASCIIraster2array(drn_cond_path, drn_cond_array[:,:])
                else:
                    drn_elev_array[:,:] = self.drn_elev[l]
                    drn_cond_array[:,:] = self.drn_cond[l]
                for i in range(self.nrow):
                    for j in range(self.ncol):
                        if drn_elev_array[i,j]<>0:
                            if drn_elev_array[i,j]<0:
                                if isinstance(self.botm[l], float):
                                    drn_elev_tmp = self.botm[l]
                                else:
                                    drn_elev_tmp = self.botm[i][j][l] #- (botm_array[i][j][l]-botm_array[i][j][l-1])/10
                            else:
                                drn_elev_tmp = drn_elev_array[i,j]
                            layer_row_column_elevation_cond[0].append([l+1, i+1, j+1, drn_elev_tmp, drn_cond_array[i,j]])
                l += 1

        # GHB
        if self.ghb_yn == 1:
            print '\nGHB package initialization'
            l = 0
            layer_row_column_head_cond = [[]]
            for d in self.ghb_cond:
                ghb_check = 1
                ghb_head_array = np.zeros((self.nrow,self.ncol))
                ghb_cond_array = np.zeros((self.nrow,self.ncol))
                if isinstance(d, str):
                    ghb_head_path = os.path.join(MF_ws, self.ghb_head[l])
                    ghb_head_array[:,:] = self.MM_PROCESS.convASCIIraster2array(ghb_head_path, ghb_head_array[:,:])
                    ghb_cond_path = os.path.join(MF_ws, self.ghb_cond[l])
                    ghb_cond_array[:,:] = self.MM_PROCESS.convASCIIraster2array(ghb_cond_path, ghb_cond_array[:,:])
                else:
                    ghb_head_array[:,:] = self.ghb_head[l]
                    ghb_cond_array[:,:] = self.ghb_cond[l]
                for i in range(self.nrow):
                    for j in range(self.ncol):
                        if ghb_head_array[i,j]<>0:
                            ghb_head_tmp = ghb_head_array[i,j]
                            layer_row_column_head_cond[0].append([l+1, i+1, j+1, ghb_head_tmp, ghb_cond_array[i,j]])
                l += 1

        # average for 1st SS stress period
        self.nper        = self.nper
        self.perlen      = list(self.perlen)
        self.nstp        = list(self.nstp)
        self.tsmult      = list(self.tsmult)
        self.Ss_tr       = list(self.Ss_tr)
        if self.dum_sssp1 == 1:
            if self.uzf_yn == 1 and isinstance(finf_input,tuple):
                finf_array = np.asarray(finf_array)
                finf_SS = np.zeros((self.nrow,self.ncol))
                for n in range(self.nper):
                    finf_SS += finf_array[n,:,:]
                finf_SS = finf_SS/self.nper
                finf_array = list(finf_array)
                finf_array.insert(0, finf_SS)
            if self.rch_yn == 1 and isinstance(rch_input,tuple):
                rch_array = np.asarray(rch_array)
                rch_SS = np.zeros((self.nrow,self.ncol))
                for n in range(self.nper):
                    rch_SS += rch_array[n,:,:]
                rch_SS = rch_SS/self.nper
                rch_array = list(rch_array)
                rch_array.insert(0, rch_SS)
            if self.wel_yn == 1 and isinstance(wel_input,tuple):
                wel_array = np.asarray(wel_array)
                wel_SS = np.zeros((self.nrow,self.ncol))
                for n in range(self.nper):
                    wel_SS += wel_array[n,:,:]
                wel_SS = wel_SS/self.nper
                wel_array = list(wel_array)
                wel_array.insert(0, wel_SS)
            self.nper +=  1
            self.perlen.insert(0,1)
            self.nstp.insert(0,1)
            self.tsmult.insert(0,1)
            self.Ss_tr.insert(0, True)

        # WEL
        if self.wel_yn == 1:
            layer_row_column_Q = []
            for n in range(self.nper):
                layer_row_column_Q.append([])
                for r in range(self.nrow):
                    for c in range(self.ncol):
                        if np.abs(self.ibound[r][c][:]).sum() != 0:
                            if isinstance(wel_array, float):
                                layer_row_column_Q[n].append([1,r+1,c+1,-wel_array*self.delr[c]*self.delc[r]])
                            else:
                                layer_row_column_Q[n].append([1,r+1,c+1,-(wel_array[n][r][c])*self.delr[c]*self.delc[r]])

        # 2 - create the modflow packages files
        print '\nMODFLOW files writing'
        ncells_package = []
        # MFfile initialization
        mfmain = mf.modflow(modelname = self.modelname, exe_name = self.exe_name, namefile_ext = self.namefile_ext, version = self.version, model_ws = MF_ws)
        # dis package
        dis = mf.mfdis(model = mfmain, nrow = self.nrow, ncol = self.ncol, nlay = self.nlay, nper = self.nper, delr = self.delr, delc = self.delc, laycbd = self.laycbd, top = self.top, botm = self.botm, perlen = self.perlen, nstp = self.nstp, tsmult = self.tsmult, itmuni = self.itmuni, lenuni = self.lenuni, steady = self.Ss_tr, extension = self.ext_dis)
        dis.write_file()
        # bas package
        bas = mf.mfbas(model = mfmain, ibound = self.ibound, strt = strt, hnoflo = self.hnoflo, extension = self.ext_bas)
        bas.write_file()
        # layer package
        if self.version != 'mfnwt':
            # lpf package
            lpf = mf.mflpf(model = mfmain, hdry = self.hdry, laytyp = self.laytyp, layavg = self.layavg, chani = self.chani, layvka = self.layvka, laywet = self.laywet, hk = hk, vka = vka, ss = ss, sy = sy, storagecoefficient = self.storagecoefficient, constantcv = self.constantcv, thickstrt = self.thickstrt, nocvcorrection = self.nocvcorrection, novfc = self.novfc, extension = self.ext_lpf)
            lpf.write_file()
            cb = lpf.ilpfcb
        else:
            upw = mf.mfupw(model = mfmain, hdry = self.hdry, iphdry = self.iphdry, laytyp = self.laytyp, layavg = self.layavg, chani = self.chani, layvka = self.layvka, laywet = self.laywet, hk = hk, vka = vka, ss = ss, sy = sy, extension = self.ext_upw)
            upw.write_file()
            cb = upw.iupwcb
        # wel package
        if self.wel_yn == 1:
            if layer_row_column_Q <> None:
                wel = mf.mfwel(model = mfmain, iwelcb = cb, layer_row_column_Q = layer_row_column_Q, extension = self.ext_wel)
                wel.write_file()
        # drn package
        if self.drn_yn == 1:
            if drn_check == 1:
                drn = mf.mfdrn(model = mfmain, idrncb = cb, layer_row_column_elevation_cond = layer_row_column_elevation_cond, extension = self.ext_drn)
                drn.write_file()
        # ghb package
        if self.ghb_yn == 1:
            ghb = mf.mfghb(model = mfmain, ighbcb = cb, layer_row_column_head_cond = layer_row_column_head_cond, extension = self.ext_ghb)
            ghb.write_file()
        # uzf package
        if self.uzf_yn == 1:
            uzf = mf.mfuzf1(model = mfmain, nuztop = self.nuztop, iuzfopt = self.iuzfopt, irunflg = self.irunflg, ietflg = self.ietflg, iuzfcb1 = self.iuzfcb1, iuzfcb2 = self.iuzfcb2, ntrail2 = self.ntrail2, nsets = self.nsets, nuzgag = self.nuzgag, surfdep = self.surfdep, iuzfbnd = self.iuzfbnd, vks = vks, eps = eps, thts = thts, thti = thti, row_col_iftunit_iuzopt = self.row_col_iftunit_iuzopt, finf = finf_array, extension = self.ext_uzf, uzfbud_ext = self.uzfbud_ext)
            uzf.write_file()
        # rch package
        if self.rch_yn == 1:
            rch = mf.mfrch(model = mfmain, irchcb = cb, nrchop = self.nrchop, rech = rch_array, extension = self.ext_rch)
            rch.write_file()
        # output control package
        oc = mf.mfoc(model = mfmain, ihedfm = self.ihedfm, iddnfm = self.iddnfm, item2 = [[0,1,1,1]], item3 = [[0,0,1,0]], extension = [self.ext_oc,self.ext_heads,self.ext_ddn,self.ext_cbc])
        oc.write_file()
        # solver
        if self.version != 'mfnwt':
            # pcg package
            pcg = mf.mfpcg(model = mfmain, mxiter = 150, iter1=75, hclose = self.hclose, rclose = self.rclose, npcond = 1, relax = 0.99, nbpol=0, iprpcg = 0, mutpcg = 1, damp = 0.6, extension = self.ext_pcg)
            pcg.write_file()
            # sip
        #    sip = mf.mfsip(mfmain, hclose=1e-3)
        #    sip.write_file()
            # sor
        #    sor = mf.mfsor(mfmain, hclose=1e-3)
        #    sor.write_file()
        else:
            nwt = mf.mfnwt(model = mfmain, headtol = self.headtol, fluxtol = self.fluxtol, maxiterout = self.maxiterout, thickfact = self.thickfact, linmeth = self.linmeth, iprnwt = self.iprnwt, ibotav = self.ibotav, options = self.options, extension = self.ext_nwt)
            nwt.write_file()

        self.h_MF_fn = os.path.join(MF_ws, self.modelname + "." + self.ext_heads)
        self.cbc_MF_fn = os.path.join(MF_ws, self.modelname + "." + self.ext_cbc)
        if self.uzf_yn == 1:
            self.cbc_MFuzf_fn = os.path.join(MF_ws, self.modelname + ".uzfbt1")

        # run MODFLOW and read the heads back into Python
        print '\nMODFLOW run'
        mfmain.write_name_file()
        mfmain.run_model(pause = False, report = report)

        def readh():
            """
            Extract heads
            """
            try:
                h = mfrdbin.mfhdsread(mfmain, 'LF95').read_all(self.h_MF_fn)
            except:
                h5_MF.close()
                raise ValueError, '\nMODFLOW error!\nCheck the MODFLOW list file in folder:\n%s' % MF_ws
            if len(h[1])<sum(self.nstp):
                h5_MF.close()
                raise ValueError, '\nMODFLOW error!\nCheck the MODFLOW list file in folder:\n%s' % MF_ws
            print ''
            return h

        def readcbc():
            """
            Extract cell-by-cell budget
            """
            cbc = mfrdbin.mfcbcread(mfmain, 'LF95').read_all(self.cbc_MF_fn)
            h5_MF.create_dataset('cbc_nam', data = np.asarray(cbc[2][1]))
            print ''
            return cbc

        def readcbcuzf():
            """
            Extract cell-by-cell uzf budget
            """
            if self.uzf_yn == 1:
                cbc_uzf = mfrdbin.mfcbcread(mfmain, 'LF95').read_all(self.cbc_MFuzf_fn)
                h5_MF.create_dataset('cbc_uzf_nam', data = np.asarray(cbc_uzf[2][1]))
                return cbc_uzf

        h5_MF = h5py.File(self.h5_MF_fn, 'w')
        print '\nPlease wait. Storing heads and cbc terms into HDF5 file\n%s\n' % (self.h5_MF_fn)
        if self.dum_sssp1 == 1:
            if chunks == 1:
                # TODO implement swapaxes
                h = readh()
                h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1][1:]),   chunks = (1,self.nrow,self.ncol,self.nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf
                del h
                cbc = readcbc()
                h5_MF.create_dataset(name = 'cbc',   data = np.asarray(cbc[1][1:]), chunks = (1,len(h5_MF['cbc_nam']),self.nrow,self.ncol,self.nlay),   compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
                del cbc
                if self.uzf_yn == 1:
                    cbc_uzf = readcbcuzf()
                    h5_MF.create_dataset(name = 'cbc_uzf', data = np.asarray(cbc_uzf[1][1:]), chunks = (1,len(h5_MF['cbc_uz_nam']),self.nrow,self.ncol,self.nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
                    del cbc_uzf
            else:
                h = readh()
                h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1][1:]))
                del h
                cbc = readcbc()
                data = np.asarray(cbc[1][1:])
                del cbc
                data = np.swapaxes(data,1,2)
                data = np.swapaxes(data,2,3)
                h5_MF.create_dataset(name = 'cbc',   data = data)
                del data
                if self.uzf_yn == 1:
                    cbc_uzf = readcbcuzf()
                    data = np.asarray(cbc_uzf[1][1:])
                    del cbc_uzf
                    data = np.swapaxes(data,1,2)
                    data = np.swapaxes(data,2,3)
                    h5_MF.create_dataset(name = 'cbc_uzf', data = data)
                    del data
            self.nper = self.nper - 1
            self.perlen = self.perlen[1:]
            self.tsmult = self.tsmult[1:]
            self.nstp = self.nstp[1:]
            self.Ss_tr = self.Ss_tr[1:]
        else:
            if chunks == 1:
                # TODO implement swapaxes
                h = readh()
                h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1]), chunks = (1,self.nrow,self.ncol,self.nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
                del h
                cbc = readcbc()
                h5_MF.create_dataset(name = 'cbc', data = np.asarray(cbc[1]), chunks = (1,len(h5_MF['cbc_nam']),self.nrow,self.ncol,self.nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
                del cbc
                if self.uzf_yn == 1:
                    cbc_uzf = readcbcuzf()
                    h5_MF.create_dataset(name = 'cbc_uzf', data = np.asarray(cbc_uzf[1]), chunks = (1,len(h5_MF['cbc_uzf_nam']),self.nrow,self.ncol,self.nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
                    del cbc_uzf
            else:
                h = readh()
                h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1]))
                del h
                cbc = readcbc()
                data = np.asarray(cbc[1])
                del cbc
                data = np.swapaxes(data,1,2)
                data = np.swapaxes(data,2,3)
                h5_MF.create_dataset(name = 'cbc', data = data)
                del data
                if self.uzf_yn == 1:
                    cbc_uzf = readcbcuzf()
                    data = np.asarray(cbc_uzf[1])
                    del cbc_uzf
                    data = np.swapaxes(data,1,2)
                    data = np.swapaxes(data,2,3)
                    h5_MF.create_dataset(name = 'cbc_uzf', data = data)
                    del data
        h4MM = np.zeros((len(self.perlen),self.nrow,self.ncol), dtype = np.float)
        h_MF = h5_MF['heads'][:,:,:,:]
        iuzfbnd = np.asarray(self.iuzfbnd)
        for l in range(self.nlay):
            for t in range(len(self.perlen)):
                mask = np.ma.make_mask(iuzfbnd == l+1)
                h4MM[t,:,:] += h_MF[t,:,:,l]*mask
        h5_MF.create_dataset(name = 'heads4MM', data = h4MM)
        exf4MM = np.zeros((len(self.perlen),self.nrow,self.ncol), dtype = np.float)
        if self.uzf_yn == 1:
            cbc_uzf_nam = []
            for c in h5_MF['cbc_uzf_nam']:
                cbc_uzf_nam.append(c.strip())
            imfEXF   = cbc_uzf_nam.index('SURFACE LEAKAGE')
            exf_MF = h5_MF['cbc_uzf'][:,:,:,imfEXF,:]
            for l in range(self.nlay):
                for t in range(len(self.perlen)):
                    mask = np.ma.make_mask(iuzfbnd == l+1)
                    exf4MM[t,:,:] += exf_MF[t,:,:,l]*mask
            h5_MF.create_dataset(name = 'exf4MM', data = exf4MM)
        h5_MF.close()
        # to delete MF binary files and save disk space
        if self.MFout_yn == 0:
            os.remove(h_MF_fn)
            os.remove(cbc_MF_fn)
            os.remove(cbc_MFuzf_fn)