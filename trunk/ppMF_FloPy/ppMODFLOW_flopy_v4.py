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
__version__ = "0.4"
__date__ = "2012"

import sys, os, traceback
import numpy as np
import matplotlib as mpl
import h5py
import mf
import mfreadbinaries as mfrdbin
import MARMITESprocess_v4 as MMproc

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
                        raise SystemExit("\nFATAL ERROR!\nMM doesn't accept nstp < perlen!")
                    if self.Ss_tr[i] == 'TR':
                        self.Ss_tr[i] = False
                    elif self.Ss_tr[i] == 'SS':
                        self.Ss_tr[i] = True
                    else:
                        raise SystemExit('\nFATAL ERROR!\nVariable Ss_tr from the DIS package is not correct, check the MODFLOW manual')
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
            self.thick = []
            for i in range(self.nlay):
                self.thick.append(inputFile[l].split()[i])
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
            raise SystemExit("\nFATAL ERROR!\nUnexpected error in the MODFLOW input file:\n", sys.exc_info()[0], '\n%s' % traceback.print_exc(file=sys.stdout))
        del inputFile

        self.MM_PROCESS = MMproc.PROCESS(MM_ws           = MM_ws,
                                MF_ws                    = MF_ws,
                                nrow                     = self.nrow,
                                ncol                     = self.ncol,
                                xllcorner                = self.xllcorner,
                                yllcorner                = self.yllcorner,
                                cellsizeMF               = self.delr[0],
                                hnoflo                   = self.hnoflo
                                )

        # 1 - reaf asc file and convert in np.array
        print "\nImporting ESRI ASCII files to initialize the MODFLOW packages..."

        self.top     = self.MM_PROCESS.checkarray(self.top)
        self.strt    = self.MM_PROCESS.checkarray(self.strt)
        self.thick   = self.MM_PROCESS.checkarray(self.thick)
        if self.nlay < 2:
            if isinstance(self.thick, list):
                self.thick = (np.asarray(self.thick)).reshape((self.nrow, self.ncol, 1))
        else:
            self.thick = np.asarray(self.thick)
        top_tmp = np.asarray(self.top)
        botm_tmp = []
        for l in range(self.nlay):
            if l == 0:
                thick_tmp = self.thick[:,:,l]
            else:
                thick_tmp += self.thick[:,:,l]
            botm_tmp.append(top_tmp - thick_tmp)
        del self.thick, top_tmp, thick_tmp
        self.botm = list(np.swapaxes(botm_tmp,0,1))
        self.botm = list(np.swapaxes(self.botm,1,2))
        del botm_tmp
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
                    self.ghbcells[l] = (np.asarray(ghb[:,:,l]) > 0.0).sum()
            else:
                self.ghbcells = (np.asarray(ghb[:,:]) > 0.0).sum()
            del ghb
        if self.drn_yn == 1:
            drn = np.asarray(self.MM_PROCESS.checkarray(self.drn_cond, dtype = np.float))
            if self.nlay > 1:
                self.drncells = np.zeros((self.nlay), dtype = np.int)
                for l in range(self.nlay):
                    self.drncells[l] = (np.asarray(drn[:,:,l]) > 0.0).sum()
            else:
                self.drncells = [(np.asarray(drn[:,:]) > 0.0).sum()]
            del drn

#####################################

    def uzf_obs(self, obs):
        self.nuzgag += len(obs)
        n = 0
        for o in range(len(obs)):
                self.row_col_iftunit_iuzopt.append([[obs.get( obs.keys()[o])['i']+1, obs.get(obs.keys()[o])['j']+1, 200+n, 2]])
                self.uzfbud_ext.append(obs.keys()[o] + '.' + self.ext_uzf)
                n += 1

####################################

    def ppMFtime(self, MM_ws, MF_ws, inputDate_fn, inputZON_dSP_RF_veg_fn, inputZON_dSP_RFe_veg_fn, inputZON_dSP_PT_fn, input_dSP_LAI_veg_fn, inputZON_dSP_PE_fn, inputZON_dSP_E0_fn, NMETEO, NVEG, NSOIL, inputZON_dSP_RF_irr_fn = None, inputZON_dSP_RFe_irr_fn = None, inputZON_dSP_PT_irr_fn = None, input_dSP_crop_irr_fn = None, NFIELD = None):

        ''' RF analysis to define automatically nper/perlen/nstp
        Daily RF>0 creates a nper
        Succesive days with RF=0 are averaged (PT, PE, E0) to a user-input maximum # of days'''

        def ExportResults1(SP, outFileExport):
            """
            Write the processed data in a open txt file readable by  MARMITES
            INPUT:      output flux time series and open file
            """
            for i in range(len(SP)):
                out_line =  '%14.9G' %SP[i],'\n'
                for l in out_line:
                    outFileExport.write(l)

        #####################################

        print'\nComputing MODFLOW time discretization based on rainfall analysis in each METEO zone...'
        # READ date of input files (RF and PT)
        inputDate_fn=os.path.join(MM_ws, inputDate_fn)
        if os.path.exists(inputDate_fn):
            inputDate_tmp = np.loadtxt(inputDate_fn, dtype = str)
            self.inputDate = inputDate_tmp[:,0]
            self.JD = np.asarray(inputDate_tmp[:,2], dtype = np.int)
            del inputDate_tmp
            self.inputDate = mpl.dates.datestr2num(self.inputDate)
            for i in range(1,len(self.inputDate)):
                #__________________Check date consistency________________#
                difDay=self.inputDate[i]-self.inputDate[i-1]
                if (difDay !=1.0):
                    print 'DifDay = ' + str(difDay)
                    raise SystemExit, '\nFATAL ERROR!\nDates of the input data are not sequencial, check your daily time step!\nError at date %s ' % str(self.inputDate[i])
        else:
            raise ValueError, "\nFATAL ERROR!\nThe file %s doesn't exist!!!" % inputDate_fn

        inputFileRF_veg = MMproc.readFile(MM_ws, inputZON_dSP_RF_veg_fn)
        RF_veg_d = np.zeros([NMETEO, len(self.JD)])
        inputFileRFe_veg = MMproc.readFile(MM_ws, inputZON_dSP_RFe_veg_fn)
        RFe_veg_d = np.zeros([NMETEO, NVEG, len(self.JD)])
        inputFilePT = MMproc.readFile(MM_ws, inputZON_dSP_PT_fn)
        PT_veg_d = np.zeros([NMETEO, NVEG, len(self.JD)])
        inputFileLAI_veg = MMproc.readFile(MM_ws, input_dSP_LAI_veg_fn)
        self.LAI_veg_d = np.zeros([NMETEO, NVEG, len(self.JD)], dtype = float)
        inputFilePE = MMproc.readFile(MM_ws, inputZON_dSP_PE_fn)
        PE_d = np.zeros([NMETEO, NSOIL, len(self.JD)])
        inputFileE0 = MMproc.readFile(MM_ws, inputZON_dSP_E0_fn)
        E0_d = np.zeros([NMETEO, len(self.JD)])
        if NFIELD <> None:
            inputFileRF_irr = MMproc.readFile(MM_ws, inputZON_dSP_RF_irr_fn)
            RF_irr_d = np.zeros([NMETEO, NFIELD, len(self.JD)], dtype = float)
            inputFileRFe_irr = MMproc.readFile(MM_ws, inputZON_dSP_RFe_irr_fn)
            RFe_irr_d = np.zeros([NMETEO, NFIELD, len(self.JD)], dtype = float)
            inputFilePT_irr = MMproc.readFile(MM_ws, inputZON_dSP_PT_irr_fn)
            PT_irr_d = np.zeros([NMETEO, NFIELD, len(self.JD)], dtype = float)
            inputFilecrop_irr = MMproc.readFile(MM_ws, input_dSP_crop_irr_fn)
            self.crop_irr_d = np.zeros([NMETEO, NFIELD, len(self.JD)], dtype = float)
        for n in range(NMETEO):
            for t in range(len(self.JD)):
                RF_veg_d[n,t] = float(inputFileRF_veg[t+len(self.JD)*n].strip())
                E0_d[n,t] = float(inputFileE0[t+len(self.JD)*n].strip())
                for v in range(NVEG):
                    PT_veg_d[n,v,t] = float(inputFilePT[t+(n*NVEG+v)*len(self.JD)].strip())
                    RFe_veg_d[n,v,t] = float(inputFileRFe_veg[t+(n*NVEG+v)*len(self.JD)].strip())
                    self.LAI_veg_d[n,v,t] = float(inputFileLAI_veg[t+v*len(self.JD)].strip())
                for s in range(NSOIL):
                    PE_d[n,s,t] = float(inputFilePE[t+(n*NSOIL+s)*len(self.JD)].strip())
                if NFIELD <> None:
                    for f in range(NFIELD):
                        RF_irr_d[n,f,t] = float(inputFileRF_irr[t+(n*NFIELD+f)*len(self.JD)].strip())
                        RFe_irr_d[n,f,t] = float(inputFileRFe_irr[t+(n*NFIELD+f)*len(self.JD)].strip())
                        PT_irr_d[n,f,t] = float(inputFilePT_irr[t+(n*NFIELD+f)*len(self.JD)].strip())
                        self.crop_irr_d[n,f,t] = float(inputFilecrop_irr[t+f*len(self.JD)].strip())
        del inputDate_fn, inputFileRF_veg, inputFilePT, inputFilePE, inputFileRFe_veg, inputFileE0
        if NFIELD <> None:
            del inputZON_dSP_RF_irr_fn, inputZON_dSP_RFe_irr_fn, inputZON_dSP_PT_irr_fn, input_dSP_crop_irr_fn, inputFileRF_irr, inputFileRFe_irr, inputFilePT_irr, inputFilecrop_irr

        if self.timedef == 0:
            # summing RF and avering other flux when RF = 0
            try:
                if isinstance(self.nper, str):
                    perlenmax = int(self.nper.split()[1].strip())
            except:
                raise SystemExit('\nFATAL ERROR!\nError in nper format of the MODFLOW ini file!\n')
            RF_veg_stp = []
            RFe_veg_stp = []
            PT_veg_stp = []
            LAI_veg_stp = []
            PE_stp = []
            E0_stp = []
            RF_veg_stp_tmp = []
            RFe_veg_stp_tmp = []
            PT_veg_stp_tmp = []
            LAI_veg_stp_tmp = []
            PE_stp_tmp=[]
            E0_stp_tmp = []
            if NFIELD <> None:
                RF_irr_stp = []
                RFe_irr_stp = []
                PT_irr_stp = []
                crop_irr_stp = []
                RF_irr_stp_tmp = []
                RFe_irr_stp_tmp = []
                PT_irr_stp_tmp = []
                crop_irr_stp_tmp = []
            self.nper = 0
            self.perlen = []
            perlen_tmp = 0
            c_ = 0
            for n in range(NMETEO):
                RF_veg_stp.append([])
                RF_veg_stp_tmp.append(0.0)
                PT_veg_stp.append([])
                PT_veg_stp_tmp.append([])
                RFe_veg_stp.append([])
                RFe_veg_stp_tmp.append([])
                LAI_veg_stp.append([])
                LAI_veg_stp_tmp.append([])
                for v in range(NVEG):
                    PT_veg_stp[n].append([])
                    PT_veg_stp_tmp[n].append(0.0)
                    RFe_veg_stp[n].append([])
                    RFe_veg_stp_tmp[n].append(0.0)
                    LAI_veg_stp[n].append([])
                    LAI_veg_stp_tmp[n].append(0.0)
                PE_stp.append([])
                PE_stp_tmp.append([])
                for v in range(NSOIL):
                    PE_stp[n].append([])
                    PE_stp_tmp[n].append(0.0)
                E0_stp.append([])
                E0_stp_tmp.append(0.0)
                if NFIELD <> None:
                    RF_irr_stp.append([])
                    RF_irr_stp_tmp.append([])
                    RFe_irr_stp.append([])
                    RFe_irr_stp_tmp.append([])
                    PT_irr_stp.append([])
                    PT_irr_stp_tmp.append([])
                    crop_irr_stp.append([])
                    crop_irr_stp_tmp.append([])
                    for f in range(NFIELD):
                        RF_irr_stp[n].append([])
                        RF_irr_stp_tmp[n].append(0.0)
                        RFe_irr_stp[n].append([])
                        RFe_irr_stp_tmp[n].append(0.0)
                        PT_irr_stp[n].append([])
                        PT_irr_stp_tmp[n].append(0.0)
                        crop_irr_stp[n].append([])
                        crop_irr_stp_tmp[n].append(0.0)
            if perlenmax < 2:
                raise SystemExit('\nFATAL ERROR!\nCorrect perlenmax in the MODFLOW ini file (perlenmax must be higher than 1) or select the daily option.')
            for t in range(len(self.JD)):
                if NFIELD <> None:
                    val_tmp = (RF_veg_d[:,t].sum()+RF_irr_d[:,:,t].sum()).sum()
                else:
                    val_tmp = RF_veg_d[:,t].sum()
                if val_tmp > 0.0:
                    if c_ == 1:
                        for n in range(NMETEO):
                            RF_veg_stp[n].append(RF_veg_stp_tmp[n]/perlen_tmp)
                            RF_veg_stp_tmp[n] = 0.0
                            for v in range(NVEG):
                                PT_veg_stp[n][v].append(PT_veg_stp_tmp[n][v]/perlen_tmp)
                                PT_veg_stp_tmp[n][v] = 0.0
                                RFe_veg_stp[n][v].append(RFe_veg_stp_tmp[n][v]/perlen_tmp)
                                RFe_veg_stp_tmp[n][v] = 0.0
                                LAI_veg_stp[n][v].append(LAI_veg_stp_tmp[n][v]/perlen_tmp)
                                LAI_veg_stp_tmp[n][v] = 0.0
                            for s in range(NSOIL):
                                PE_stp[n][s].append(PE_stp_tmp[n][s]/perlen_tmp)
                                PE_stp_tmp[n][s] = 0.0
                            E0_stp[n].append(E0_stp_tmp[n]/perlen_tmp)
                            E0_stp_tmp[n] = 0.0
                            if NFIELD <> None:
                                for f in range(NFIELD):
                                    RF_irr_stp[n][f].append(RF_irr_stp_tmp[n][f]/perlen_tmp)
                                    RF_irr_stp_tmp[n][f] = 0.0
                                    RFe_irr_stp[n][f].append(RFe_irr_stp_tmp[n][f]/perlen_tmp)
                                    RFe_irr_stp_tmp[n][f] = 0.0
                                    PT_irr_stp[n][f].append(PT_irr_stp_tmp[n][f]/perlen_tmp)
                                    PT_irr_stp_tmp[n][f] = 0.0
                                    crop_irr_stp[n][f].append(np.ceil(crop_irr_stp_tmp[n][f]/perlen_tmp))
                                    crop_irr_stp_tmp[n][f] = 0.0
                        self.perlen.append(perlen_tmp)
                        perlen_tmp = 0
                    for n in range(NMETEO):
                        RF_veg_stp[n].append(RF_veg_d[n,t])
                        for v in range(NVEG):
                            PT_veg_stp[n][v].append(PT_veg_d[n,v,t])
                            RFe_veg_stp[n][v].append(RFe_veg_d[n,v,t])
                            LAI_veg_stp[n][v].append(self.LAI_veg_d[n,v,t])
                        for s in range(NSOIL):
                            PE_stp[n][s].append(PE_d[n,s,t])
                        E0_stp[n].append(E0_d[n,t])
                        if NFIELD <> None:
                            for f in range(NFIELD):
                                RF_irr_stp[n][f].append(RF_irr_d[n,f,t])
                                RFe_irr_stp[n][f].append(RFe_irr_d[n,f,t])
                                PT_irr_stp[n][f].append(PT_irr_d[n,f,t])
                                crop_irr_stp[n][f].append(self.crop_irr_d[n,f,t])
                    self.nper += 1
                    self.perlen.append(1)
                    c_ = 0
                else:
                    if perlen_tmp < perlenmax:
                        for n in range(NMETEO):
                            RF_veg_stp_tmp[n]  += RF_veg_d[n,t]
                            for v in range(NVEG):
                                PT_veg_stp_tmp[n][v] += PT_veg_d[n,v,t]
                                LAI_veg_stp_tmp[n][v] += self.LAI_veg_d[n,v,t]
                            for s in range(NSOIL):
                                PE_stp_tmp[n][s]  += PE_d[n,s,t]
                            E0_stp_tmp[n]  += E0_d[n,t]
                            if NFIELD <> None:
                                for f in range(NFIELD):
                                    RF_irr_stp_tmp[n][f] += RF_irr_d[n,f,t]
                                    PT_irr_stp_tmp[n][f] += PT_irr_d[n,f,t]
                                    crop_irr_stp_tmp[n][f] += self.crop_irr_d[n,f,t]
                        if c_ == 0:
                            self.nper += 1
                        perlen_tmp += 1
                        c_ = 1
                    else:
                        for n in range(NMETEO):
                            RF_veg_stp[n].append(RF_veg_stp_tmp[n]/perlen_tmp)
                            RF_veg_stp_tmp[n] = 0.0
                            for v in range(NVEG):
                                PT_veg_stp[n][v].append(PT_veg_stp_tmp[n][v]/perlen_tmp)
                                PT_veg_stp_tmp[n][v] = 0.0
                                RFe_veg_stp[n][v].append(RFe_veg_stp_tmp[n][v]/perlen_tmp)
                                RFe_veg_stp_tmp[n][v] = 0.0
                                LAI_veg_stp[n][v].append(LAI_veg_stp_tmp[n][v]/perlen_tmp)
                                LAI_veg_stp_tmp[n][v] = 0.0
                            for s in range(NSOIL):
                                PE_stp[n][s].append(PE_stp_tmp[n][s]/perlen_tmp)
                                PE_stp_tmp[n][s] = 0.0
                            E0_stp[n].append(E0_stp_tmp[n]/perlen_tmp)
                            E0_stp_tmp[n] = 0.0
                            if NFIELD <> None:
                                for f in range(NFIELD):
                                    RF_irr_stp[n][f].append(RF_irr_stp_tmp[n][f]/perlen_tmp)
                                    RF_irr_stp_tmp[n][f] = 0.0
                                    RFe_irr_stp[n][f].append(RFe_irr_stp_tmp[n][f]/perlen_tmp)
                                    RFe_irr_stp_tmp[n][f] = 0.0
                                    PT_irr_stp[n][f].append(PT_irr_stp_tmp[n][f]/perlen_tmp)
                                    PT_irr_stp_tmp[n][f] = 0.0
                                    crop_irr_stp[n][f].append(np.ceil(crop_irr_stp_tmp[n][f]/perlen_tmp))
                                    crop_irr_stp_tmp[n][f] = 0.0
                        self.perlen.append(perlen_tmp)
                        for n in range(NMETEO):
                            RF_veg_stp_tmp[n]  += RF_veg_d[n,t]
                            for v in range(NVEG):
                                PT_veg_stp_tmp[n][v] += PT_veg_d[n,v,t]
                                LAI_veg_stp_tmp[n][v] += self.LAI_veg_d[n,v,t]
                            for s in range(NSOIL):
                                PE_stp_tmp[n][s]  += PE_d[n,s,t]
                            E0_stp_tmp[n]  += E0_d[n,t]
                            if NFIELD <> None:
                                for f in range(NFIELD):
                                    RF_irr_stp_tmp[n][f]  += RF_irr_d[n,f,t]
                                    PT_irr_stp_tmp[n][f] += PT_irr_d[n,f,t]
                                    crop_irr_stp_tmp[n][f] += self.crop_irr_d[n,f,t]
                        self.nper += 1
                        perlen_tmp = 1
                        c_ = 1
            if c_ == 1:
                for n in range(NMETEO):
                    RF_veg_stp[n].append(RF_veg_stp_tmp[n]/perlen_tmp)
                    for v in range(NVEG):
                        PT_veg_stp[n][v].append(PT_veg_stp_tmp[n][v]/perlen_tmp)
                        RFe_veg_stp[n][v].append(RFe_veg_stp_tmp[n][v]/perlen_tmp)
                        LAI_veg_stp[n][v].append(LAI_veg_stp_tmp[n][v]/perlen_tmp)
                    for s in range(NSOIL):
                        PE_stp[n][s].append(PE_stp_tmp[n][s]/perlen_tmp)
                    E0_stp[n].append(E0_stp_tmp[n]/perlen_tmp)
                    if NFIELD <> None:
                        for f in range(NFIELD):
                            RF_irr_stp[n][f].append(RF_irr_stp_tmp[n][f]/perlen_tmp)
                            RFe_irr_stp[n][f].append(RFe_irr_stp_tmp[n][f]/perlen_tmp)
                            PT_irr_stp[n][f].append(PT_irr_stp_tmp[n][f]/perlen_tmp)
                            crop_irr_stp[n][f].append(np.ceil(crop_irr_stp_tmp[n][f]/perlen_tmp))
                self.perlen.append(perlen_tmp)
            del perlen_tmp
            del RF_veg_d, RFe_veg_d, PT_veg_d, PE_d, E0_d, RF_veg_stp_tmp, PT_veg_stp_tmp, PE_stp_tmp, E0_stp_tmp, c_
            self.perlen = np.asarray(self.perlen, dtype = np.int)
            self.nstp = np.ones(self.nper, dtype = np.int)
            self.tsmult = self.nstp
            self.Ss_tr = []
            for n in range(self.nper):
                self.Ss_tr.append(False)
        elif self.timedef>0:
            RF_veg_stp = np.zeros([NMETEO,sum(self.nstp)])
            PT_veg_stp = np.zeros([NMETEO,NVEG,sum(self.nstp)])
            RFe_veg_stp = np.zeros([NMETEO,NVEG,sum(self.nstp)])
            LAI_veg_stp = np.zeros([NMETEO,NVEG,sum(self.nstp)])
            PE_stp = np.zeros([NMETEO,NSOIL,sum(self.nstp)])
            E0_stp = np.zeros([NMETEO,sum(self.nstp)])
            if NFIELD <> None:
                RF_irr_stp = np.zeros([NMETEO,NFIELD,(self.nstp)], dtype = float)
                RFe_irr_stp = np.zeros([NMETEO,NFIELD,(self.nstp)], dtype = float)
                PT_irr_stp = np.zeros([NMETEO,NFIELD,(self.nstp)], dtype = float)
                crop_irr_stp = np.zeros([NMETEO,NFIELD,(self.nstp)], dtype = float)
            stp = 0
            for per in range(self.nper):
                for stp in range(self.nstp[per]):
                    tstart = sum(self.perlen[0:per])+stp*(self.perlen[per]/self.nstp[per])
                    tend = sum(self.perlen[0:per])+(self.perlen[per]/self.nstp[per])*(1+stp)
                    for n in range(NMETEO):
                        RF_veg_stp[n,sum(self.nstp[0:per])+stp] += RF_veg_d[n,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                        for v in range(NVEG):
                            PT_veg_stp[n,v,sum(nstp[0:per])+stp] += PT_veg_d[n,v,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                            RFe_veg_stp[n,v,sum(nstp[0:per])+stp] += RFe_veg_d[n,v,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                            LAI_veg_stp[n,v,sum(nstp[0:per])+stp] += self.LAI_veg_d[n,v,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                        for s in range(NSOIL):
                            PE_stp[n,s,sum(nstp[0:per])+stp] += PE_d[n,s,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                        E0_stp[n,sum(nstp[0:per])+stp] += E0_d[n,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                        if NFIELD <> None:
                            for f in range(NFIELD):
                                RF_irr_stp[n,f,sum(self.nstp[0:per])+stp] += RF_irr_d[n,f,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                                RFe_irr_stp[n,f,sum(nstp[0:per])+stp] += RFe_irr_d[n,f,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                                PT_irr_stp[n,f,sum(nstp[0:per])+stp] += PT_irr_d[n,f,tstart:tend].sum()/(self.perlen[per]/self.nstp[per])
                                crop_irr_stp[n,f,sum(nstp[0:per])+stp] += np.ceil(self.crop_irr_d[n,f,tstart:tend].sum()/(self.perlen[per]/self.nstp[per]))

        self.inputZON_SP_RF_veg_fn = "inputZONRF_veg_stp.txt"
        self.inputZONRF_veg_fn = os.path.join(MM_ws, self.inputZON_SP_RF_veg_fn)
        inputZONRF_veg = open(self.inputZONRF_veg_fn, 'w')
        inputZONRF_veg.write('#\n')

        self.inputZON_SP_PT_fn = "inputZONPT_veg_stp.txt"
        self.inputZONPT_fn = os.path.join(MM_ws, self.inputZON_SP_PT_fn)
        inputZONPT = open(self.inputZONPT_fn, 'w')
        inputZONPT.write('#\n')

        self.inputZON_SP_RFe_veg_fn = "inputZONRFe_veg_stp.txt"
        self.inputZONRFe_veg_fn = os.path.join(MM_ws, self.inputZON_SP_RFe_veg_fn)
        inputZONRFe_veg = open(self.inputZONRFe_veg_fn, 'w')
        inputZONRFe_veg.write('#\n')

        self.inputZON_SP_LAI_veg_fn = "inputZONLAI_veg_stp.txt"
        self.inputLAI_veg_fn = os.path.join(MM_ws, self.inputZON_SP_LAI_veg_fn)
        inputLAI_veg = open(self.inputLAI_veg_fn, 'w')
        inputLAI_veg.write('#\n')

        self.inputZON_SP_PE_fn = "inputZONPE_stp.txt"
        self.inputZONPE_fn = os.path.join(MM_ws, self.inputZON_SP_PE_fn)
        inputZONPE = open(self.inputZONPE_fn, 'w')
        inputZONPE.write('#\n')

        self.inputZON_SP_E0_fn = "inputZONE0_stp.txt"
        self.inputZONE0_fn = os.path.join(MM_ws, self.inputZON_SP_E0_fn)
        inputZONE0 = open(self.inputZONE0_fn, 'w')
        inputZONE0.write('#\n')

        if NFIELD <> None:
            self.inputZON_SP_RF_irr_fn = "inputZONRF_irr_stp.txt"
            self.inputZONRF_irr_fn = os.path.join(MM_ws, self.inputZON_SP_RF_irr_fn)
            inputZONRF_irr = open(self.inputZONRF_irr_fn, 'w')
            inputZONRF_irr.write('#\n')

            self.inputZON_SP_RFe_irr_fn = "inputZONRFe_irr_stp.txt"
            self.inputZONRFe_irr_fn = os.path.join(MM_ws, self.inputZON_SP_RFe_irr_fn)
            inputZONRFe_irr = open(self.inputZONRFe_irr_fn, 'w')
            inputZONRFe_irr.write('#\n')

            self.inputZON_SP_PT_irr_fn = "inputZONPT_irr_stp.txt"
            self.inputZONPT_irr_fn = os.path.join(MM_ws, self.inputZON_SP_PT_irr_fn)
            inputZONPT_irr = open(self.inputZONPT_irr_fn, 'w')
            inputZONPT_irr.write('#\n')

            self.input_SP_crop_irr_fn = "inputZONcrop_irr_stp.txt"
            self.inputcrop_irr_fn = os.path.join(MM_ws, self.input_SP_crop_irr_fn)
            inputcrop_irr = open(self.inputcrop_irr_fn, 'w')
            inputcrop_irr.write('#\n')

        try:
            for n in range(NMETEO):
                ExportResults1(RF_veg_stp[n], inputZONRF_veg)
                ExportResults1(E0_stp[n], inputZONE0)
                #VEG
                if NVEG>0:
                    for v in range(NVEG):
                        ExportResults1(PT_veg_stp[n][v], inputZONPT)
                        ExportResults1(RFe_veg_stp[n][v], inputZONRFe_veg)
                # SOIL
                if NSOIL>0:
                    for s in range(NSOIL):
                        ExportResults1(PE_stp[n][s], inputZONPE)
                # CROP
                if NFIELD <> None:
                    for f in range(NFIELD):
                        ExportResults1(RF_irr_stp[n][f], inputZONRF_irr)
                        ExportResults1(RFe_irr_stp[n][f], inputZONRFe_irr)
                        ExportResults1(PT_irr_stp[n][f], inputZONPT_irr)
            if NVEG>0:
                for v in range(NVEG):
                    ExportResults1(LAI_veg_stp[0][v], inputLAI_veg)
            if NFIELD <> None:
                for f in range(NFIELD):
                    ExportResults1(crop_irr_stp[0][f], inputcrop_irr)
        except:
            raise SystemExit("\nFATAL ERROR!\nError in exporting output files in MF time processing.")

        print 'Found %d days converted into %d stress periods.' % (sum(self.perlen), self.nper)

        inputZONRF_veg.close()
        inputZONRFe_veg.close()
        inputLAI_veg.close()
        inputZONPT.close()
        inputZONPE.close()
        inputZONE0.close()
        if NFIELD <> None:
            inputZONRF_irr.close()
            inputZONRFe_irr.close()
            inputZONPT_irr.close()
            inputcrop_irr.close()

#####################################

    def ppMF(self, MF_ws, MF_ini_fn, finf_MM = "", finf_user = None, wel_MM = "", wel_user = None, report = None, verbose = 1, chunks = 0, numDays = -1, MMsoil = 0, silent = 0):

        if verbose == 0:
            print '--------------'

        if type(finf_MM) == np.ndarray:
            finf_input = finf_MM
            wel_input = wel_MM
        else:
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
            if silent == 0: print '\nUZF1 package initialization'
            if isinstance(finf_input,float):
                finf_array = finf_input
                if silent == 0: print 'Infiltration input: %s' % str(finf_input)
            else:
                finf_array = []
                if silent == 0: print 'Infiltration input: %s' % finf_input[0]
                try:
                    h5_finf = h5py.File(finf_input[0], 'r')
                    for n in range(self.nper):
                        finf_array.append(h5_finf[finf_input[1]][n])
                    h5_finf.close()
                except:
                    finf_array = self.finf_user
                    if silent == 0: print 'WARNING!\nNo valid UZF1 package file(s) provided, running MODFLOW using user-input UZF1 infiltration value: %.3G' % self.finf_user
                    finf_input = self.finf_user
            if silent == 0: print "Done!"

        # WELL
        # TODO add well by user to simulate extraction by borehole
        if self.wel_yn == 1:
            if silent == 0: print '\nWEL package initialization'
            if isinstance(wel_input,float):
                wel_array = wel_input
                if silent == 0: print 'Discharge input: %s' % str(wel_input)
            else:
                wel_array = []
                if silent == 0: print 'Discharge input: %s' % wel_input[0]
                try:
                    h5_wel = h5py.File(wel_input[0], 'r')
                    for n in range(self.nper):
                        wel_array.append(h5_wel[wel_input[1]][n])
                    h5_wel.close()
                except:
                    if silent == 0: print 'WARNING!\nNo valid WEL package file(s) provided, running MODFLOW using user-input well value: %.3G' % self.wel_user
                    wel_array = self.wel_user
            # implement a well in every active cell
            layer_row_column_Q = []
            iuzfbnd = np.asarray(self.iuzfbnd)
            wel_dum = 0
            if isinstance(wel_array, float):
                for n in range(self.nper):
                    layer_row_column_Q.append([])
                    for r in range(self.nrow):
                        for c in range(self.ncol):
                            if np.abs(self.ibound[r][c][:]).sum() != 0:
                                if wel_array > 0.0:
                                    layer_row_column_Q[n].append([iuzfbnd[r,c],r+1,c+1,-wel_array*self.delr[c]*self.delc[r]])
                                else:
                                    layer_row_column_Q[n].append([iuzfbnd[r,c],r+1,c+1,0.0])
                                    wel_dum = 1
                            if wel_dum == 1:
                                break
                        if wel_dum == 1:
                            break
                    if wel_dum == 1:
                        break
            else:
                if sum(sum(sum(wel_array))) > 0.0:
                    for n in range(self.nper):
                        layer_row_column_Q.append([])
                        if sum(sum(wel_array[n]))>0.0:
                            for r in range(self.nrow):
                                for c in range(self.ncol):
                                    if np.abs(self.ibound[r][c][:]).sum() != 0:
                                        if wel_array[n][r][c]>0.0:
                                            layer_row_column_Q[n].append([iuzfbnd[r,c],r+1,c+1,-(wel_array[n][r][c])*self.delr[c]*self.delc[r]])
                        else:
                            for r in range(self.nrow):
                                for c in range(self.ncol):
                                    if np.abs(self.ibound[r][c][:]).sum() != 0:
                                        layer_row_column_Q[n].append([iuzfbnd[r,c],r+1,c+1,0.0])
                                        wel_dum = 1
                                    if wel_dum == 1:
                                        break
                                if wel_dum == 1:
                                    break
                            wel_dum = 0
                else:
                    for r in range(self.nrow):
                        for c in range(self.ncol):
                            if np.abs(self.ibound[r][c][:]).sum() != 0:
                                layer_row_column_Q = [[iuzfbnd[r,c],r+1,c+1,0.0]]
                                wel_dum = 1
                            if wel_dum == 1:
                                break
                        if wel_dum == 1:
                            break
            del iuzfbnd
            if silent == 0: print "Done!"

        # RCH
        if self.rch_yn == 1:
            if silent == 0: print '\nRCH package initialization'
            if isinstance(rch_input,float):
                rch_array = rch_input
                if silent == 0: print 'recharge input: %s' % str(rch_input)
            else:
                rch_array = []
                if silent == 0: print 'recharge input: %s' % rch_input[0]
                try:
                    h5_rch = h5py.File(rch_input[0], 'r')
                    for n in range(self.nper):
                        rch_array.append(h5_rch[rch_input[1]][n])
                    h5_rch.close()
                except:
                    rch_array = rch_user
                    if silent == 0: print 'WARNING!\nNo valid RCH package file(s) provided, running MODFLOW using user-input recharge value: %.3G' % rch_user
                    rch_input = rch_user
            if silent == 0: print "Done!"

        # DRAIN
        if self.drn_yn == 1:
            if silent == 0: print '\nDRN package initialization'
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
                                    drn_elev_tmp = self.botm[l] + 0.01
                                else:
                                    drn_elev_tmp = self.botm[i][j][l] + 0.01 #- (botm_array[i][j][l]-botm_array[i][j][l-1])/10
                            else:
                                drn_elev_tmp = drn_elev_array[i,j] + 0.01
                            layer_row_column_elevation_cond[0].append([l+1, i+1, j+1, drn_elev_tmp, drn_cond_array[i,j]])
                            del drn_elev_tmp
                del drn_elev_array, drn_cond_array
                l += 1
            if silent == 0: print "Done!"

        # GHB
        if self.ghb_yn == 1:
            if silent == 0: print '\nGHB package initialization'
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
                            del ghb_head_tmp
                del ghb_head_array, ghb_cond_array
                l += 1
            if silent == 0: print "Done!"

        if type(self.perlen) == np.ndarray: self.perlen = list(self.perlen)
        if type(self.nstp)   == np.ndarray: self.nstp   = list(self.nstp)
        if type(self.tsmult) == np.ndarray: self.tsmult = list(self.tsmult)
        if type(self.Ss_tr)  == np.ndarray: self.Ss_tr  = list(self.Ss_tr)

        # 1 - create the modflow packages files
        if silent == 0: print '\nMODFLOW files writing'
        ncells_package = []
        # MFfile initialization
        mfmain = mf.modflow(modelname = self.modelname, exe_name = self.exe_name, namefile_ext = self.namefile_ext, version = self.version, model_ws = MF_ws, silent = silent)
        # dis package
        dis = mf.mfdis(model = mfmain, nrow = self.nrow, ncol = self.ncol, nlay = self.nlay, nper = self.nper, delr = self.delr, delc = self.delc, laycbd = self.laycbd, top = self.top, botm = self.botm, perlen = self.perlen, nstp = self.nstp, tsmult = self.tsmult, itmuni = self.itmuni, lenuni = self.lenuni, steady = self.Ss_tr, extension = self.ext_dis)
        dis.write_file()
        # bas package
        bas = mf.mfbas(model = mfmain, ibound = self.ibound, strt = self.strt, hnoflo = self.hnoflo, extension = self.ext_bas)
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
        del hk, vka, ss, sy
        # wel package
        if self.wel_yn == 1:
            if layer_row_column_Q <> None:
                wel = mf.mfwel(model = mfmain, iwelcb = cb, layer_row_column_Q = layer_row_column_Q, extension = self.ext_wel)
                wel.write_file()
                del layer_row_column_Q
        # drn package
        if self.drn_yn == 1:
            if drn_check == 1:
                drn = mf.mfdrn(model = mfmain, idrncb = cb, layer_row_column_elevation_cond = layer_row_column_elevation_cond, extension = self.ext_drn)
                drn.write_file()
                del layer_row_column_elevation_cond
        # ghb package
        if self.ghb_yn == 1:
            ghb = mf.mfghb(model = mfmain, ighbcb = cb, layer_row_column_head_cond = layer_row_column_head_cond, extension = self.ext_ghb)
            ghb.write_file()
            del layer_row_column_head_cond
        # uzf package
        if self.uzf_yn == 1:
            uzf = mf.mfuzf1(model = mfmain, nuztop = self.nuztop, iuzfopt = self.iuzfopt, irunflg = self.irunflg, ietflg = self.ietflg, iuzfcb1 = self.iuzfcb1, iuzfcb2 = self.iuzfcb2, ntrail2 = self.ntrail2, nsets = self.nsets, nuzgag = self.nuzgag, surfdep = self.surfdep, iuzfbnd = self.iuzfbnd, vks = vks, eps = eps, thts = thts, thti = thti, row_col_iftunit_iuzopt = self.row_col_iftunit_iuzopt, finf = finf_array, extension = self.ext_uzf, uzfbud_ext = self.uzfbud_ext)
            uzf.write_file()
            del thti, thts, eps, vks, finf_array
        # rch package
        if self.rch_yn == 1:
            rch = mf.mfrch(model = mfmain, irchcb = cb, nrchop = self.nrchop, rech = rch_array, extension = self.ext_rch)
            rch.write_file()
            del rch_array
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
        if silent == 0: print '\nMODFLOW run'
        mfmain.write_name_file()
        mfmain.run_model(pause = False, report = report)
        if silent == 0: print "Done!"

        def readh():
            """
            Extract heads
            """
            try:
                h = mfrdbin.mfhdsread(mfmain, 'LF95').read_all(self.h_MF_fn)
            except:
                h5_MF.close()
                raise ValueError, '\nMODFLOW error!\nCheck the MODFLOW list file in folder:\n%s' % MF_ws
            if type(self.nstp) == np.int32:
                nstp_tmp = self.nstp
            else:
                nstp_tmp = sum(self.nstp)
            if len(h[1]) < nstp_tmp:
                h5_MF.close()
                raise ValueError, '\nMODFLOW error!\nCheck the MODFLOW list file in folder:\n%s' % MF_ws
            if silent == 0: print ''
            return h
            del h

        def readcbc(cbc_fn, h5_ds_fn):
            """
            Extract cell-by-cell budget
            """
            cbc = mfrdbin.mfcbcread(mfmain, 'LF95').read_all(cbc_fn)
            if MMsoil == 1:
                h5_MF.create_dataset(h5_ds_fn, data = np.asarray(cbc[2][0]))
            else:
                h5_MF.create_dataset(h5_ds_fn, data = np.asarray(cbc[2][1]))
            if silent == 0: print ''
            return cbc
            del cbc

        # Create HDF5 arrays to store MF output
        h5_MF = h5py.File(self.h5_MF_fn, 'w')
        if silent == 0: print '\nStoring heads and cbc terms into HDF5 file\n%s\n' % (self.h5_MF_fn)
        if chunks == 1:
            # TODO implement swapaxes
            h = readh()
            h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1]), chunks = (1,self.nrow,self.ncol,self.nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
            del h
            cbc = readcbc(self.cbc_MF_fn, 'cbc_nam')
            h5_MF.create_dataset(name = 'cbc', data = np.asarray(cbc[1]), chunks = (1,len(h5_MF['cbc_nam']),self.nrow,self.ncol,self.nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
            del cbc
            if self.uzf_yn == 1:
                cbc_uzf = readcbc(self.cbc_MFuzf_fn, 'cbc_uzf_nam')
                h5_MF.create_dataset(name = 'cbc_uzf', data = np.asarray(cbc_uzf[1]), chunks = (1,len(h5_MF['cbc_uzf_nam']),self.nrow,self.ncol,self.nlay), compression = 'gzip', compression_opts = 5, shuffle = True)  # 'lzf'
                del cbc_uzf
        else:
            h = readh()
            h5_MF.create_dataset(name = 'heads', data = np.asarray(h[1]))
            del h
            cbc = readcbc(self.cbc_MF_fn, 'cbc_nam')
            data = np.asarray(cbc[1])
            del cbc
            data = np.swapaxes(data,1,2)
            data = np.swapaxes(data,2,3)
            try:
                h5_MF.create_dataset(name = 'cbc', data = data)
            except:
                h5_MF.create_dataset(name = 'cbc', shape = (self.nper, self.nrow, self.ncol, h5_MF['cbc_nam'].shape[0], self.nlay), dtype = np.float)
                for x in range(h5_MF['cbc_nam'].shape[0]):
                    h5_MF['cbc'][:,:,:,x,:] = data[:,:,:,x,:]
            del data
            if self.uzf_yn == 1:
                cbc_uzf = readcbc(self.cbc_MFuzf_fn, 'cbc_uzf_nam')
                data = np.asarray(cbc_uzf[1])
                del cbc_uzf
                data = np.swapaxes(data,1,2)
                data = np.swapaxes(data,2,3)
                try:
                    h5_MF.create_dataset(name = 'cbc_uzf', data = data)
                except:
                    h5_MF.create_dataset(name = 'cbc_uzf', shape = (self.nper, self.nrow, self.ncol, h5_MF['cbc_uzf_nam'].shape[0], self.nlay), dtype = np.float)
                    for x in range(h5_MF['cbc_uzf_nam'].shape[0]):
                        h5_MF['cbc_uzf'][:,:,:,x,:] = data[:,:,:,x,:]
                del data
        h5_MF.close()
        # to delete MF binary files and save disk space
        if self.MFout_yn == 0:
            os.remove(h_MF_fn)
            os.remove(cbc_MF_fn)
            os.remove(cbc_MFuzf_fn)