# -*- coding: utf-8 -*-
"""MODFLOW configuration reader and time-discretization pre-processor.

Phase-1 port (2026): the MODFLOW-NWT model construction and run (runMF)
were removed together with the MM-MF Picard loop (MARMITES_code_review.md,
section 6). This class now only (1) parses the MF ini file, (2) imports
the spatial arrays, and (3) computes the stress-period discretization
(ppMFtime). Phase 3 replaces it by clsMF6 (flopy.mf6 construction +
MODFLOW 6 API coupling driver).
"""

__author__ = "Alain P. Francés <frances.alain@gmail.com>"
__version__ = "0.4.0.dev0"

import os
import sys
import numpy as np
import matplotlib as mpl
import MARMITESprocess_v3 as MMproc

#####################################
class clsMF():
    def __init__(self, cUTIL, MM_ws, MM_ws_out, MF_ws, MF_ini_fn, xllcorner, yllcorner, numDays = -1, stdout = None, report = None):

        ''' read input file (called _input.ini in the MARMITES workspace
        the first character on the first line has to be the character used to comment
        the file can contain any comments as the user wish, but the sequence of the input has to be respected'''

        #global thick_tmp, thick_tmp
        self.cUTIL = cUTIL
        self.MM_ws = MM_ws
        self.MM_ws_out = MM_ws_out
        self.MF_ws = MF_ws
        self.MF_ini_fn = MF_ini_fn
        self.xllcorner = xllcorner
        self.yllcorner = yllcorner
        self.numDays = numDays

        inputFile = self.cUTIL.readFile(self.MF_ws, self.MF_ini_fn)
        l=0
        try:
            self.h5_MF_fn = os.path.join(self.MF_ws, '_h5_MF.h5')
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
            nper = int(inputFile[l].strip())
            if nper<1:
                nper = 1
                print("WARNING! nper must be >0! nper was fixed to 1 (daily time step)")
            self.perlen = []
            self.nstp = []
            self.tsmult = []
            self.Ss_tr = []
            self.nper = nper
            l += 1
            self.itmuni = int(inputFile[l].strip())
            l += 1
            self.lenuni = int(inputFile[l].strip())
            l += 1
            self.laycbd = []
            laycbd_tmp =  inputFile[l].split()
            if len(laycbd_tmp) > 1:
                for i in range(self.nlay):
                    self.laycbd.append(int(laycbd_tmp[i]))
            else:
                for i in range(self.nlay):
                    self.laycbd.append(int(laycbd_tmp[0]))
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
            self.elev = [inputFile[l].strip()]
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
                if self.iupwcb<1:
                    self.iupwcb = 53
                    print("\nWARNING! UZF cbc itunit cannot be <1, fixed to 53!")
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
                print('FATAL ERROR!\nMODFLOW version should be mf2005 or mfnwt!')
                print('Value %s provided in the MF ini file.' % self.version)
                sys.exit()
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
            # OPTIONNAL PACKAGES
            # uzf
            l += 1
            self.uzf_yn = int(inputFile[l].strip())
            if self.uzf_yn == 1 :
                l += 1
                self.ext_uzf = str(inputFile[l].strip())
                l += 1
                specifithtr_tmp =  inputFile[l].split()
                self.specifythtr = int(specifithtr_tmp[0].strip())
                self.thtr = [specifithtr_tmp[1].strip()]
                l += 1
                self.specifythti = int(inputFile[l].strip())
                l += 1
                self.nosurfleak = int(inputFile[l].strip())
                l += 1
                self.nuztop = int(inputFile[l].strip())
                l += 1
                self.iuzfopt = int(inputFile[l].strip())
                l += 1
                self.irunflg = int(inputFile[l].strip())
                l += 1
                self.ietflg = int(inputFile[l].strip())
                if self.ietflg>0:
                    l += 1
                    #TO DO: read ET option in MM
                else:
                    l += 6
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
                if self.iuzfopt ==1:
                    self.vks = [inputFile[l].strip()]
                else:
                    self.vks = None
                l += 1
                self.eps = [inputFile[l].strip()]
                l += 1
                self.thts = [inputFile[l].strip()]
                l += 1
                self.thti = [inputFile[l].strip()]
                if self.nuzgag > 0:
                    l += 1
                    self.iuzrow =  inputFile[l].split()
                    l += 1
                    self.iuzcol = inputFile[l].split()
                    l += 1
                    self.iftunit = inputFile[l].split()
                    l += 1
                    self.iuzopt = inputFile[l].split()
                    l += 1
                self.perc_user = float(inputFile[l].strip())
            else:
                l += 24
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
            # NOT MODFLOW - define main, hydrogeological layers from MODFLOW layers
            l += 1
            self.Mnlay = int(inputFile[l].strip())
            l += 1
            self.Mlay = []
            Mlay_tmp = inputFile[l].split()
            for i in range(self.nlay):
                self.Mlay.append(int(Mlay_tmp[i]))
            if max(self.Mlay) != self.Mnlay:
                print('\nERROR!\nMnlay must be a list with nlay values that must be incremental and with maximum igual to Mnlay')
                sys.exit()
            # NOT MODFLOW - define wich layers are plotted on the output time series
            l += 1
            self.h_plt = []
            h_plt_tmp = inputFile[l].split()
            for i in range(self.nlay):
                self.h_plt.append(int(h_plt_tmp[i]))
            l += 1
            self.h_lbl = []
            h_lbl_tmp = inputFile[l].split()
            for i in range(self.nlay):
                self.h_lbl.append(h_lbl_tmp[i])
        except Exception:
            self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nUnexpected error in the MODFLOW input file:\n", stdout = stdout, report = report)
        del inputFile

        self.cPROCESS = MMproc.clsPROCESS(
                                cUTIL           = self.cUTIL,
                                MM_ws           = self.MM_ws,
                                MM_ws_out       = self.MM_ws_out,
                                MF_ws           = self.MF_ws,
                                nrow            = self.nrow,
                                ncol            = self.ncol,
                                nlay            = self.nlay,
                                xllcorner       = self.xllcorner,
                                yllcorner       = self.yllcorner,
                                cellsizeMF      = self.delr[0],
                                hnoflo          = self.hnoflo
                                )

        # 1 - reaf asc file and convert in np.array
        print("\nImporting ESRI ASCII files to initialize the MODFLOW packages...")

        self.elev     = self.cPROCESS.checkarray(self.elev, stdout = stdout, report = report)
        self.strt     = self.cPROCESS.checkarray(self.strt, stdout = stdout, report = report)
        self.thick    = self.cPROCESS.checkarray(self.thick, stdout = stdout, report = report)
        self.ibound   = self.cPROCESS.checkarray(self.ibound, dtype = int, stdout = stdout, report = report)
        if self.nlay < 2:
            if isinstance(self.ibound, (list, np.ndarray)):
                self.ibound = (np.asarray(self.ibound)).reshape((1, self.nrow, self.ncol))
            if isinstance(self.strt, (list, np.ndarray)):
                self.strt = (np.asarray(self.strt)).reshape((1, self.nrow, self.ncol))                
        else:
            self.ibound = np.asarray(self.ibound)
            self.strt = np.asarray(self.strt)
        self.thick = self.cPROCESS.float2array(self.thick)
        elev_tmp = np.ma.masked_values(np.asarray(self.elev), self.hnoflo, atol = 0.09)
        botm_tmp = []
        for l in range(self.nlay):
            if l == 0:
                thick_tmp = self.thick[l,:,:]*np.abs(self.ibound)[l,:,:]
            else:
                thick_tmp += self.thick[l,:,:]*np.abs(self.ibound)[l,:,:]
            botm_tmp.append(np.ma.masked_values(elev_tmp - thick_tmp, self.hnoflo, atol = 0.09))
        del elev_tmp, thick_tmp
#        self.botm = list(np.swapaxes(botm_tmp,0,1))
#        self.botm = list(np.swapaxes(self.botm,1,2))
        self.botm = np.asarray(botm_tmp)
        del botm_tmp
        if self.nlay < 2:
            if isinstance(self.botm, list):
                self.botm = np.ma.masked_values((np.asarray(self.botm)).reshape((1, self.nrow, self.ncol)), self.hnoflo, atol = 0.09)
        if self.uzf_yn == 1:
            self.iuzfbnd = self.cPROCESS.checkarray(self.iuzfbnd, dtype = int, stdout = stdout, report = report)
        if self.ghb_yn == 1:
            ghb = np.asarray(self.cPROCESS.checkarray(self.ghb_cond, dtype = np.float32, stdout = stdout, report = report))
            if self.nlay > 1:
                self.ghbcells = np.zeros((self.nlay), dtype = int)
                for l in range(self.nlay):
                    self.ghbcells[l] = (np.asarray(ghb[l,:,:]) > 0.0).sum()
            else:
                self.ghbcells = [(np.asarray(ghb[:,:]) > 0.0).sum()]
            del ghb
        if self.drn_yn == 1:
            drn = np.asarray(self.cPROCESS.checkarray(self.drn_cond, dtype = np.float32, stdout = stdout, report = report))
            if self.nlay > 1:
                self.drncells = np.zeros((self.nlay), dtype = int)
                for l in range(self.nlay):
                    self.drncells[l] = (np.asarray(drn[l,:,:]) > 0.0).sum()
            else:
                self.drncells = [(np.asarray(drn[:,:]) > 0.0).sum()]
            del drn

        self.array_ini(MF_ws = self.MF_ws, stdout = stdout, report = report)

####################################

    def array_ini(self, MF_ws, stdout = None, report = None):

        self.hk_actual      = self.cPROCESS.checkarray(self.hk, stdout = stdout, report = report)
        self.vka_actual     = self.cPROCESS.checkarray(self.vka, stdout = stdout, report = report)
        self.ss_actual      = self.cPROCESS.checkarray(self.ss, stdout = stdout, report = report)
        self.sy_actual      = self.cPROCESS.checkarray(self.sy, stdout = stdout, report = report)
        if np.sum(np.asarray(self.sy_actual)>1.0)>0.0:
            self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nSy value > 1.0!!!\nCorrect it (it should be expressed in [L^3/L^3]), and verify also thts, thtr and thti.", stdout = stdout, report = report)
        # uzf        
        if self.uzf_yn == 1:
            sy_actual = self.cPROCESS.float2array(self.sy_actual)
            if self.iuzfopt == 1:
                self.vks_actual     = self.cPROCESS.checkarray(self.vks, stdout = stdout, report = report)
            else:
                self.vks_actual = None
            self.eps_actual     = self.cPROCESS.checkarray(self.eps, stdout = stdout, report = report)
            self.thtr_actual = None
            if self.specifythtr > 0:
                self.thtr_actual     = self.cPROCESS.checkarray(self.thtr, stdout = stdout, report = report)
                thtr_tmp  = np.asarray(self.thtr_actual)
            else:
                self.thtr_actual = 0.0
                thtr_tmp = 0.0
            try:
                self.thts_actual = float(self.thts[0])
            except Exception:
                self.thts_actual = self.thts[0]
            if isinstance(self.thts_actual, float) and self.thts_actual < 0:
                self.thts_actual = np.abs(self.thts_actual)
                for l in range(self.nlay):
                    if l == 0:
                        if len(sy_actual) > self.nlay:
                            self.thts_actual += sy_actual[l,:,:]*np.abs(self.ibound)[l,:,:] + thtr_tmp
                        else:
                            self.thts_actual += sy_actual[l]*np.abs(self.ibound)[l,:,:] + thtr_tmp
                    else:
                        if len(sy_actual) > self.nlay:
                            self.thts_actual += sy_actual[l,:,:]*np.abs(self.ibound)[l,:,:]*abs(np.abs(self.ibound)[l-1,:,:]-1)
                        else:
                            self.thts_actual += sy_actual[l]*np.abs(self.ibound)[l,:,:]*abs(np.abs(self.ibound)[l-1,:,:]-1)
                if self.thtr_actual is not None:
                    del thtr_tmp
            else:
                self.thts_actual = self.cPROCESS.checkarray(self.thts, stdout = stdout, report = report)
            try:
                self.thti_actual = float(self.thti[0])
            except Exception:
                self.thti_actual = self.thti[0]
            if isinstance(self.thti_actual, float) and self.thti_actual < 0:
                self.thti_actual = self.thts_actual/(np.abs(self.thti_actual))
            else:
                self.thti_actual = self.cPROCESS.checkarray(self.thti, stdout = stdout, report = report)
            for l in range(self.nlay):
                try:
                    sy_tmp = np.asarray(sy_actual[l,:,:])
                except Exception:
                     sy_tmp = np.asarray(sy_actual[l])
                if (sy_tmp*np.abs(self.ibound)[l,:,:]+np.asarray(self.thtr_actual)>np.asarray(self.thts_actual)).sum()> 0.0:
                    self.thts_actual = sy_tmp*np.abs(self.ibound)[l,:,:]+2.0*np.asarray(self.thtr_actual)
                    print('\nWARNING!\nSy + THTR > THTS in layer %d! Corrected: THTS = Sy + 2.0*THTR' % l)
            if (np.asarray(self.thti_actual)<np.asarray(self.thtr_actual)).sum()>0.0 or (np.asarray(self.thti_actual)>np.asarray(self.thts_actual)).sum()>0.0:
                self.thti_actual = np.asarray(self.thtr_actual) + (np.asarray(self.thts_actual)-np.asarray(self.thtr_actual))/4.0
                print('\nWARNING!\nTHTI > THTS or THTI< THTR!Corrected: THTI = THTR + (THTS-THTR)/4.0')
            del sy_actual
            
        # DRAIN
        if self.drn_yn == 1:
            print('\nDRN package initialization')
            l = 0
            self.drn_elev_array = np.zeros((self.nlay,self.nrow,self.ncol))
            self.drn_cond_array = np.zeros((self.nlay,self.nrow,self.ncol))
            self.layer_row_column_elevation_cond = {}
            self.layer_row_column_elevation_cond[0] = []
            for d in self.drn_cond:
                if isinstance(d, str):
                    drn_elev_path = os.path.join(self.MF_ws, self.drn_elev[l])
                    self.drn_elev_array[l,:,:] = self.cPROCESS.convASCIIraster2array(drn_elev_path, self.drn_elev_array[l,:,:])
                    drn_cond_path = os.path.join(self.MF_ws, self.drn_cond[l])
                    self.drn_cond_array[l,:,:] = self.cPROCESS.convASCIIraster2array(drn_cond_path, self.drn_cond_array[l,:,:])
                else:
                    self.drn_elev_array[l,:,:] = self.drn_elev[l]
                    self.drn_cond_array[l,:,:] = self.drn_cond[l]
                for i in range(self.nrow):
                    for j in range(self.ncol):
                        if self.drn_elev_array[l,i,j]!=0:
                            if self.drn_elev_array[l,i,j]<0 and self.drn_elev_array[l,i,j] != self.hnoflo:
                                if isinstance(self.botm[l], float):
                                    drn_elev_tmp = self.botm[l] + 0.01
                                else:
                                    drn_elev_tmp = self.botm[l,i,j] + 0.01 #- (botm_array[l][i][j]-botm_array[i][j][l-1])/10
                            else:
                                drn_elev_tmp = self.drn_elev_array[l,i,j] + 0.01
                            if self.drn_cond_array[l,i,j] != self.hnoflo:
                                self.layer_row_column_elevation_cond[0].append([l, i, j, drn_elev_tmp, self.drn_cond_array[l,i,j]])
                            self.drn_elev_array[l,i,j] = drn_elev_tmp
                            del drn_elev_tmp
                l += 1
            print("Done!")

        # GHB
        if self.ghb_yn == 1:
            print('\nGHB package initialization')
            l = 0
            self.ghb_head_array = np.zeros((self.nlay,self.nrow,self.ncol))
            self.ghb_cond_array = np.zeros((self.nlay,self.nrow,self.ncol))
            self.layer_row_column_head_cond = {}
            self.layer_row_column_head_cond[0] = []
            for d in self.ghb_cond:
                if isinstance(d, str):
                    ghb_head_path = os.path.join(self.MF_ws, self.ghb_head[l])
                    self.ghb_head_array[l,:,:] = self.cPROCESS.convASCIIraster2array(ghb_head_path, self  .ghb_head_array[l,:,:])
                    ghb_cond_path = os.path.join(self.MF_ws, self.ghb_cond[l])
                    self.ghb_cond_array[l,:,:] = self.cPROCESS.convASCIIraster2array(ghb_cond_path, self.ghb_cond_array[l,:,:])
                else:
                    self.ghb_head_array[l,:,:] = self.ghb_head[l]
                    self.ghb_cond_array[l,:,:] = self.ghb_cond[l]
                for i in range(self.nrow):
                    for j in range(self.ncol):
                        if self.ghb_head_array[l,i,j] != 0 and self.ghb_head_array[l,i,j] != self.hnoflo:
                            ghb_head_tmp = self.ghb_head_array[l,i,j]
                            if self.ghb_cond_array[l,i,j] != self.hnoflo:
                                self.layer_row_column_head_cond[0].append([l, i, j, ghb_head_tmp, self.ghb_cond_array[l,i,j]])
                            del ghb_head_tmp
                l += 1
            print("Done!")

####################################


    def ppMFtime(self, inputDate_fn, inputZON_dSP_P_veg_fn, inputZON_dSP_Pe_veg_fn, inputZON_dSP_PT_fn, input_dSP_LAI_veg_fn, inputZON_dSP_PE_fn, inputZON_dSP_Eo_fn, NMETEO, NVEG, NSOIL, inputZON_dSP_P_irr_fn = None, inputZON_dSP_Pe_irr_fn = None, inputZON_dSP_PT_irr_fn = None, input_dSP_crop_irr_fn = None, NFIELD = None, stdout = None, report = None):

        ''' P analysis to define automatically nper/perlen/nstp
        Daily P>0 creates a nper
        Succesive days with P=0 are averaged (PT, PE, Eo) to a user-input maximum # of days'''

        #global inputcrop_irr, inputZON_PT_irr, inputZON_Pe_irr, inputZON_P_irr, inputFilecrop_irr, inputFilePT_irr, PT_irr_d, inputFilePe_irr, Pe_irr_d, inputFileP_irr, P_irr_d

        def ExportResults1(SP, outFileExport):
            """
            Write the processed data in a open txt file readable by  MARMITES
            INPUT:      output flux time series and open file
            """
            for sp in range(len(SP)):
                out_line =  '%14.9G' %SP[sp],'\n'
                for l in out_line:
                    outFileExport.write(l)

        #####################################

        # READ date of input files (P and PT)
        inputDate_fn=os.path.join(self.MM_ws, inputDate_fn)
        if os.path.exists(inputDate_fn):
            inputDate_tmp = np.loadtxt(inputDate_fn, dtype = str)
            self.inputDate = inputDate_tmp[:,0]
            self.JD = np.asarray(inputDate_tmp[:,2], dtype = int)
            del inputDate_tmp
            self.inputDate = mpl.dates.datestr2num(self.inputDate)
            for t in range(1,len(self.inputDate)):
                #__________________Check date consistency________________#
                difDay=self.inputDate[t]-self.inputDate[t-1]
                if (difDay !=1.0):
                    print('DifDay = ' + str(difDay))
                    self.cUTIL.ErrorExit(msg = '\nFATAL ERROR!\nDates of the input data are not sequencial, check your daily time step!\nError at date %s ' % str(self.inputDate[t]), stdout = stdout, report = report)
        else:
            self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe file %s doesn't exist!!!" % inputDate_fn, stdout = stdout, report = report)

        inputFileP_veg = self.cUTIL.readFile(self.MM_ws, inputZON_dSP_P_veg_fn)
        P_veg_d = np.zeros([NMETEO, len(self.JD)], dtype = np.float32)
        inputFilePe_veg = self.cUTIL.readFile(self.MM_ws, inputZON_dSP_Pe_veg_fn)
        Pe_veg_d = np.zeros([NMETEO, NVEG, len(self.JD)], dtype = np.float32)
        inputFilePT = self.cUTIL.readFile(self.MM_ws, inputZON_dSP_PT_fn)
        PT_veg_d = np.zeros([NMETEO, NVEG, len(self.JD)], dtype = np.float32)
        inputFileLAI_veg = self.cUTIL.readFile(self.MM_ws, input_dSP_LAI_veg_fn)
        self.LAI_veg_d = np.zeros([NMETEO, NVEG, len(self.JD)], dtype = np.float32)
        inputFilePE = self.cUTIL.readFile(self.MM_ws, inputZON_dSP_PE_fn)
        PE_d = np.zeros([NMETEO, NSOIL, len(self.JD)], dtype = np.float32)
        inputFileEo = self.cUTIL.readFile(self.MM_ws, inputZON_dSP_Eo_fn)
        Eo_d = np.zeros([NMETEO, len(self.JD)], dtype = np.float32)
        if NFIELD is not None:
            inputFileP_irr = self.cUTIL.readFile(self.MM_ws, inputZON_dSP_P_irr_fn)
            P_irr_d = np.zeros([NMETEO, NFIELD, len(self.JD)], dtype = np.float32)
            inputFilePe_irr = self.cUTIL.readFile(self.MM_ws, inputZON_dSP_Pe_irr_fn)
            Pe_irr_d = np.zeros([NMETEO, NFIELD, len(self.JD)], dtype = np.float32)
            inputFilePT_irr = self.cUTIL.readFile(self.MM_ws, inputZON_dSP_PT_irr_fn)
            PT_irr_d = np.zeros([NMETEO, NFIELD, len(self.JD)], dtype = np.float32)
            inputFilecrop_irr = self.cUTIL.readFile(self.MM_ws, input_dSP_crop_irr_fn)
            self.crop_irr_d = np.zeros([NMETEO, NFIELD, len(self.JD)], dtype = np.float32)
        for n in range(NMETEO):
            for t in range(len(self.JD)):
                P_veg_d[n,t] = float(inputFileP_veg[t+len(self.JD)*n].strip())
                Eo_d[n,t] = float(inputFileEo[t+len(self.JD)*n].strip())
                for v in range(NVEG):
                    PT_veg_d[n,v,t] = float(inputFilePT[t+(n*NVEG+v)*len(self.JD)].strip())
                    Pe_veg_d[n,v,t] = float(inputFilePe_veg[t+(n*NVEG+v)*len(self.JD)].strip())
                    self.LAI_veg_d[n,v,t] = float(inputFileLAI_veg[t+v*len(self.JD)].strip())
                for s in range(NSOIL):
                    PE_d[n,s,t] = float(inputFilePE[t+(n*NSOIL+s)*len(self.JD)].strip())
                if NFIELD is not None:
                    for f in range(NFIELD):
                        P_irr_d[n,f,t] = float(inputFileP_irr[t+(n*NFIELD+f)*len(self.JD)].strip())
                        Pe_irr_d[n,f,t] = float(inputFilePe_irr[t+(n*NFIELD+f)*len(self.JD)].strip())
                        PT_irr_d[n,f,t] = float(inputFilePT_irr[t+(n*NFIELD+f)*len(self.JD)].strip())
                        self.crop_irr_d[n,f,t] = float(inputFilecrop_irr[t+f*len(self.JD)].strip())
        del inputFileP_veg, inputFilePT, inputFilePE, inputFilePe_veg, inputFileEo
        if NFIELD is not None:
            del inputZON_dSP_P_irr_fn, inputZON_dSP_Pe_irr_fn, inputZON_dSP_PT_irr_fn, input_dSP_crop_irr_fn, inputFileP_irr, inputFilePe_irr, inputFilePT_irr, inputFilecrop_irr

        print('\nComputing MODFLOW time discretization based on rainfall analysis in each METEO zone...')
        # summing P and avering other flux when P = 0
        perlenmax = self.nper
        P_veg_stp = []
        Pe_veg_stp = []
        PT_veg_stp = []
        LAI_veg_stp = []
        PE_stp = []
        Eo_stp = []
        P_veg_stp_tmp = []
        Pe_veg_stp_tmp = []
        PT_veg_stp_tmp = []
        LAI_veg_stp_tmp = []
        PE_stp_tmp=[]
        Eo_stp_tmp = []
        if NFIELD is not None:
            P_irr_stp = []
            Pe_irr_stp = []
            PT_irr_stp = []
            crop_irr_stp = []
            P_irr_stp_tmp = []
            Pe_irr_stp_tmp = []
            PT_irr_stp_tmp = []
            crop_irr_stp_tmp = []
        self.nper = 0
        self.perlen = []
        perlen_tmp = 0
        c_ = 0
        for n in range(NMETEO):
            P_veg_stp.append([])
            P_veg_stp_tmp.append(0.0)
            PT_veg_stp.append([])
            PT_veg_stp_tmp.append([])
            Pe_veg_stp.append([])
            Pe_veg_stp_tmp.append([])
            LAI_veg_stp.append([])
            LAI_veg_stp_tmp.append([])
            for v in range(NVEG):
                PT_veg_stp[n].append([])
                PT_veg_stp_tmp[n].append(0.0)
                Pe_veg_stp[n].append([])
                Pe_veg_stp_tmp[n].append(0.0)
                LAI_veg_stp[n].append([])
                LAI_veg_stp_tmp[n].append(0.0)
            PE_stp.append([])
            PE_stp_tmp.append([])
            for v in range(NSOIL):
                PE_stp[n].append([])
                PE_stp_tmp[n].append(0.0)
            Eo_stp.append([])
            Eo_stp_tmp.append(0.0)
            if NFIELD is not None:
                P_irr_stp.append([])
                P_irr_stp_tmp.append([])
                Pe_irr_stp.append([])
                Pe_irr_stp_tmp.append([])
                PT_irr_stp.append([])
                PT_irr_stp_tmp.append([])
                crop_irr_stp.append([])
                crop_irr_stp_tmp.append([])
                for f in range(NFIELD):
                    P_irr_stp[n].append([])
                    P_irr_stp_tmp[n].append(0.0)
                    Pe_irr_stp[n].append([])
                    Pe_irr_stp_tmp[n].append(0.0)
                    PT_irr_stp[n].append([])
                    PT_irr_stp_tmp[n].append(0.0)
                    crop_irr_stp[n].append([])
                    crop_irr_stp_tmp[n].append(0.0)
        for t in range(len(self.JD)):
            if NFIELD is not None:
                val_tmp = (P_veg_d[:,t].sum()+P_irr_d[:,:,t].sum()).sum()
            else:
                val_tmp = P_veg_d[:,t].sum()
            if val_tmp > 0.0:
                if c_ == 1:
                    for n in range(NMETEO):
                        P_veg_stp[n].append(P_veg_stp_tmp[n]/perlen_tmp)
                        P_veg_stp_tmp[n] = 0.0
                        for v in range(NVEG):
                            PT_veg_stp[n][v].append(PT_veg_stp_tmp[n][v]/perlen_tmp)
                            PT_veg_stp_tmp[n][v] = 0.0
                            Pe_veg_stp[n][v].append(Pe_veg_stp_tmp[n][v]/perlen_tmp)
                            Pe_veg_stp_tmp[n][v] = 0.0
                            LAI_veg_stp[n][v].append(LAI_veg_stp_tmp[n][v]/perlen_tmp)
                            LAI_veg_stp_tmp[n][v] = 0.0
                        for s in range(NSOIL):
                            PE_stp[n][s].append(PE_stp_tmp[n][s]/perlen_tmp)
                            PE_stp_tmp[n][s] = 0.0
                        Eo_stp[n].append(Eo_stp_tmp[n]/perlen_tmp)
                        Eo_stp_tmp[n] = 0.0
                        if NFIELD is not None:
                            for f in range(NFIELD):
                                P_irr_stp[n][f].append(P_irr_stp_tmp[n][f]/perlen_tmp)
                                P_irr_stp_tmp[n][f] = 0.0
                                Pe_irr_stp[n][f].append(Pe_irr_stp_tmp[n][f]/perlen_tmp)
                                Pe_irr_stp_tmp[n][f] = 0.0
                                PT_irr_stp[n][f].append(PT_irr_stp_tmp[n][f]/perlen_tmp)
                                PT_irr_stp_tmp[n][f] = 0.0
                                crop_irr_stp[n][f].append(np.ceil(crop_irr_stp_tmp[n][f]/perlen_tmp))
                                crop_irr_stp_tmp[n][f] = 0.0
                    self.perlen.append(perlen_tmp)
                    perlen_tmp = 0
                for n in range(NMETEO):
                    P_veg_stp[n].append(P_veg_d[n,t])
                    for v in range(NVEG):
                        PT_veg_stp[n][v].append(PT_veg_d[n,v,t])
                        Pe_veg_stp[n][v].append(Pe_veg_d[n,v,t])
                        LAI_veg_stp[n][v].append(self.LAI_veg_d[n,v,t])
                    for s in range(NSOIL):
                        PE_stp[n][s].append(PE_d[n,s,t])
                    Eo_stp[n].append(Eo_d[n,t])
                    if NFIELD is not None:
                        for f in range(NFIELD):
                            P_irr_stp[n][f].append(P_irr_d[n,f,t])
                            Pe_irr_stp[n][f].append(Pe_irr_d[n,f,t])
                            PT_irr_stp[n][f].append(PT_irr_d[n,f,t])
                            crop_irr_stp[n][f].append(self.crop_irr_d[n,f,t])
                self.nper += 1
                self.perlen.append(1)
                c_ = 0
            else:
                if perlen_tmp < perlenmax:
                    for n in range(NMETEO):
                        P_veg_stp_tmp[n]  += P_veg_d[n,t]
                        for v in range(NVEG):
                            PT_veg_stp_tmp[n][v] += PT_veg_d[n,v,t]
                            LAI_veg_stp_tmp[n][v] += self.LAI_veg_d[n,v,t]
                        for s in range(NSOIL):
                            PE_stp_tmp[n][s]  += PE_d[n,s,t]
                        Eo_stp_tmp[n]  += Eo_d[n,t]
                        if NFIELD is not None:
                            for f in range(NFIELD):
                                P_irr_stp_tmp[n][f] += P_irr_d[n,f,t]
                                PT_irr_stp_tmp[n][f] += PT_irr_d[n,f,t]
                                crop_irr_stp_tmp[n][f] += self.crop_irr_d[n,f,t]
                    if c_ == 0:
                        self.nper += 1
                    perlen_tmp += 1
                    c_ = 1
                else:
                    for n in range(NMETEO):
                        P_veg_stp[n].append(P_veg_stp_tmp[n]/perlen_tmp)
                        P_veg_stp_tmp[n] = 0.0
                        for v in range(NVEG):
                            PT_veg_stp[n][v].append(PT_veg_stp_tmp[n][v]/perlen_tmp)
                            PT_veg_stp_tmp[n][v] = 0.0
                            Pe_veg_stp[n][v].append(Pe_veg_stp_tmp[n][v]/perlen_tmp)
                            Pe_veg_stp_tmp[n][v] = 0.0
                            LAI_veg_stp[n][v].append(LAI_veg_stp_tmp[n][v]/perlen_tmp)
                            LAI_veg_stp_tmp[n][v] = 0.0
                        for s in range(NSOIL):
                            PE_stp[n][s].append(PE_stp_tmp[n][s]/perlen_tmp)
                            PE_stp_tmp[n][s] = 0.0
                        Eo_stp[n].append(Eo_stp_tmp[n]/perlen_tmp)
                        Eo_stp_tmp[n] = 0.0
                        if NFIELD is not None:
                            for f in range(NFIELD):
                                P_irr_stp[n][f].append(P_irr_stp_tmp[n][f]/perlen_tmp)
                                P_irr_stp_tmp[n][f] = 0.0
                                Pe_irr_stp[n][f].append(Pe_irr_stp_tmp[n][f]/perlen_tmp)
                                Pe_irr_stp_tmp[n][f] = 0.0
                                PT_irr_stp[n][f].append(PT_irr_stp_tmp[n][f]/perlen_tmp)
                                PT_irr_stp_tmp[n][f] = 0.0
                                crop_irr_stp[n][f].append(np.ceil(crop_irr_stp_tmp[n][f]/perlen_tmp))
                                crop_irr_stp_tmp[n][f] = 0.0
                    self.perlen.append(perlen_tmp)
                    for n in range(NMETEO):
                        P_veg_stp_tmp[n]  += P_veg_d[n,t]
                        for v in range(NVEG):
                            PT_veg_stp_tmp[n][v] += PT_veg_d[n,v,t]
                            LAI_veg_stp_tmp[n][v] += self.LAI_veg_d[n,v,t]
                        for s in range(NSOIL):
                            PE_stp_tmp[n][s]  += PE_d[n,s,t]
                        Eo_stp_tmp[n]  += Eo_d[n,t]
                        if NFIELD is not None:
                            for f in range(NFIELD):
                                P_irr_stp_tmp[n][f]  += P_irr_d[n,f,t]
                                PT_irr_stp_tmp[n][f] += PT_irr_d[n,f,t]
                                crop_irr_stp_tmp[n][f] += self.crop_irr_d[n,f,t]
                    self.nper += 1
                    perlen_tmp = 1
                    c_ = 1
        if c_ == 1:
            for n in range(NMETEO):
                P_veg_stp[n].append(P_veg_stp_tmp[n]/perlen_tmp)
                for v in range(NVEG):
                    PT_veg_stp[n][v].append(PT_veg_stp_tmp[n][v]/perlen_tmp)
                    Pe_veg_stp[n][v].append(Pe_veg_stp_tmp[n][v]/perlen_tmp)
                    LAI_veg_stp[n][v].append(LAI_veg_stp_tmp[n][v]/perlen_tmp)
                for s in range(NSOIL):
                    PE_stp[n][s].append(PE_stp_tmp[n][s]/perlen_tmp)
                Eo_stp[n].append(Eo_stp_tmp[n]/perlen_tmp)
                if NFIELD is not None:
                    for f in range(NFIELD):
                        P_irr_stp[n][f].append(P_irr_stp_tmp[n][f]/perlen_tmp)
                        Pe_irr_stp[n][f].append(Pe_irr_stp_tmp[n][f]/perlen_tmp)
                        PT_irr_stp[n][f].append(PT_irr_stp_tmp[n][f]/perlen_tmp)
                        crop_irr_stp[n][f].append(np.ceil(crop_irr_stp_tmp[n][f]/perlen_tmp))
            self.perlen.append(perlen_tmp)
        del perlen_tmp
        del P_veg_d, Pe_veg_d, PT_veg_d, PE_d, Eo_d, P_veg_stp_tmp, PT_veg_stp_tmp, PE_stp_tmp, Eo_stp_tmp, c_
        self.perlen = np.asarray(self.perlen, dtype = int)
        self.nstp = np.ones(self.nper, dtype = int)
        self.tsmult = self.nstp
        self.Ss_tr = []
        for n in range(self.nper):
            self.Ss_tr.append(False)

        self.inputZON_SP_P_veg_fn = "inputZON_P_veg_stp.txt"
        self.inputZON_P_veg_fn = os.path.join(self.MM_ws, self.inputZON_SP_P_veg_fn)
        inputZON_P_veg = open(self.inputZON_P_veg_fn, 'w')
        inputZON_P_veg.write('#\n')

        self.inputZON_SP_PT_fn = "inputZON_PT_veg_stp.txt"
        self.inputZON_PT_fn = os.path.join(self.MM_ws, self.inputZON_SP_PT_fn)
        inputZON_PT = open(self.inputZON_PT_fn, 'w')
        inputZON_PT.write('#\n')

        self.inputZON_SP_Pe_veg_fn = "inputZON_Pe_veg_stp.txt"
        self.inputZON_Pe_veg_fn = os.path.join(self.MM_ws, self.inputZON_SP_Pe_veg_fn)
        inputZON_Pe_veg = open(self.inputZON_Pe_veg_fn, 'w')
        inputZON_Pe_veg.write('#\n')

        self.inputZON_SP_LAI_veg_fn = "inputZON_LAI_veg_stp.txt"
        self.inputLAI_veg_fn = os.path.join(self.MM_ws, self.inputZON_SP_LAI_veg_fn)
        inputLAI_veg = open(self.inputLAI_veg_fn, 'w')
        inputLAI_veg.write('#\n')

        self.inputZON_SP_PE_fn = "inputZON_PE_stp.txt"
        self.inputZON_PE_fn = os.path.join(self.MM_ws, self.inputZON_SP_PE_fn)
        inputZON_PE = open(self.inputZON_PE_fn, 'w')
        inputZON_PE.write('#\n')

        self.inputZON_SP_Eo_fn = "inputZON_Eo_stp.txt"
        self.inputZON_Eo_fn = os.path.join(self.MM_ws, self.inputZON_SP_Eo_fn)
        inputZON_Eo = open(self.inputZON_Eo_fn, 'w')
        inputZON_Eo.write('#\n')

        if NFIELD is not None:
            self.inputZON_SP_P_irr_fn = "inputZON_P_irr_stp.txt"
            self.inputZON_P_irr_fn = os.path.join(self.MM_ws, self.inputZON_SP_P_irr_fn)
            inputZON_P_irr = open(self.inputZON_P_irr_fn, 'w')
            inputZON_P_irr.write('#\n')

            self.inputZON_SP_Pe_irr_fn = "inputZON_Pe_irr_stp.txt"
            self.inputZON_Pe_irr_fn = os.path.join(self.MM_ws, self.inputZON_SP_Pe_irr_fn)
            inputZON_Pe_irr = open(self.inputZON_Pe_irr_fn, 'w')
            inputZON_Pe_irr.write('#\n')

            self.inputZON_SP_PT_irr_fn = "inputZON_PT_irr_stp.txt"
            self.inputZON_PT_irr_fn = os.path.join(self.MM_ws, self.inputZON_SP_PT_irr_fn)
            inputZON_PT_irr = open(self.inputZON_PT_irr_fn, 'w')
            inputZON_PT_irr.write('#\n')

            self.input_SP_crop_irr_fn = "inputZON_crop_irr_stp.txt"
            self.inputcrop_irr_fn = os.path.join(self.MM_ws, self.input_SP_crop_irr_fn)
            inputcrop_irr = open(self.inputcrop_irr_fn, 'w')
            inputcrop_irr.write('#\n')

        try:
            for n in range(NMETEO):
                ExportResults1(P_veg_stp[n], inputZON_P_veg)
                ExportResults1(Eo_stp[n], inputZON_Eo)
                #VEG
                if NVEG>0:
                    for v in range(NVEG):
                        ExportResults1(PT_veg_stp[n][v], inputZON_PT)
                        ExportResults1(Pe_veg_stp[n][v], inputZON_Pe_veg)
                # SOIL
                if NSOIL>0:
                    for s in range(NSOIL):
                        ExportResults1(PE_stp[n][s], inputZON_PE)
                # CROP
                if NFIELD is not None:
                    for f in range(NFIELD):
                        ExportResults1(P_irr_stp[n][f], inputZON_P_irr)
                        ExportResults1(Pe_irr_stp[n][f], inputZON_Pe_irr)
                        ExportResults1(PT_irr_stp[n][f], inputZON_PT_irr)
            if NVEG>0:
                for v in range(NVEG):
                    ExportResults1(LAI_veg_stp[0][v], inputLAI_veg)
            if NFIELD is not None:
                for f in range(NFIELD):
                    ExportResults1(crop_irr_stp[0][f], inputcrop_irr)
        except Exception:
            self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nError in exporting output files in MF time processing.", stdout = stdout, report = report)

        print('Found %d days converted into %d stress periods.' % (sum(self.perlen), self.nper))

        inputZON_P_veg.close()
        inputZON_Pe_veg.close()
        inputLAI_veg.close()
        inputZON_PT.close()
        inputZON_PE.close()
        inputZON_Eo.close()
        if NFIELD is not None:
            inputZON_P_irr.close()
            inputZON_Pe_irr.close()
            inputZON_PT_irr.close()
            inputcrop_irr.close()

#####################################

    # NOTE Phase 1 (2026): runMF() (MODFLOW-NWT/mf2005 construction and run
    # via flopy.modflow) was removed with the Picard loop. Phase 3 adds
    # clsMF6 (flopy.mf6 + MODFLOW 6 API). See MARMITES_code_review.md.

##################
if __name__ == "__main__":
    print('\nWARNING!\nStart MARMITES-MODFLOW models using the script startMARMITES_v3.py\n')

# EOF            