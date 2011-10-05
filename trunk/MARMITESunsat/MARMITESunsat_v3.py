# -*- coding: utf-8 -*-
"""
MARMITES is a distributed depth-wise lumped-parameter model for
spatio-temporal assessment of water fluxes in the unsaturated zone.
MARMITES is a french word to
design a big cooking pot used by sorcerers for all kinds of experiments!
The main objective of MARMITES development is to partition the rainfall
in the several fluxes from the unsaturated and saturated zone (source
reservoir). It applies the concepts enunciated by Lubcynski (2010):                        ## TO BE UPDATED## 20100823
ET = Es + I + ETu + ETg
ETu = Eu + Tu
ETg = Eg + Tg
MARMITES is a modular algorithm that computes on a daily temporal
basis: interception, surface storage and runoff, evaporation from the
unsaturated and saturated zone (from MODFLOW), soil moisture storage,
aquifer recharge.
Input driving forces are rainfall and potential evapotranspiration daily
time series.
Calibration is made against soil moisture (DPLWATFLUX) and hydraulic
heads (MODFLOW).
It is link to MODFLOW2000 (mf2k) and the calibration is done
through PEST (Doherty 2010).
MARMITES imports and exports mf2k packages for spatial and temporal
discretization, intial heads, recharge, etc.
Thanks to Toews 2007 and flopy.                                         ## TO BE UPDATED## 20100823
For more documentation, see Frances et al 2010.                         ## TO BE UPDATED## 20100823
References:
Lub 2010
Frances et al 2010                                                      ## TO BE UPDATED## 20101015
"""

__author__ = "Alain Francés <frances08512@itc.nl>"
__version__ = "1.0"
__date__ = "November 2010"

import os
import numpy as np

class UNSAT:

    """
    unsat class: compute the water balance in soil (depth-wise, each reservoir
    correspond to a soil horizon (A, B, C and rock)
    ## TO BE UPDATED""
    _______________________________________________________________________________

    INPUTS
            PARAMETERS
                MAXIL       Interception looses
                Sm          Max. Soil moiture storage
                Sfc         Soil moiture storage at field capacity
                Sr          Residual Soil moiture storage
                Su_ini          Inicial soil moisture
                Ks          Saturated hydraulic condutivity
                Ss_max       Max. ponding capacity
            STATE VARIABLES
                RF           Daily rainfall
                PET         Daily evapotranspiration
    OUTPUTS
            RFe              Daily Excess rainfall
            ETu             Daily evapotranspiration
            Su               Daily soil moisture
            Rp              Daily percolation
            Ss            Daily ponding
            Ro              Daily runoff

    Provide daily values for the LINRES module
    ______________________________________________________________________________

    """

    def __init__(self, hnoflo):
        self.hnoflo = hnoflo

##        @ARTICLE{Shah2007,
##          author = {Shah, Nirjhar and Nachabe, Mahmood and Ross, Mark},
##          title = {Extinction Depth and Evapotranspiration from Ground Water under Selected
##        	Land Covers},
##          journal = {Ground Water},
##          year = {2007},
##          volume = {45},
##          pages = {329-338},
##          number = {3},
##        }
        # table 5, values for bare soil
        # dll [cm], y0 [], b [cm^-1], ext_d [cm]
        self.paramEg = {'sand'           : {'dll':16.0,'y0':0.000,'b':0.171, 'ext_d':50.0},
                       'loamy sand'      : {'dll':21.0,'y0':0.002,'b':0.130, 'ext_d':70.0},
                       'sandy loam'      : {'dll':30.0,'y0':0.004,'b':0.065, 'ext_d':130.0},
                       'sandy clay loam' : {'dll':30.0,'y0':0.006,'b':0.046, 'ext_d':200.0},
                       'sandy clay'      : {'dll':20.0,'y0':0.005,'b':0.042, 'ext_d':210.0},
                       'loam'            : {'dll':33.0,'y0':0.004,'b':0.028, 'ext_d':260.0},
                       'silty clay'      : {'dll':37.0,'y0':0.007,'b':0.046, 'ext_d':330.0},
                       'clay loam'       : {'dll':33.0,'y0':0.008,'b':0.027, 'ext_d':400.0},
                       'silt loam'       : {'dll':38.0,'y0':0.006,'b':0.019, 'ext_d':420.0},
                       'silt'            : {'dll':31.0,'y0':0.007,'b':0.021, 'ext_d':430.0},
                       'silty clay loam' : {'dll':40.0,'y0':0.007,'b':0.021, 'ext_d':450.0},
                       'clay'            : {'dll':45.0,'y0':0.006,'b':0.019, 'ext_d':620.0}
}

#####################

    def flux(self, RFe, PET, PE, E0, Zr_elev, VEGarea, HEADS, TopSoilLay, BotSoilLay, Tl, nsl, Sm, Sfc, Sr, Ks, Ss_max, Ss_ratio, Su_ini, Rp_ini, Ss_ini, EXF, dtwt, st, i, j, n, kTu_min, kTu_n, dt, dti):

        def surfwater(s_tmp,Sm,Ss_max, E0, i, j, n, dt):
            '''
            Ponding and surface runoff function
            '''
            if s_tmp > Sm:
                Ss_tmp = s_tmp-Sm
                if Ss_tmp > Ss_max or np.abs(Ss_tmp - Ss_max) < 1.0E-7:
                    Ro_tmp = (Ss_tmp-Ss_max)/dt
                    Ss_tmp = Ss_max
                else:
                    Ro_tmp = 0.0
                if Ss_tmp > E0 or np.abs(Ss_tmp - E0) < 1.0E-7:
                    Es_tmp = E0
                    Ss_tmp = Ss_tmp - E0*dt
                else:
                    Es_tmp = Ss_tmp/dt
                    Ss_tmp = 0.0
            else:
                Ss_tmp = 0.0
                Ro_tmp = 0.0
                Es_tmp = 0.0
            return Ss_tmp, Ro_tmp, Es_tmp

        ##################

        def perc(s_tmp, Sm, Sfc, Ks, s_lp1, Sm_sp1, i, j, n, dt):
            '''
            Percolation function
            '''
            if Sfc > s_tmp or np.abs(s_tmp-Sfc) < 1.0E-7:
                rp_tmp = 0.0
            elif s_tmp > Sfc and ( Sm > Sfc or np.abs(s_tmp-Sm) < 1.0E-7):
                Sg = (s_tmp-Sfc)/(Sm-Sfc) # gravitational water
                if (Ks*Sg*dt) > (s_tmp-Sfc):
                    rp_tmp = (s_tmp-Sfc)/dt
                else:
                    rp_tmp = Ks*Sg
            elif s_tmp > Sm:
                if (Ks*dt) > (Sm-Sfc):
                    rp_tmp = (Sm-Sfc)/dt
                else:
                    rp_tmp = Ks
            # verify the vol. available in deeper soil layer
#            if (rp_tmp*dt) > (Sm_sp1 - s_lp1):
#                rp_tmp = (Sm_sp1 - s_lp1)/dt
            return rp_tmp

        ##################

        def evp(s_tmp,Sfc,Sr, pet, i, j, n, dt):
            '''
            Actual evapotranspiration function
            '''
            if (s_tmp - Sr) <= 1E-9:
                evp_tmp = 0.0
            elif ((s_tmp - Sr) > 1E-9) and ((s_tmp < Sfc) or np.abs(s_tmp -Sfc) < 1.0E-7):
                Se = (s_tmp - Sr)/(Sfc - Sr)
                if (pet*Se*dt - (s_tmp - Sr)) > 1.0E-7:
                    evp_tmp = (s_tmp - Sr)/dt
                else:
                    evp_tmp = pet*Se
            elif (s_tmp - Sfc) > 1.0E-7:
                if (pet*dt-(Sfc - Sr)) > 1.0E-7:
                    evp_tmp = (Sfc - Sr)/dt
                else:
                    evp_tmp = pet
            return evp_tmp

        ##################
        if EXF < 0.0:
            print 'WARNING!\nEXF < 0.0, value %.6f corrected to 0.0.' % EXF
            EXF = 0.0
        SAT = np.zeros([nsl], dtype = bool)

         # Su_tmp from previous time step
        # if first time step, use Su_ini
        Su_tmp = np.zeros([nsl])
        Rexf_tmp = np.zeros([nsl])
        for l in range(nsl):
            Su_tmp[l] = Su_ini[l]*1.0

        # INFILTRATION
        # first soil layer
        Su_tmp[0] += RFe + Ss_ini
        # other soil layer, percolation from previous time step
        for l in range(1,nsl):
            Su_tmp[l] += Rp_ini[l-1]*dti

        # GW EXF
        Su_tmp[nsl-1] += EXF*dt

        # SOIL EXF
        llst = range(nsl)
        llst.reverse()
        for l in llst[:-1]:
            if Su_tmp[l] > Sm[l]*Tl[l]:
                Rexf_tmp[l] += Su_tmp[l] - Sm[l]*Tl[l]
                Su_tmp[l-1] += Rexf_tmp[l]
                Su_tmp[l] = Sm[l]*Tl[l]

        # SAT = True indicates saturation overland flow
        if EXF > 0.0:
            for l in llst:
                if Su_tmp[l] >= (Sm[l]*Tl[l]):
                    HEADS += Tl[l]
                    dtwt -= Tl[l]
                    SAT[l] = True
                else:
                    break

        # Ss and Ro
        surfwater_tmp = surfwater(Su_tmp[0],Sm[0]*Tl[0],Ss_max, E0*Ss_ratio, i, j, n, dt)
        Ss_tmp = surfwater_tmp[0]
        Ro_tmp = surfwater_tmp[1]
        Es_tmp = surfwater_tmp[2]
        Su_tmp[0] -= (Ss_tmp + dt*(Ro_tmp + Es_tmp))

        Rexf_tmp /= dt

        Rp_tmp = np.zeros([nsl])
        Tu_tmpZr = np.zeros([nsl,len(Zr_elev)])
        Tu_tmp = np.zeros([nsl])
        Eu_tmp = np.zeros([nsl])
        Su_pc_tmp = np.zeros([nsl])
        # Groundwater transpiration
        Tg_tmp_Zr = np.zeros([nsl, len(Zr_elev)])
        Tg_tmp = 0.0
        # soil layers
        for l in range(nsl):
            # Rp
            if l < (nsl-1):
                if SAT[l+1] == False:
                    Rp_tmp[l] = perc(Su_tmp[l],Sm[l]*Tl[l],Sfc[l]*Tl[l],Ks[l], Su_tmp[l+1], Sm[l+1]*Tl[l+1], i, j, n, dt)
            elif EXF == 0.0:
                Rp_tmp[l] = perc(Su_tmp[l],Sm[l]*Tl[l],Sfc[l]*Tl[l], Ks[l],0.0,1.0E6, i, j, n, dt)
            Su_tmp[l] -= Rp_tmp[l]*dt
            if SAT[l] == False:
                # Eu
                if Ss_tmp == 0.0:
                    if PE > 0.0:
                        Eu_tmp[l] = evp(Su_tmp[l],Sfc[l]*Tl[l], Sr[l]*Tl[l], PE, i, j, n, dt)
                        Su_tmp[l] -= Eu_tmp[l]*dt
                        PE -= Eu_tmp[l]
                # Tu
                for z in range(len(Zr_elev)):
                    if PET[z] > 0.0 :
                        if BotSoilLay[l] > Zr_elev[z]:
                            Tu_tmpZr[l,z] = evp(Su_tmp[l],Sfc[l]*Tl[l], Sr[l]*Tl[l], PET[z], i, j, n, dt)
                        elif TopSoilLay[l] > Zr_elev[z] :
                            PETc = PET[z]*(TopSoilLay[l]-Zr_elev[z])/Tl[l]
                            Tu_tmpZr[l,z] = evp(Su_tmp[l],Sfc[l]*Tl[l], Sr[l]*Tl[l], PETc, i, j, n, dt)
                        Tu_tmp[l] += Tu_tmpZr[l,z]*VEGarea[z]/100
                        PET[z] -= Tu_tmpZr[l,z]
                Su_tmp[l] -= Tu_tmp[l]*dt
            elif SAT[l] == True:
                # Tg
                for z in range(len(Zr_elev)):
                    if PET[z] > 0.0 :
                        if HEADS > Zr_elev[z]:
                            Tg_tmp_Zr[l,z] = evp(Sm[l]*Tl[l],Sfc[l]*Tl[l], Sr[l]*Tl[l], PET[z], i, j, n, dt)
                            Tg_tmp += Tg_tmp_Zr[l,z]*VEGarea[z]/100
                            PET[z] -= Tg_tmp_Zr[l,z]

        for l in range(nsl):
            Su_pc_tmp[l] = Su_tmp[l]/Tl[l]

        # GW evaporation, equation 17 of Shah et al 2007, see ref in the __init__
        if Ss_tmp == 0.0:
            if PE > 0.0:
                dtwt /= 10
                y0    = self.paramEg[st]['y0']
                b     = self.paramEg[st]['b']
                dll   = self.paramEg[st]['dll']
                ext_d = self.paramEg[st]['ext_d']
                if dtwt <= dll:
                    Eg_tmp = PE
                elif dtwt < ext_d:
                    Eg_tmp = PE*(y0 + np.exp(-b*(dtwt-dll)))
                else:
                    Eg_tmp = 0.0
                dtwt *= 10
            else:
                Eg_tmp = 0.0
        else:
            Eg_tmp = 0.0

        # Groundwater transpiration
        kTu_max = 1.0
        for z in range(len(Zr_elev)):
            if HEADS > Zr_elev[z]:
                for l in range(nsl):
                    if Su_pc_tmp[l] <> self.hnoflo:
                        if Su_pc_tmp[l] >= Sr[l]:
                            if Su_pc_tmp[l] <= Sm[l]:
                                kTu = kTu_min[z]+(kTu_max-kTu_min[z])*np.power(1-np.power(np.abs(Su_pc_tmp[l]-Sm[l])/(Sm[l]-Sr[l]),kTu_n[z]),1/kTu_n[z])
                            else:
                                kTu = kTu_max
                        else:
                            kTu = kTu_min[z]
                        Tg_tmp_Zr[l,z] = Tu_tmpZr[l,z]*(1/kTu-1)
                        if Tg_tmp_Zr[l,z] > PET[z]:
                            Tg_tmp_Zr[l,z] = PET[z]
                        PET[z] -= Tg_tmp_Zr[l,z]
                        Tg_tmp += (Tg_tmp_Zr[l,z]*VEGarea[z]/100)

        return Es_tmp, Ss_tmp, Ro_tmp, Rp_tmp, Eu_tmp, Tu_tmp, Su_tmp, Su_pc_tmp, Eg_tmp, Tg_tmp, HEADS, dtwt, SAT, Rexf_tmp
        del Es_tmp, Ss_tmp, Ro_tmp, Rp_tmp, Eu_tmp, Tu_tmp, Su_tmp, Su_pc_tmp, Su_ini, Eg_tmp, Tg_tmp, dtwt

#####################

    def run(self, i, j, n,
                  nsl, st, Sm, Sfc, Sr, Su_ini, Ss_ini, Rp_ini, botm_l0, TopSoilLay, BotSoilLay, Tl, Ks, Ss_max, Ss_ratio, HEADS, EXF,
                  RF, E0, PETveg, RFeveg, PEsoil, VEGarea, Zr,
                  nstp, perlen, dti, hdry,
                  kTu_min, kTu_n):

        # Output initialisation
        Ttotal = len(RF)
        # PET for the vegetation patchwork
        PET_tot = np.zeros([len(PETveg[0])], dtype = float)
        # PE for the remaining bare soil
        PE_tot = np.zeros([Ttotal], dtype = float)
        # RFe for the vegetation patchwork
        RFe_tot = np.zeros([len(RFeveg[0])], dtype = float)
        # vegetation interception (weigthed average patchwork)
        INTER_tot = np.zeros([Ttotal], dtype = np.float)
        # Ponding storage
        Ss = np.zeros([Ttotal], dtype = np.float)
        # Ponding storage change
        dSs = np.zeros([Ttotal], dtype = np.float)
        # Surface runoff (hortonian and saturated overland flow)
        Ro = np.zeros([Ttotal], dtype = np.float)
        # Percolation
        Rp = np.zeros([Ttotal,nsl], dtype = np.float)
        # Exfiltration
        Rexf = np.zeros([Ttotal,nsl], dtype = np.float)
        #Soil moisture storage
        Su = np.zeros([Ttotal,nsl], dtype = np.float)
        #Soil moisture storage in volumetric percent of saturation
        Su_pc = np.zeros([Ttotal,nsl], dtype = np.float)
        #Soil moisture storage changes
        dSu = np.zeros([Ttotal,nsl], dtype = np.float)
        #Total soil moisture storage changes
        dSu_tot = np.zeros([Ttotal], dtype = np.float)
        #Actual evaporation from unsaturated zone
        Eu = np.zeros([Ttotal,nsl], dtype = np.float)
        #Actual transpiration from unsaturated zone
        Tu = np.zeros([Ttotal,nsl], dtype = np.float)
        #Total actual evapotranspiration from unsaturated zone
        ETu_tot = np.zeros([Ttotal], dtype = np.float)
        #Saturation overland flow
        SAT = np.zeros([Ttotal, nsl], dtype = np.int)
        #Actual evapotranspiration from surface
        Es = np.zeros([Ttotal], dtype = np.float)
        # HEADS
        HEADS_corr = np.zeros([Ttotal], dtype = np.float)
        # dtwt
        dtwt = np.zeros([Ttotal], dtype = np.float)
        # uzthick
        uzthick = np.zeros([Ttotal], dtype = np.float)
        # MASS BALANCE
        MB = np.zeros([Ttotal], dtype = np.float)
        # MASS BALANCE each soil layer
        MB_l = np.zeros([Ttotal, nsl], dtype = np.float)
        # bare soil GW evaporation
        Eg = np.zeros([Ttotal], dtype = np.float)
        # GW transpiration
        Tg = np.zeros([Ttotal], dtype = np.float)
        # GW evapotranspiration
        ETg = np.zeros([Ttotal], dtype = np.float)
        # output arrays
        nflux1 = 20
        nflux2 = 9
        results1 = np.zeros([Ttotal,nflux1])
        results2 = np.zeros([Ttotal,nsl,nflux2])

        Zr_elev = []
        for z in range(len(Zr)):
            Zr_elev.append(TopSoilLay[0] - Zr[z]*1000.0)

        dt = float(perlen)/float(nstp)

        # PROCESSING THE WHOLE DATA SET
        for t in range(int(nstp)):    # t: current time step

            # Preprocessing of PET/PE/INTER and RFe
            # for different vegetation
            SOILarea = 100
            for v in range(len(PETveg)):
                if VEGarea[v] != self.hnoflo:
                    RFe_tot[t] += RFeveg[v,t]*VEGarea[v]/100
                    PET_tot[t] += PETveg[v,t]*VEGarea[v]/100
                    SOILarea   -= VEGarea[v]
            RFe_tot[t]   += RF[t]*SOILarea/100
            INTER_tot[t]  = RF[t] - RFe_tot[t]
            PE_tot[t]     = PEsoil[t]*SOILarea/100

            # handle drycell
            if HEADS[t] > (hdry-1E3):
                HEADS_tmp = botm_l0*1000
            else:
                HEADS_tmp = HEADS[t]*1000.0


            # dtwt and uzthick
            uzthick[t] = BotSoilLay[nsl-1] - HEADS_tmp
            if EXF[t] <= 0.0:  #BotSoilLay[nsl-1] > HEADS_tmp:
                dtwt[t] = TopSoilLay[0] - HEADS_tmp
            elif EXF[t] > 0.0:
                dtwt[t] = sum(Tl)
                HEADS_tmp = BotSoilLay[nsl-1]

            # for the first TS, S_ini is expressed in % and has to be converted in mm
            if n == 0 and t == 0:
                Su_ini = Su_ini * Tl

            # fluxes
            Es_tmp, Ss_tmp, Ro_tmp, Rp_tmp, Eu_tmp, Tu_tmp, Su_tmp, Su_pc_tmp, Eg_tmp, Tg_tmp, HEADS_tmp, dtwt_tmp, SAT_tmp, Rexf_tmp = self.flux(RFe_tot[t], PETveg[:,t], PE_tot[t], E0[t], Zr_elev, VEGarea, HEADS_tmp, TopSoilLay, BotSoilLay, Tl, nsl, Sm, Sfc, Sr, Ks, Ss_max, Ss_ratio, Su_ini, Rp_ini, Ss_ini, EXF[t], dtwt[t], st, i, j, n, kTu_min, kTu_n, dt, dti)

            # fill the output arrays
            Ss[t] = Ss_tmp
            Ro[t]   = Ro_tmp
            Es[t]   = Es_tmp
            Eg[t]   = Eg_tmp
            Tg[t]   = Tg_tmp
            HEADS_corr[t] = HEADS_tmp
            dtwt[t] = dtwt_tmp
            Su[t,:]    = Su_tmp[:]
            Su_pc[t,:] = Su_pc_tmp[:]
            Rp[t,:]    = Rp_tmp[:]
            Rexf[t,:]  = Rexf_tmp[:]
            Eu[t,:]    = Eu_tmp[:]
            Tu[t,:]    = Tu_tmp[:]
            SAT[t,:]   = SAT_tmp[:]
            ETg[t] = Eg[t] + Tg[t]
            dSs_MB = (Ss_ini - Ss[t])/dt
            if (Ss[t] > Ss_ini):
                dSs[t] = (Ss[t] - Ss_ini)/dt
            # compute the water mass balance (MB) in the unsaturated zone
            Rp_in_MB  = 0.0
            Rp_out_MB = 0.0
            Eu_MB    = 0.0
            Tu_MB    = 0.0
            for l in range(nsl):
                Eu_MB     += Eu[t,l]
                Tu_MB     += Tu[t,l]
                if l < (nsl-1):
                    Rp_in_MB  += Rp_ini[l]*dti
                Rp_out_MB += Rp[t,l]*dt
                dSu[t,l] = (Su_ini[l] - Su[t,l])/dt
                dSu_tot[t] += dSu[t,l]
            dRp_tot = (Rp_in_MB - Rp_out_MB)/dt
            ETu_tot[t] = Eu_MB + Tu_MB

            # MASS BALANCE COMPUTING
            if nsl > 1:
                # surficial soil layer
                l = 0
                MB_l[t,l] = (RFe_tot[t] + dSu[t,l] + dSs_MB + Rexf[t,l+1]) - (Rp[t,l] + Eu[t,l] + Tu[t,l] + Ro_tmp + Es_tmp)
                # intermediate soil layers
                llst = range(1,nsl-1)
                for l in llst:
                    MB_l[t,l] = (Rp_ini[l-1]*dti - Rp[t,l]*dt)/dt + dSu[t,l] + Rexf[t,l+1] - (Eu[t,l] + Tu[t,l] + Rexf[t,l])
                # last soil layer
                l = nsl-1
                MB_l[t,l] = (Rp_ini[l-1]*dti - Rp[t,l]*dt)/dt + dSu[t,l] + EXF[t] - (Eu[t,l] + Tu[t,l] + Rexf[t,l])
            else:
                # only one soil layer
                l = 0
                MB_l[t,l] = (RFe_tot[t] + dSu[t,l] + dSs_MB + EXF[t]) - (Rp[t,l] + Eu[t,l] + Tu[t,l] + Ro_tmp + Es_tmp)
                # last soil layer
            # total mass balance for the unsaturated zone
            MB[t] = RFe_tot[t] + dSs_MB + dSu_tot[t] + dRp_tot + EXF[t] - (Ro[t] + Es[t] + Eu_MB + Tu_MB)

            # export list
            # indexes of the HDF5 output arrays
            # index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iEXF':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSs':13, 'iETg':14, 'iETu':15, 'idSu':16, 'iHEADScorr':17, 'idtwt':18, 'iuzthick':19}
            results1[t,:] = [RF[t], PET_tot[t], PE_tot[t], RFe_tot[t], Ss[t], Ro[t], EXF[t], Es[t], MB[t], INTER_tot[t], E0[t], Eg[t], Tg[t], dSs[t], ETg[t], ETu_tot[t], dSu_tot[t], HEADS_corr[t]/1000.0, -dtwt[t]/1000.0, uzthick[t]/1000.0]
            # index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'iRexf':4, 'idSu':5, 'iSu':6, 'iSAT':7, 'iMB_l':8}
            for l in range(nsl):
                results2[t,l,:] = [Eu[t,l], Tu[t,l], Su_pc[t,l], Rp[t,l], Rexf[t,l], dSu[t,l], Su[t,l], SAT[t,l], MB_l[t,l]]

            # initial values for next time step
            dti = dt*1.0
            Ss_ini = Ss[t]*1.0
            Su_ini = Su[t,:]*1.0
            Rp_ini = Rp[t,:]*1.0

        return results1, results2

        del RF, PET_tot, PE_tot, RFe_tot, Ss, Ro, EXF, Es, MB, INTER_tot, E0, Eg, Tg, dSs, ETg, ETu_tot, dSu_tot, HEADS, HEADS_corr, dtwt
        del Eu, Tu, Su_pc, Rp, dSu, Su
        del results1, results2

class SATFLOW:
    """
    SATFLOW: SATurated FLOW
    Calculate water level fluctuations in function of recharge
    _________________________
    INPUTS
            PARAMETERS
                hi
                h0          base level
                RC          Saturated recession constant
                STO         Aquifer Storage capacity
            STATE VARIABLES
                R           Daily recharge
    OUTPUTS
            h               Daily water level
     _________________________
    """

    def run(self, R, hi, h0, RC, STO):
        h1 = np.zeros([len(R)], dtype=np.float)
        h_tmp = (hi*1000.0 + R[0]/STO - hi*1000.0/RC)
        h1[0] = h_tmp
        for t in range(1,len(R)):
            h_tmp = h1[t-1] + R[t]/STO -h1[t-1]/RC
            h1[t] = h_tmp

        h = np.zeros([len(R)], dtype=np.float)
        for t in range(0,len(R)):
            h_tmp = (h1[t] + h0*1000.0)/1000.0
            h[t] = h_tmp

        return h
        del R, h1, h_tmp, h

#EOF#