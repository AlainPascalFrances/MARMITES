# -*- coding: utf-8 -*-
"""
DLPWATFLUX stands for Distributed depth-wise Lumped-Parameter model for
spatio-temporal assessment of WATer FLUXes in the unsaturated zone.
The familiar and friendly name of DLPWATFLUX is MARMITES, a french word to
design a big cooking pot used by sorcerers for all kinds of experiments!
The main objective of DPLWATFLUX development is to partition the rainfall
in the several fluxes from the unsaturated and saturated zone (source
reservoir). It applies the concepts enunciated by Lubcynski (2010):                        ## TO BE UPDATED## 20100823
ET = Es + I + ETu + ETg
ETu = Eu + Tu
ETg = Eg + Tg
DLPWATFLUX is a modular algorithm that computes on a daily temporal
basis: interception, surface storage and runoff, evaporation from the
unsaturated and saturated zone (from MODFLOW), soil moisture storage,
aquifer recharge.
Input driving forces are rainfall and potential evapotranspiration daily
time series.
Calibration is made against soil moisture (DPLWATFLUX) and hydraulic
heads (MODFLOW).
It is link to MODFLOW2000 (mf2k) and the calibration is done
through PEST (Doherty 2010).
DLPWATFLUX imports and exports mf2k packages for spatial and temporal
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
        self.paramEg = {'sand'           : {'dll':160.0,'y0':0.000,'b':0.0171, 'ext_d':500.0},
                       'loamy sand'      : {'dll':210.0,'y0':0.002,'b':0.0130, 'ext_d':700.0},
                       'sandy loam'      : {'dll':300.0,'y0':0.004,'b':0.0065, 'ext_d':1300.0},
                       'sandy clay loam' : {'dll':300.0,'y0':0.006,'b':0.0046, 'ext_d':2000.0},
                       'sandy clay'      : {'dll':200.0,'y0':0.005,'b':0.0042, 'ext_d':2100.0},
                       'loam'            : {'dll':330.0,'y0':0.004,'b':0.0028, 'ext_d':2600.0},
                       'silty clay'      : {'dll':370.0,'y0':0.007,'b':0.0046, 'ext_d':3300.0},
                       'clay loam'       : {'dll':330.0,'y0':0.008,'b':0.0027, 'ext_d':4000.0},
                       'silt loam'       : {'dll':380.0,'y0':0.006,'b':0.0019, 'ext_d':4200.0},
                       'silt'            : {'dll':310.0,'y0':0.007,'b':0.0021, 'ext_d':4300.0},
                       'silty clay loam' : {'dll':400.0,'y0':0.007,'b':0.0021, 'ext_d':4500.0},
                       'clay'            : {'dll':450.0,'y0':0.006,'b':0.0019, 'ext_d':6200.0}
}

#####################

    def flux(self, RFe, PET, PE, Su, Rp, Zr_elev, VEGarea, HEADS, HEADS_ini, Dltop, Dlbot, Dl, Dl_ini, nsl, Sm, Sfc, Sr, Ks, Ss_max, Ss_ratio, t, Su_ini, Rp_ini, DRN, DRN_ini, Ss_ini, E0, dtwt, st, i, j, n, kTu_min, kTu_n, dt, dti):

        def surfwater(s_tmp,Sm,Ss_max, E0, i, j, n, dt):
            '''
            Ponding and surface runoff function
            '''
            if (s_tmp-Sm) > 1E-9:
                Ss_tmp = s_tmp-Sm
                if (Ss_tmp-Ss_max) >= 1E-9:
                    Ro_tmp = (Ss_tmp-Ss_max)/dt
                    Ss_tmp = Ss_max
                else:
                    Ro_tmp = 0.0
                if (Ss_tmp - E0) >= 1E-9:
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

        def perc(s_tmp,Sm,Sfc,Ks, s_lp1, Sm_sp1, i, j, n, dt):
            '''
            Percolation function
            '''
            global rp_tmp
            if (s_tmp-Sfc) <= 1E-9:
                rp_tmp = 0.0
            elif (s_tmp-Sfc) > 1E-9 and (s_tmp-Sm) <= 1E-9:
                Sg = (s_tmp-Sfc)/(Sm-Sfc) # gravitational water
                if (Ks*Sg*dt-(s_tmp-Sfc)) > 1E-9:
                    rp_tmp = (s_tmp-Sfc)/dt
                else:
                    rp_tmp = Ks*Sg
            elif (s_tmp-Sm) > 1E-9:
                if (Ks*dt-(Sm-Sfc)) > 1E-9:
                    rp_tmp = (Sm-Sfc)/dt
                else:
                    rp_tmp = Ks
            # verify the vol. available in deeper soil layer
            if (rp_tmp*dt - (Sm_sp1 - s_lp1)) > 1E-9:
                rp_tmp = (Sm_sp1 - s_lp1)/dt
            return rp_tmp

        ##################

        def evp(s_tmp,Sfc,Sr, pet, i, j, n, dt):
            '''
            Actual evapotranspiration function
            '''
            # Percent. of soil saturation
            global evp_tmp
            Se = (s_tmp - Sr)/(Sfc - Sr)
            if (s_tmp - Sr) <= 1E-9:
                evp_tmp = 0.0
            elif ((s_tmp - Sr) > 1E-9) and ((s_tmp - Sfc) <= 1E-9):
                if (pet*Se*dt - (s_tmp - Sr)) > 1E-9:
                    evp_tmp = (s_tmp - Sr)/dt
                else:
                    evp_tmp = pet*Se
            elif (s_tmp - Sfc) > 1E-9:
                if (pet*dt-(Sfc - Sr)) > 1E-9:
                    evp_tmp = (Sfc - Sr)/dt
                else:
                    evp_tmp = pet
            return evp_tmp

        ##################

##        def corr(flx, S_i, dt, Sr, Dl):
##            flxdt = flx*dt
##            S_f = S_i - flxdt
##            if S_f <= 0:
##                S_f = Sr*Dl
##                flx = (S_i-S_f)/dt
##            return flx, S_f
##
##        ##################

        # Su_tmp from previous time step
        # if first time step, use Su_ini
        Su_tmp = np.zeros([nsl])
        for l in range(nsl):
            Su_tmp[l] = Su_ini[l]

        llst = range(nsl-1)
        llst.reverse()

        llst_full = range(nsl)
        llst_full.reverse()

        # percolation from previous time step
        # if first time step, use Rp_ini
        for l in llst_full[:-1]:
            if t == 0:
                Rpprev = Rp_ini[l-1]*dti
            else:
                Rpprev = Rp[t-1,l-1]*dti
            if Dl[l] > 0.0:
                if Rpprev - (Sm[l]*Dl[l]-Su_tmp[l]) > 1E-9:
                    Rpprev -= (Sm[l]*Dl[l]-Su_tmp[l])
                    Su_tmp[l] = Sm[l]*Dl[l]
                    Su_tmp[l-1] += Rpprev
                else:
                    Su_tmp[l] += Rpprev
            else:
                Su_tmp[l-1] += Rpprev

        DRN_corr = 0.0
        DRN_ini_TS = DRN
        if DRN == 0.0:
##            if DRN_ini > 0.0:
##                DRN_corr = (Sfc[nsl-1]-Sr[nsl-1])*Dl[nsl-1]/dt
##                Su_tmp[nsl-1] = Sfc[nsl-1]*Dl[nsl-1]
##                Su_ini[nsl-1] = Su_tmp[nsl-1]
##            else:
                if (HEADS_ini - HEADS) > 1E-9:
                    DRN_corr = Sfc[nsl-1]*(HEADS_ini - HEADS)/dt
                    Su_tmp[nsl-1] += DRN_corr*dt
                    Su_ini[nsl-1] = Su_tmp[nsl-1]
        else:
        #if DRN > 0.0:
            # GW seepage
            DRN += Su_tmp[nsl-1]-Sr[nsl-1]*Dl_ini[nsl-1]
            for l in llst:
                if (DRN - (Sm[l]*Dl[l]-Su_tmp[l])) > 1E-9:
                    DRN -= (Sm[l]*Dl[l]-Su_tmp[l])
                    Su_tmp[l] = Sm[l]*Dl[l]
                else:
                    Su_tmp[l] += DRN
                    DRN = 0.0
            Su_tmp[nsl-1] = None

        # initialization first soil layer
        Su_tmp[0] += RFe*dt + DRN + Ss_ini

        # Ss and Ro
        surfwater_tmp = surfwater(Su_tmp[0],Sm[0]*Dl[0],Ss_max, E0*Ss_ratio, i, j, n, dt)
        Ss_tmp = surfwater_tmp[0]
        Ro_tmp = surfwater_tmp[1]
        Es_tmp = surfwater_tmp[2]
        Su_tmp[0] -= (Ss_tmp + dt*(Ro_tmp + Es_tmp))

        SAT = np.zeros([nsl], dtype = bool)
        for l in llst_full:
            if DRN_ini_TS > 0.0  and Su_tmp[l] == (Sm[l]*Dl[l]): # > 1E-9:
                HEADS += Dl[l]
                dtwt -= Dl[l]
                SAT[l] = True
            else:
                break

        Rp_tmp = np.zeros([nsl])
        R_tmp = 0.0
        Tu_tmpZr = np.zeros([nsl,len(Zr_elev)])
        Tu_tmp = np.zeros([nsl])
        Eu_tmp = np.zeros([nsl])
        Su_pc_tmp = np.zeros([nsl])
        # Groundwater transpiration
        Tg_tmpZr = np.zeros([nsl, len(Zr_elev)])
        Tg_tmp = 0.0
        kTu_max = 1.0
        # soil layers and vadose zone layer
        for l in range(nsl):
            # Rp
            if l < (nsl-1):
                if Dl[l+1] > 0.0:
                    Rp_tmp[l] = perc(Su_tmp[l],Sm[l]*Dl[l],Sfc[l]*Dl[l],Ks[l], Su_tmp[l+1], Sm[l+1]*Dl[l+1], i, j, n, dt)
                    Su_tmp[l] -= Rp_tmp[l]*dt
                    #Rp_tmp[l], Su_tmp[l] = corr(Rp_tmp[l], Su_tmp[l], dt, Sr[l], Dl[l])
            else:
                if Dl[l] > 0.0 and DRN_corr == 0.0:
                    Rp_tmp[l] = perc(Su_tmp[l],Sm[l]*Dl[l],Sfc[l]*Dl[l],Ks[l],0.0,1.0E6, i, j, n, dt)
                    Su_tmp[l] -= Rp_tmp[l]*dt
                    #Rp_tmp[l], Su_tmp[l] = corr(Rp_tmp[l], Su_tmp[l], dt, Sr[l], Dl[l])
                    R_tmp = Rp_tmp[l]
            if SAT[l] == False:
                if (Dlbot[l] - HEADS) > 1E-9:
                    Dl_tmp = Dl[l]
                elif (Dltop[l] - HEADS) > 1E-9:
                    Dl_tmp = Dltop[l] - HEADS
                if Dl[l] > 0.0:
                    # Eu
                    if Ss_tmp == 0.0:
                        if PE > 0.0:
                            Eu_tmp[l] = evp(Su_tmp[l],Sfc[l]*Dl_tmp,Sr[l]*Dl_tmp, PE, i, j, n, dt)
                            #Eu_tmp[l], Su_tmp[l] = corr(Eu_tmp[l], Su_tmp[l], dt, Sr[l], Dl_tmp)
                            Su_tmp[l] -= Eu_tmp[l]*dt
                            PE -= Eu_tmp[l]
                    # Tu
                    for z in range(len(Zr_elev)):
                        if PET[z] > 0.0:
                            if HEADS > Zr_elev[z]:
                                Tu_tmpZr[l,z] = evp(Su_tmp[l],Sfc[l]*Dl_tmp,Sr[l]*Dl_tmp, PET[z], i, j, n, dt)
                            elif Dltop[l] > Zr_elev[z]:
                                PETc = PET[z]*(Dltop[l]-Zr_elev[z])/Dl_tmp
                                Tu_tmpZr[l,z] = evp(Su_tmp[l],Sfc[l]*Dl_tmp,Sr[l]*Dl_tmp, PETc, i, j, n, dt)
                            Tu_tmp[l] += Tu_tmpZr[l,z]*VEGarea[z]/100.0
                            PET[z] -= Tu_tmpZr[l,z]
                            #Tu_tmp[l], Su_tmp[l] = corr(Tu_tmp[l], Su_tmp[l], dt, Sr[l], Dl_tmp)
                    Su_tmp[l] -= Tu_tmp[l]*dt
            elif SAT[l] == True:
               if Dl[l] > 0.0:
                for z in range(len(Zr_elev)):
                    Tg_tmpZr[l,z] = evp(Sm[l]*Dl[l],Sfc[l]*Dl[l],Sr[l]*Dl[l], PET[z], i, j, n, dt)
                    Tg_tmp += Tg_tmpZr[l,z]*VEGarea[z]/100.0
                    PET[z] -= Tg_tmpZr[l,z]

        for l in range(nsl):
            if Dl[l] > 0.0:
                Su_pc_tmp[l] = Su_tmp[l]/Dl[l]
            else:
                Su_pc_tmp[l] = Sm[l]
                Su_tmp[l] = Sm[l]*Dl[l]

        # GW evaporation, equation 17 of Shah et al 2007, see ref in the __init__
        if Ss_tmp == 0.0:
            if PE > 0.0:
                y0  = self.paramEg[st]['y0']
                b   = self.paramEg[st]['b']
                dll = self.paramEg[st]['dll']
                if dtwt<=dll:
                    Eg_tmp = PE
                elif dtwt > 0.0:
                    Eg_tmp = PE*(y0 + np.exp(-b*(dtwt-dll)))
                else:
                    Eg_tmp = 0.0
                if Eg_tmp/PE < 0.5/100.0:
                    Eg_tmp = 0.0
            else:
                Eg_tmp = 0.0
        else:
            Eg_tmp = 0.0

        for z in range(len(Zr_elev)):
            if ((HEADS - Dlbot[nsl-1]) > 1E-9) and ((HEADS - Zr_elev[z]) > 1e-9):
                for l in range(nsl):
                    if Su_pc_tmp[l] >= Sr[l]:
                        if Su_pc_tmp[l] <= Sm[l]:
                            kTu = kTu_min[z]+(kTu_max-kTu_min[z])*np.power(1-np.power(np.abs(Su_pc_tmp[l]-Sm[l])/(Sm[l]-Sr[l]),kTu_n[z]),1/kTu_n[z])
                        else:
                            kTu = kTu_max
                    else:
                        kTu = kTu_min[z]
                    Tg_tmpZr[l,z] = Tu_tmpZr[l,z]*(1/kTu-1)
                    if Tg_tmpZr[l,z] > PET[z]:
                        Tg_tmpZr[l,z] = PET[z]
                    PET[z] -= Tg_tmpZr[l,z]
                    Tg_tmp += (Tg_tmpZr[l,z]*VEGarea[z]/100)   # *(Dl[l]/sum(Dl))

        return Es_tmp, Ss_tmp, Ro_tmp, Rp_tmp, Eu_tmp, Tu_tmp, Su_tmp, Su_pc_tmp, Su_ini, R_tmp, Eg_tmp, Tg_tmp, HEADS, DRN_corr, dtwt
        del Es_tmp, Ss_tmp, Ro_tmp, Rp_tmp, Eu_tmp, Tu_tmp, Su_tmp, Su_pc_tmp, Su_ini, R_tmp, Eg_tmp, Tg_tmp, dtwt

#####################

    def run(self, i, j, n,
                  nsl, st, slprop, Sm, Sfc, Sr, Su_ini, Ss_ini, Rp_ini, D, Dl_ini, Ks, Ss_max, Ss_ratio,
                  TopAquif, HEADS, HEADS_ini, DRN, DRN_ini, cf,
                  RF, E0, PETveg, RFeveg, PEsoil, VEGarea, Zr,
                  nstp, perlen, dti, hdry,
                  kTu_min, kTu_n):

        # Output initialisation
        Ttotal=len(RF)
        # PET for the vegetation patchwork
        PET_tot=np.zeros([len(PETveg[0])], dtype=float)
        # PE for the remaining bare soil
        PE_tot=np.zeros([Ttotal], dtype=float)
        # RFe for the vegetation patchwork
        RFe_tot=np.zeros([len(RFeveg[0])], dtype=float)
        # vegetation interception (weigthed average patchwork)
        INTER_tot=np.zeros([Ttotal], dtype=np.float)
        # Ponding storage
        Ss=np.zeros([Ttotal], dtype=np.float)
        # Ponding storage change
        dSs=np.zeros([Ttotal], dtype=np.float)
        # Surface runoff (hortonian and saturated overland flow)
        Ro=np.zeros([Ttotal], dtype=np.float)
        #Deep percolation
        Rp=np.zeros([Ttotal,nsl], dtype=np.float)
        #Soil moisture storage
        Su=np.zeros([Ttotal,nsl], dtype=np.float)
        #Soil moisture storage in volumetric percent of saturation
        Su_pc=np.zeros([Ttotal,nsl], dtype=np.float)
        #Soil moisture storage changes
        dSu=np.zeros([Ttotal,nsl], dtype=np.float)
        #Total soil moisture storage changes
        dSu_tot=np.zeros([Ttotal], dtype=np.float)
        #Actual evaporation from unsaturated zone
        Eu=np.zeros([Ttotal,nsl], dtype=np.float)
        #Actual transpiration from unsaturated zone
        Tu=np.zeros([Ttotal,nsl], dtype=np.float)
        #Total actual evapotranspiration from unsaturated zone
        ETu_tot=np.zeros([Ttotal], dtype=np.float)
        #Actual evapotranspiration from surface
        Es=np.zeros([Ttotal], dtype=np.float)
        # GW seepage into soil zone
        SEEPAGE=np.zeros([Ttotal], dtype=int)
        # HEADS
        HEADScorr=np.zeros([Ttotal], dtype=int)
        #MASS BALANCE
        MB=np.zeros([Ttotal], dtype=float)
        # bare soil GW evaporation
        Eg=np.zeros([Ttotal], dtype=np.float)
        # GW transpiration
        Tg=np.zeros([Ttotal], dtype=np.float)
        # GW recharge
        R=np.zeros([Ttotal], dtype=np.float)
        # GW net recharge
        Rn=np.zeros([Ttotal], dtype=np.float)
        # GW evapotranspiration
        ETg=np.zeros([Ttotal], dtype=np.float)
        # output arrays
        nflux1 = 22
        nflux2 = 6
        results1 = np.zeros([Ttotal,nflux1])
        results2 = np.zeros([Ttotal,nsl,nflux2])

#        if D < 0.05: #correct soil thickness less than 5cm
#            D=0.05
        # conversion from m to mm
        D = D*1000.0
        TopAquif = TopAquif*1000.0
        # elevation of bottom of soil reservoir
        Dbot = TopAquif
        # topography elevation
        Dtop = TopAquif + D

        # thickness of soil layers
        Dl = np.zeros([nsl], dtype=float)
        for l in range(nsl-1):
            Dl[l] = D*slprop[l]
        # elevation of top and bottom of soil layers
        Dltop = np.zeros([nsl], dtype=float)
        Dlbot = np.zeros([nsl], dtype=float)
        for l in range(nsl-1):
            if l==0:
                Dltop[l] = Dtop
                Dlbot[l] = Dtop-Dl[l]
            else:
                Dltop[l] = Dlbot[l-1]
                Dlbot[l] = Dltop[l]-Dl[l]
        Dltop[nsl-1] = Dlbot[nsl-2]

        Zr_elev = []
        for z in range(len(Zr)):
            Zr_elev.append(Dtop - Zr[z]*1000.0)

        dt = float(perlen)/float(nstp)

        # PROCESSING THE WHOLE DATA SET
        for t in range(int(nstp)):    # t: current time step

            # Preprocessing of PET/PE/INTER and RFe
            # for different vegetation
            SOILarea = 100
            for v in range(len(PETveg)):
                if VEGarea[v] <> self.hnoflo:
                    RFe_tot[t] += RFeveg[v,t]*VEGarea[v]/100
                    PET_tot[t] += PETveg[v,t]*VEGarea[v]/100
                    SOILarea   -= VEGarea[v]
            RFe_tot[t]   += RF[t]*SOILarea/100
            INTER_tot[t]  = RF[t] - RFe_tot[t]
            PE_tot[t]     = PEsoil[t]*SOILarea/100

            # handle drycell
            if HEADS[t] > (hdry-1E3):
                HEADS_tmp = Dbot - 100000.0
            else:
                HEADS_tmp = HEADS[t]*1000.0
            if HEADS_ini > (hdry-1E3):
                HEADS_ini = Dbot - 100000.0
            else:
                HEADS_ini *= 1000.0
            dtwt = Dtop-HEADS_tmp

            DRN_ini = DRN_ini * dti

            # heads below soil bottom
            if DRN[t] == 0.0:
                DRN_tmp = 0
                # bottom boundary
                Dlbot[nsl-1] = HEADS_tmp
                Dl[nsl-1] = Dltop[nsl-1] - HEADS_tmp
                if n == 0 and t == 0:
                    Dl_ini = Dl
                    for l in range(nsl):
                        if Dl[l] > 0.0:
                            Su_ini[l] = Su_ini[l]*Dl[l]
                        else:
                            Su_ini[l] = Sm[l]*Dl[l]
                # fluxes
                Es_tmp, Ss_tmp, Ro_tmp, Rp_tmp, Eu_tmp, Tu_tmp, Su_tmp, Su_pc_tmp, Su_ini, R_tmp, Eg_tmp, Tg_tmp, HEADS_tmp, DRN_corr, dtwt = self.flux(RFe_tot[t], PETveg[:,t], PE_tot[t], Su, Rp, Zr_elev, VEGarea, HEADS_tmp, HEADS_ini, Dltop, Dlbot, Dl, Dl_ini, nsl, Sm, Sfc, Sr, Ks, Ss_max, Ss_ratio, t, Su_ini, Rp_ini, DRN_tmp, DRN_ini, Ss_ini, E0[t], dtwt, st, i, j, n, kTu_min, kTu_n, dt, dti)
            # heads above soil bottom
            else:
                DRN_tmp = DRN[t]*cf*dt  # cbc are expressed in volume/unit time
                # bottom boundary
                Dlbot[nsl-1] = TopAquif
                Dl[nsl-1] = 0.0
                if n == 0 and t == 0:
                    Dl_ini = Dl
                    for l in range(nsl):
                        if Dl[l] > 0.0:
                            Su_ini[l] = Su_ini[l]*Dl[l]
                        else:
                            Su_ini[l] = Sm[l]*Dl[l]
                # fluxes
                Es_tmp, Ss_tmp, Ro_tmp, Rp_tmp, Eu_tmp, Tu_tmp, Su_tmp, Su_pc_tmp, Su_ini, R_tmp, Eg_tmp, Tg_tmp, HEADS_tmp, DRN_corr, dtwt = self.flux(RFe_tot[t], PETveg[:,t], PE_tot[t], Su, Rp, Zr_elev, VEGarea, HEADS_tmp, HEADS_ini, Dltop, Dlbot, Dl, Dl_ini, nsl, Sm, Sfc, Sr, Ks, Ss_max, Ss_ratio, t, Su_ini, Rp_ini, DRN_tmp, DRN_ini, Ss_ini, E0[t], dtwt, st, i, j, n, kTu_min, kTu_n, dt, dti)
            # fill the output arrays
            Ss[t] = Ss_tmp
            Ro[t]   = Ro_tmp
            Es[t]   = Es_tmp
            Eg[t]   = Eg_tmp
            Tg[t]   = Tg_tmp
            R[t]    = R_tmp
            HEADScorr[t] = HEADS_tmp
            for l in range(nsl):
                Su[t,l]   = Su_tmp[l]
                Su_pc[t,l] = Su_pc_tmp[l]
                Rp[t,l]  = Rp_tmp[l]
                Eu[t,l]  = Eu_tmp[l]
                Tu[t,l]  = Tu_tmp[l]
            ETg[t] = Eg[t] + Tg[t]
            Rn[t]  = R[t] - ETg[t] - DRN_corr
            if (Ss[t]-Ss_ini) > 1E-9:
                dSs[t] = (Ss[t] - Ss_ini)/dt
            # compute the water mass balance (MB) in the unsaturated zone
            Rp_in_MB  = 0.0
            Rp_out_MB = 0.0
            Eu_MB    = 0.0
            Tu_MB    = 0.0
            for l in range(nsl):
                Eu_MB += Eu[t,l]
                Tu_MB += Tu[t,l]
                Rp_out_MB += Rp[t,l]
                if t==0:
                    dSu[t,l] = (Su[t,l]-Su_ini[l])/dt
                    if l<(nsl-1):
                        Rp_in_MB += Rp_ini[l]
                else:
                    dSu[t,l] = (Su[t,l]-Su[t-1,l])/dt
                    if l<(nsl-1):
                        Rp_in_MB += Rp[t-1,l]
                dSu_tot[t] += dSu[t,l]
            ETu_tot[t] = Eu_MB + Tu_MB
            MB[t] = RFe_tot[t] + Ss_ini + dti*Rp_in_MB + DRN_tmp - Ss[t] - dt*(Ro[t] + Es[t] + Eu_MB + Tu_MB + Rp_out_MB + dSu_tot[t])
            # index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iSs':4, 'iRo':5, 'iDRN':6, 'iEs':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'iR':13, 'iRn':14, 'idSs':15, 'iETg':16, 'iETu':17, 'idSu':18, 'iHEADScorr':19}
            results1[t,:] = [RF[t], PET_tot[t], PE_tot[t], RFe_tot[t], Ss[t], Ro[t], DRN_tmp/dt, Es[t], MB[t], INTER_tot[t], E0[t], Eg[t], Tg[t], R[t], Rn[t], dSs[t], ETg[t], ETu_tot[t], dSu_tot[t], HEADScorr[t]/1000.0, -dtwt/1000, Dl[nsl-1]]
            # index_S = {'iEu':0, 'iTu':1,'iSu_pc':2, 'iRp':3, 'idSu':4, 'iSu':5}
            for l in range(nsl):
                results2[t,l,:] = [Eu[t,l], Tu[t,l], Su_pc[t,l], Rp[t,l], dSu[t,l], Su[t,l]]
            dti = dt
            Su_ini = Su[t,:]
            Ss_ini = Ss[t]
            DRN_ini = DRN_tmp
            Dl_ini = Dl
        return results1, results2, Dl

        del RF, PET_tot, PE_tot, RFe_tot, Ss, Ro, DRN, Es, MB, INTER_tot, E0, Eg, Tg, R, Rn, dSs, ETg, ETu_tot, dSu_tot, HEADS, HEADScorr, dtwt
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
        h1=np.zeros([len(R)], dtype=np.float)
        h_tmp = (hi*1000.0 + R[0]/STO -hi*1000.0/RC)
        h1[0]=h_tmp
        for t in range(1,len(R)):
            h_tmp = h1[t-1] + R[t]/STO -h1[t-1]/RC
            h1[t]=h_tmp

        h=np.zeros([len(R)], dtype=np.float)
        for t in range(0,len(R)):
            h_tmp= (h1[t] + h0*1000.0)/1000.0
            h[t]=h_tmp

        return h
        del R, h1, h_tmp, h

#EOF#