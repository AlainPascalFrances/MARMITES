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
It is link to MODFLOW2000 (mf2k) through PEST (Doherty 2010).
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
import pylab

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
                Si          Inicial soil moisture
                Ks          Saturated hydraulic condutivity
                SUSTm     Max. surface storage
            STATE VARIABLES
                RF           Daily rainfall
                PET         Daily evapotranspiration
    OUTPUTS
            RFe              Daily Excess rainfall
            ETu             Daily evapotranspiration
            S               Daily soil moisture
            Rp              Daily percolation
            SUST            Daily ponding
            Qs              Daily runoff

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
                       'clay'            : {'dll':450.0,'y0':0.006,'b':0.0019, 'ext_d':6200.0},
                       }

#####################

    def unsatflux(self, PET, PE, DRN, Sini, Zr_botavg, Dltop_tmp, Dlbot_tmp, Dl_tmp, Dl, Sm, Sfc, Sr, Ks, SUSTm):

        def pond(s_tmp,Dl,Sm,SUSTm):
            '''
            Ponding and surface runoff function
            '''
            countPOND_tmp = 0
            countRunoff_tmp = 0
            if (s_tmp-Sm)>0.000001:
                sust_tmp =Dl*(s_tmp-Sm)
                if SUSTm>0.0:
                    countPOND_tmp = 1
                if sust_tmp-SUSTm >= 0.000001:
                    qs_tmp =(sust_tmp-SUSTm)
                    sust_tmp = SUSTm
                    countRunoff_tmp = 1
                else:
                    qs_tmp =(0.0)
            else:
                sust_tmp =(0.0)
                qs_tmp =(0.0)
            return sust_tmp, qs_tmp, countPOND_tmp, countRunoff_tmp

        ##################

        def perc(s_tmp,Dl,Sm,Sfc,Ks, s_lp1, Dl_sp1, Sm_sp1):
            '''
            Percolation function
            '''
            # Percent. of gravitational water
            Sg=(s_tmp-Sfc)/(Sm-Sfc)
            if s_tmp-Sfc<=0.000001:
                rp_tmp= (0.0)
            elif (s_tmp-Sfc)>0.000001 and (s_tmp-Sm)<=0.000001:
                if (Ks*Sg-Dl*(s_tmp-Sfc))> 0.000001:
                    rp_tmp= Dl*(s_tmp-Sfc)
                else:
                    rp_tmp= Ks*Sg
            elif (s_tmp-Sm)>0.000001:
                if Ks-Dl*(Sm-Sfc)>0.000001:
                    rp_tmp= Dl*(Sm-Sfc)
                else:
                    rp_tmp = Ks
            if rp_tmp-(Sm_sp1-s_lp1)*Dl_sp1>0.000001:
                rp_tmp = (Sm_sp1-s_lp1)*Dl_sp1
            return rp_tmp

        ##################

        def perc1(s_tmp,Dl,Sm,Sfc,Ks, s_lp1, Dl_sp1, Sm_sp1):
            '''
            Percolation function
            # SWAT PERCOLATION pag 150 chap 2:3.2
            '''
            SWe = (s_tmp-Sfc)*Dl
            TTperc = (Sm-Sfc)*Dl/Ks
            if (s_tmp-Sfc)<=0.000001:
                rp_tmp = 0.0
            elif (s_tmp-Sfc)>0.000001:
                rp_tmp = SWe*(1-(pylab.exp(-1/TTperc)))
            return rp_tmp

        ##################

        def evp(s_tmp,pet,Dl,Sm,Sr):
            '''
            Actual evapotranspiration function
            '''
            # Percent. of soil saturation
            Se=(s_tmp-Sr)/(Sm-Sr)
            if s_tmp-Sr<=0.000001:
                evp_tmp= 0.0
            elif s_tmp-Sr>0.000001 and s_tmp-Sm<=0.000001:
                if (pet*Se-(Dl*(s_tmp-Sr)))> 0.000001:
                    evp_tmp= Dl*(s_tmp-Sr)
                else:
                    evp_tmp= pet*Se
            elif (s_tmp-Sm)>0.000001:
                if pet-Dl*(Sm-Sr)>0.000001:
                    evp_tmp= Dl*(Sm-Sr)
                else:
                    evp_tmp= pet
            return evp_tmp

        # MAIN

        if DRN>0.0:
            1



        # SUST and Qs
        PONDtmp = pond(Sini[0]/Dl_tmp[0],Dl_tmp[0],Sm[0],SUSTm)
        SUSTtmp = PONDtmp[0]
        Qstmp = PONDtmp[1]
        countPONDtmp = PONDtmp[2]
        countRunofftmp = PONDtmp[3]
        Sini[0] = Sini[0]-(SUSTtmp+Qstmp)/Dl_tmp[0]

        Rptmp = []
        Tutmp = []
        Eutmp = []
        for l in range(len(Dl_tmp)):
            if Dl_tmp[l] == 0.0: # it means that the soil layer is fully saturated by GW
                Rptmp.append(0.0)
                Eutmp.append(0.0)
                Tutmp.append(0.0)
                Sini[l]=Sm[l]
            else:
                # Rp
                if l==len(Dl_tmp)-1:
                    Rptmp.append(perc1(Sini[l],Dl_tmp[l],Sm[l],Sfc[l],Ks[l],0,1,1))
                else:
                    Rptmp.append(perc1(Sini[l],Dl_tmp[l],Sm[l],Sfc[l],Ks[l], Sini[l+1]/Dl_tmp[l+1],Dl_tmp[l+1], Sm[l+1]))
                Sini[l] = Sini[l]-Rptmp[l]/Dl_tmp[l]
                # Tu
                if l>0:
                    PET = PET - Tutmp[l-1]
                if PET>0.0:
                    if Dlbot_tmp[l] > Zr_botavg :
                        Tutmp.append(evp(Sini[l],PET,Dl_tmp[l],Sm[l],Sr[l]))
                    elif Dltop_tmp[l] > Zr_botavg:
                        PETc = PET*(Dltop_tmp[l]-Zr_botavg)/Dl_tmp[l]
                        Tutmp.append(evp(Sini[l],PETc,Dl_tmp[l],Sm[l],Sr[l]))
                    else:
                        Tutmp.append(0.0)
                else:
                    Tutmp.append(0.0)
                Sini[l]=Sini[l]-Tutmp[l]/Dl_tmp[l]
                # Eu
                if l>0:
                    PE = PE - Eutmp[l-1]
                if PE>0.0:
                    Eutmp.append(evp(Sini[l],PE,Dl_tmp[l],Sm[l],Sr[l]))
                else:
                    Eutmp.append(0.0)
                Sini[l]=Sini[l]-Eutmp[l]/Dl_tmp[l]

        return (SUSTtmp, Qstmp, countPONDtmp, countRunofftmp, Rptmp, Eutmp, Tutmp, Sini)

        ##################

    def satflux(self, PE, dtwt, st):
        '''
        Groundwater evaporation, equation 17 of Shah et al 2007, see ref in the __init__
        '''
        if PE>0.0:
            y0  = self.paramEg[st]['y0']
            b   = self.paramEg[st]['b']
            dll = self.paramEg[st]['dll']
            if dtwt<=dll:
                Egtmp = PE
            elif dtwt>0.0:
                Egtmp = PE*(y0 + pylab.exp(-b*(dtwt-dll)))
            else:
                Egtmp = 0.0
            if Egtmp/PE < 0.5/100.0:
                Egtmp = 0.0
        else:
            Egtmp =0.0

        Tgtmp = 0.0

        return (Egtmp, Tgtmp)

        ##################

    def Sini(self, t, Sini, Si, Dl_tmp, RFe_tot, nsl, Rpi, SUST, E0, SUSTprev, S, Rp, Sm):
        if t == 0:
            # if first time step, use Si and Rpi
            Sini[0] = Si[0] + RFe_tot/Dl_tmp[0]
            if nsl>1:
                for l in range(1,nsl):
                    Sini[l] = Si[l] + Rpi[l-1]/Dl_tmp[l]
            SUSTprev = 0.0
            Estmp = 0.0
        else:
            SUSTprev = SUST
            # soil reservoir water content initialisation
            # test if SUST can infiltrate in soil
            # compute Es from SUST and E0
            if E0-SUSTprev>=0.000001:
                Sini[0] = S[0,t-1] + RFe_tot/Dl_tmp[0]
                Estmp = SUSTprev
            else:
                Sini[0] = S[0,t-1] + (RFe_tot + SUSTprev-E0)/Dl_tmp[0]
                Estmp = E0
            if nsl>1:
                for l in range(1,nsl):
                    if Dl_tmp[l] < 0.000001:
                        Sini[l] = Sm[l]
                    else:
                        Sini[l] = S[l,t-1] + Rp[l-1,t-1]/Dl_tmp[l]
        return Sini, SUSTprev, Estmp

#####################

    def run(self, i, j,
                  nsl, st, slprop, Sm, Sfc, Sr, Si, Rpi, D, Ks, SUSTm,
                  ELEV, HEADS, DRN,
                  RF, E0, PETveg, RFeveg, PEsoil, VEGarea, Zr,
                  perlen, AqType, hdry):

        # Output initialisation
        Ttotal=len(RF)
        # PET for the vegetation patchwork
        PET_tot=np.zeros([len(PETveg[0])], dtype=float)
        # PE for the remaining bare soil
        PE_tot=np.zeros([Ttotal], dtype=float)
        # RFe for the vegetation patchwork
        RFe_tot=np.zeros([len(RFeveg[0])], dtype=float)
        # Surface storage initialisation
        SUST=np.zeros([Ttotal], dtype=np.float)
        # Surface runoff initialisation
        Qs=np.zeros([Ttotal], dtype=np.float)
        #Deep percolation initialisation
        Rp=np.zeros([nsl,Ttotal], dtype=np.float)                             #Soil moisture storage initialisation
        S=np.zeros([nsl,Ttotal], dtype=np.float)                              #Actual evaporation from unsaturated zone initialisation
        Eu=np.zeros([nsl,Ttotal], dtype=np.float)
        #Actual transpiration from unsaturated zone initialisation
        Tu=np.zeros([nsl,Ttotal], dtype=np.float)
        #Actual evapotranspiration from surface initialisation
        Es=np.zeros([Ttotal], dtype=np.float)
        #Flood frequency (saturated overland flow)
        countFLOOD=np.zeros([Ttotal], dtype=int)
        #soil partial saturation by GW frequency
        countSATpart=np.zeros([Ttotal], dtype=int)
        #MASS BALANCE
        MB=np.zeros(Ttotal, dtype=float)
        # hortonian overland flow
        countPOND=np.zeros([Ttotal], dtype=int)
        countRunoff=np.zeros([Ttotal], dtype=int)
        INTER=np.zeros([Ttotal], dtype=np.float)
        Eg=np.zeros([Ttotal], dtype=np.float)
        Tg=np.zeros([Ttotal], dtype=np.float)

        if D < 0.05: #correct soil thickness less than 5cm
            D=0.05
        D = D*1000
        ELEV = ELEV*1000
        # elevation of bottom of soil reservoir
        Dbot = ELEV
        # topography elevation
        Dtop = ELEV + D

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

        Zr_tmp = 0.0
        for z in range(len(Zr)):
            Zr_tmp = Zr_tmp + Zr[z]*VEGarea[z]/100
        Zr_botavg = 1000*Zr_tmp/len(Zr)

        # PROCESSING THE WHOLE DATA SET
        for t in range(int(perlen)):    # t: current time step

            # Preprocessing of PET/PE/INTER and RFe
            # for different vegetation
            SOILarea = 100
            for v in range(len(PETveg)):
                if VEGarea[v]<>self.hnoflo:
                    PET_tot[t] = PET_tot[t]+PETveg[v,t]*VEGarea[v]/100
                    RFe_tot[t] = RFe_tot[t]+RFeveg[v,t]*VEGarea[v]/100
                    SOILarea = SOILarea - VEGarea[v]
            RFe_tot[t] = RFe_tot[t] + RF[t]*SOILarea/100
            INTER[t] = RF[t] - RFe_tot[t]
            PE_tot[t] = PEsoil[t]*SOILarea/100

            # handle drycell
            if np.abs(HEADS[t]-hdry) > 1E+22:
                HEADStmp = Dbot - 1000
            else:
                HEADStmp=HEADS[t]*1000
            dtwt = Dtop-HEADStmp

            Sini = np.zeros([nsl], dtype=float)
            SUSTprev = []

            # reservoir thickness re-initialization
            Dltop_tmp = Dltop*1.0
            Dlbot_tmp = Dlbot*1.0
            Dl_tmp    = Dl*1.0

            if AqType==1:
            # AQUIFER UNCONFINED, water table can rise in the soil and above surface

                # heads below soil bottom
                if HEADStmp-Dbot<=0.000001:
                    # bottom boundary
                    Dlbot_tmp[nsl-1] = HEADStmp
                    Dl_tmp[nsl-1] = Dltop_tmp[nsl-1] - HEADStmp

                    Sini, SUSTprev, Estmp = self.Sini(t, Sini, Si, Dl_tmp, RFe_tot[t], nsl, Rpi, SUST[t-1], E0[t], SUSTprev, S, Rp, Sm)

                    SUSTtmp, Qstmp, countPONDtmp, countRunofftmp, Rptmp, Eutmp, Tutmp, Stmp = self.unsatflux(PET_tot[t], PE_tot[t], DRN[t], Sini, Zr_botavg, Dltop, Dlbot_tmp, Dl_tmp, Dl, Sm,Sfc,Sr,Ks,SUSTm)
                    Egtmp, Tgtmp = self.satflux(PET_tot[t], dtwt, st)

                # heads above soil bottom and below surface
                elif HEADStmp-Dbot>0.000001 and Dtop-HEADStmp>=0.000001:
                    countSATpart[t] = 1
                    for l in range(nsl-1):
                        if (HEADStmp-Dlbot_tmp[l]) > 0.000001:
                            if (HEADStmp-Dltop_tmp[l]) > 0.000001:
                                Dl_tmp[l] = 0.0
                                Dlbot_tmp[l] = 0.0
                                Dltop_tmp[l] = 0.0
                            else:
                                Dl_tmp[l] = Dltop_tmp[l] - HEADStmp
                                Dlbot_tmp[l] = HEADStmp

                    Sini, SUSTprev, Estmp = self.Sini(t, Sini, Si, Dl_tmp, RFe_tot[t], nsl, Rpi, SUST[t-1], E0[t], SUSTprev, S, Rp, Sm)

                    SUSTtmp, Qstmp, countPONDtmp, countRunofftmp, Rptmp, Eutmp, Tutmp, Stmp = self.unsatflux(PET_tot[t], PE_tot[t], DRN[t], Sini, Zr_botavg, Dltop, Dlbot_tmp, Dl_tmp, Dl, Sm,Sfc,Sr,Ks,SUSTm)
                    Egtmp, Tgtmp = self.satflux(PET_tot[t], dtwt, st)

                # heads above surface
                else:
                    Sini, SUSTprev, Estmp = self.Sini(t, Sini, Si, Dl_tmp, RFe_tot[t], nsl, Rpi, SUST[t-1], E0[t], SUSTprev, S, Rp, Sm)
                    countFLOOD[t] = 1
                    countPONDtmp = 0
                    countRunofftmp = 0
                    # TOBE FIXED, MB problem between atm, unsta and GW
                    SUSTtmp = HEADStmp-ELEV+SUSTprev-E0[t]
                    if SUSTtmp > SUSTm:
                        Qstmp =(SUSTtmp-SUSTm)
                        SUSTtmp = SUSTm
                    else:
                        Qstmp =0.0
                    Estmp = E0[t]
                    Eutmp=[]
                    Tutmp=[]
                    Rptmp=[]
                    Stmp=[]
                    for l in range(nsl):
                        Eutmp.append(0.0)
                        Tutmp.append(0.0)
                        Rptmp.append(0.0)
                        Stmp.append(Sm[l])
                    Egtmp, Tgtmp = self.satflux(PET_tot[t], dtwt, st)
            else:
            # TODO to be confirmed, not correct currently
            # AQUIFER CONFINED, water table cannot rise in the soil or above topography
                SUSTtmp, Qstmp, countPONDtmp, countRunofftmp, Rptmp, Eutmp, Tutmp, Stmp, Egtmp, Tgtmp = self.partition(PET_tot[t], PE_tot[t], DRN[t], Sini, Zr_botavg, Dltop, Dlbot_tmp, Dl_tmp, Sm,Sfc,Sr,Ks,SUSTm, st, dtwt)

            # fill the table and compute water balance
            SUST[t]=SUSTtmp
            Qs[t]=Qstmp
            Es[t]=Estmp
            Eg[t]=Egtmp
            Tg[t]=Tgtmp
            countPOND[t]=countPONDtmp
            countRunoff[t]=countRunofftmp
            for l in range(nsl):
                S[l,t]=Stmp[l]
                Rp[l,t]=Rptmp[l]
                Eu[l,t]=Eutmp[l]
                Tu[l,t]=Tutmp[l]

            # compute the water mass balance (MB) in the unsaturated zone
            if countFLOOD[t]==0:
                RpMB=0.0
                EuMB=0.0
                TuMB=0.0
                SMB=0.0
                for l in range(len(Dl)):
                    EuMB = EuMB + Eu[l,t]
                    TuMB = TuMB + Tu[l,t]
                    if t==0:
                        SMB=SMB+(S[l,t]-Si[l])*Dl_tmp[l]
                        if l == 0:
                            RpMB  = RpMB + Rp[l,t]
                        else:
                            RpMB  = RpMB + Rp[l,t] - Rpi[l-1]
                    else:
                        SMB=SMB+(S[l,t]-S[l,t-1])*Dl_tmp[l]
                        if l == 0:
                            RpMB  = RpMB + Rp[l,t]
                        else:
                            RpMB  = RpMB + Rp[l,t] - Rp[l-1,t-1]
                MB[t]=RF[t]+SUSTprev-INTER[t]-Qs[t]-SUST[t]-Es[t]-RpMB-EuMB-TuMB-SMB
            else:
                MB[t]=HEADStmp-ELEV-SUST[t]-Qs[t]+SUSTprev

#index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iSUST':4, 'iQs':5, 'iFLOOD':6, 'iEs':7, 'iSATpart':8, 'iMB':9, 'iINTER':10, 'iPOND':11, 'iRunoff':12, 'iE0':13, 'iEg':14, 'iTg':15, 'iR':16, 'iRn':17}
        results1 = [RF, PET_tot, PE_tot, RFe_tot, SUST, Qs, countFLOOD, Es, countSATpart, MB, INTER, countPOND, countRunoff, E0, Eg, Tg]
        results2 = [Eu, Tu, S, Rp]

        results1 = np.asarray(results1)
        results2 = np.asarray(results2)

        return results1, results2

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
        h_tmp = (hi*1000 + R[0]/STO -hi*1000/RC)
        h1[0]=h_tmp
        for t in range(1,len(R)):
            h_tmp = h1[t-1] + R[t]/STO -h1[t-1]/RC
            h1[t]=h_tmp

        h=np.zeros([len(R)], dtype=np.float)
        for t in range(0,len(R)):
            h_tmp= (h1[t] + h0*1000)/1000
            h[t]=h_tmp

        return h

#######################################################

class PROCESS:
    def __init__(self, MM_ws, MF_ws, nrow, ncol, xllcorner, yllcorner, cellsizeMF, perlen, hnoflo):
        self.MM_ws = MM_ws
        self.MF_ws = MF_ws
        self.nrow= nrow
        self.ncol= ncol
        self.xllcorner= xllcorner
        self.yllcorner=yllcorner
        self.cellsizeMF=cellsizeMF
        self.perlen=perlen
        self.hnoflo=hnoflo

    def inputEsriAscii(self, grid_fn, datatype):

        grid_fn=os.path.join(self.MM_ws,grid_fn)

        grid_out=np.zeros([self.nrow,self.ncol], dtype = datatype)

        grid_out=self.convASCIIraster2array(grid_fn,grid_out)
        return grid_out

    def convASCIIraster2array(self, filenameIN, arrayOUT):
        '''
        Read ESRI/MODFLOW ASCII file and convert to numpy array
        '''

        # Load the grid files
        if os.path.exists(filenameIN):
            fin = open(filenameIN, 'r')
        else:
            raise ValueError, "The file %s doesn't exist!!!" % filenameIN
    #        fout = open(outfile, 'w')

        # Read the header
        line = fin.readline().split()
        ncol_tmp = int(line[1])

        line = fin.readline().split()
        nrow_tmp = int(line[1])

        line = fin.readline().split()
        xllcorner_tmp = float(line[1])

        line = fin.readline().split()
        yllcorner_tmp = float(line[1])

        line = fin.readline().split()
        cellsizeEsriAscii = float(line[1])

        line = fin.readline().split()
        NODATA_value = float(line[1])

        # Process the file
#        print "\nConverting %s to np.array..." % (filenameIN)
        for row in range(nrow_tmp):
        #   if (row % 100) == 0: print ".",
            for col in range(ncol_tmp):
                # Get new data if necessary
                if col == 0: line = fin.readline().split()
                if line[col] == NODATA_value:
                    arrayOUT[row,col] = self.hnoflo
                else:
                    arrayOUT[row,col] = line[col]

        # verify grid consistency between MODFLOW and ESRI ASCII
        if arrayOUT.shape[0] != self.nrow or arrayOUT.shape[1] != self.ncol or self.cellsizeMF != cellsizeEsriAscii:
            raise BaseException, '\nError in consistency between the MODFLOW grid and the input gridof the file %s.\nCheck the cell size and the number of rows, columns and cellsize' % filenameIN

        fin.close()

        return arrayOUT

    def inputTS(self, NMETEO, NVEG, NSOIL,
                inputDate_fn, inputZON_TS_RF_fn,
                inputZON_TS_PET_fn, inputZON_TS_RFe_fn,
                inputZON_TS_PE_fn, inputZON_TS_E0_fn
                ):   #IRR_fn

        # READ date of input files (RF and PET)
        inputDate_fn=os.path.join(self.MM_ws, inputDate_fn)
        if os.path.exists(inputDate_fn):
            inputDate=np.loadtxt(inputDate_fn, dtype = str)
            inputDate = inputDate[:,0]
            inputDate = pylab.datestr2num(inputDate)
            for i in range(1,len(inputDate)):
                #__________________Check date consistency________________#
                difDay=inputDate[i]-inputDate[i-1]
                if (difDay !=1.0):
                    print 'DifDay = ' + str(difDay)
                    raise ValueError, 'The dates of the input data (RF and PET) are not sequencial, check your daily time step!\nError in date %s ' % str(inputDate[i])
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % inputDate_fn
        if len(inputDate)<>sum(self.perlen):
            raise ValueError, 'The number of time steps in MF (%i) is not the same as the number of days (%i) of the input data (RF and PET).\n' % (int(sum(self.perlen)), int(len(inputDate)))

        # READ input ESRI ASCII rasters vegetation
        gridVEGarea_fn=[]
        for v in range(NVEG):
            gridVEGarea_fn.append(os.path.join(self.MM_ws,'inputVEG' + str(v+1)+'area.asc'))
        gridVEGarea=np.zeros([NVEG,self.nrow,self.ncol], dtype=float)
        for v in range(NVEG):
            gridtmp=np.zeros([self.nrow,self.ncol], dtype=float)
            gridVEGarea[v,:,:]=self.convASCIIraster2array(gridVEGarea_fn[v],gridtmp)

        # READ RF for each zone
        RF_fn=os.path.join(self.MM_ws, inputZON_TS_RF_fn)
        if os.path.exists(RF_fn):
            RF = np.loadtxt(RF_fn)
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % RF_fn
        RFzonesTS=np.zeros([NMETEO,sum(self.perlen)], dtype=float)
        for i in range(NMETEO):
            for t in range(int(sum(self.perlen))):
                RFzonesTS[i,t]=RF[i*sum(self.perlen)+t]

        E0_fn=os.path.join(self.MM_ws, inputZON_TS_E0_fn)
        if os.path.exists(E0_fn):
            E0 = np.loadtxt(E0_fn)
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % E0_fn
        E0zonesTS=np.zeros([NMETEO,sum(self.perlen)], dtype=float)
        for i in range(NMETEO):
            for t in range(int(sum(self.perlen))):
                E0zonesTS[i,t]=E0[i*sum(self.perlen)+t]

        ### READ IRR for each zone
        ##IRR_fn=os.path.join(self.MM_ws,IRR_fn)
        ##if os.path.exists(IRR_fn):
        ##    IRR = np.loadtxt(IRR_fn)
        ##else:
        ##    raise ValueError, "\nThe file %s doesn't exist!!!" % IRR_fn
        ##IRRzonesnumb=int(IRR[0])
        ##IRRzones=np.zeros([IRRzonesnumb,sum(dis.self.perlen)], dtype=float)
        ##for i in range(IRRzonesnumb):
        ##    for j in range(sum(dis.self.perlen)):
        ##        IRRzones[i,j]=IRR[i*sum(dis.self.perlen)+j+1]


        # READ PET for each zone and each vegetation
        PET_fn = os.path.join(self.MM_ws, inputZON_TS_PET_fn)
        if os.path.exists(PET_fn):
            PETtmp = np.loadtxt(PET_fn)
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % PET_fn
        PETvegzonesTS=np.zeros([NMETEO,NVEG,int(sum(self.perlen))], dtype=float)
        for z in range(NMETEO):
            for v in range(NVEG):
                for t in range(int(sum(self.perlen))):
                    PETvegzonesTS[z,v,t]=PETtmp[z*v*int(sum(self.perlen))+t]
                    #structure is [number of zones, number of vegetation type, time]
        PETvegzonesTS=np.asarray(PETvegzonesTS)

        # READ RFe for each zone and each vegetation
        RFe_fn = os.path.join(self.MM_ws, inputZON_TS_RFe_fn)
        if os.path.exists(RFe_fn):
            RFetmp = np.loadtxt(RFe_fn)
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % RFe_fn
        RFevegzonesTS=np.zeros([NMETEO,NVEG,int(sum(self.perlen))], dtype=float)
        for z in range(NMETEO):
            for v in range(NVEG):
                for t in range(int(sum(self.perlen))):
                    RFevegzonesTS[z,v,t]=RFetmp[z*v*int(sum(self.perlen))+t]
                    #structure is [number of zones, number of vegetation type, time]
        RFevegzonesTS=np.asarray(RFevegzonesTS)


        # READ PE for each zone and each soil
        PE_fn = os.path.join(self.MM_ws, inputZON_TS_PE_fn)
        if os.path.exists(PE_fn):
            PEtmp = np.loadtxt(PE_fn)
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % PE_fn
        PEsoilzonesTS=np.zeros([NMETEO,NSOIL,int(sum(self.perlen))], dtype=float)
        for z in range(NMETEO):
            for v in range(NSOIL):
                for t in range(int(sum(self.perlen))):
                    PEsoilzonesTS[z,v,t]=PEtmp[z*v*int(sum(self.perlen))+t]
                    #structure is [number of zones, number of vegetation type, time]
        PEsoilzonesTS=np.asarray(PEsoilzonesTS)

        return gridVEGarea, RFzonesTS, E0zonesTS, PETvegzonesTS, RFevegzonesTS, PEsoilzonesTS, inputDate


    def inputSoilParam(self, SOILparam_fn, NSOIL):

        # Soils parameter initialisation
        nam_soil=[]
        nsl=[]
        st=[]
        slprop=[]
        Sm=[]
        Sfc=[]
        Sr=[]
        Si=[]
        Ks=[]

        # soil parameters file
        SOILparam_fn=os.path.join(self.MM_ws,SOILparam_fn)
        if os.path.exists(SOILparam_fn):
            fin = open(SOILparam_fn, 'r')
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % SOILparam_fn
        line = fin.readline().split()
        delimChar = line[0]
        inputFile = []
        try:
            for line in fin:
                line_tmp = line.split(delimChar)
                if not line_tmp == []:
                    if (not line_tmp[0] == '') and (not line_tmp[0] == '\n'):
                        inputFile.append(line_tmp[0])
                else:
                    raise NameError('InputFileFormat')
        except NameError:
            raise ValueError, 'Error in the input file, check format!\n'

        SOILzones=int(int(inputFile[0]))
        if SOILzones>NSOIL:
            print '\nWARNING:\n' + str(SOILzones) + ' soil parameters groups in file [' + SOILparam_fn + ']\n Only ' + str(NSOIL) + ' PE time serie(s) found.'
        # trick to initialise the reading position in the next loop
        nsl.append(0)
        for i in range(SOILzones):
            nsl.append(int(inputFile[i+1]))

        # soil parameter definition for each soil type
        nslst = SOILzones+1
        nam_soil=[]
        st=[]
        for z in range(SOILzones):
            nam_soil.append(inputFile[nslst].strip())
            nslst = nslst + 1
            st.append(inputFile[nslst].strip())
            nslst = nslst + 1
            slprop.append([])
            Sm.append([])
            Sfc.append([])
            Sr.append([])
            Si.append([])
            Ks.append([])
            for ns in range(nsl[z+1]):
                slprop[z].append(float(inputFile[nslst]))
                nslst = nslst + 1
                Sm[z].append(float(inputFile[nslst]))
                nslst = nslst + 1
                Sfc[z].append(float(inputFile[nslst]))
                nslst = nslst + 1
                Sr[z].append(float(inputFile[nslst]))
                nslst = nslst + 1
                Si[z].append(float(inputFile[nslst]))
                nslst = nslst + 1
                Ks[z].append(float(inputFile[nslst]))
                nslst = nslst + 1
            if sum(slprop[z][0:nsl[z+1]-1])>1:
                raise ValueError, '\nERROR!\nThe sum of the soil layers proportion of %s is >1!\nCorrect your soil data input!\nEXISTING MARMITES\n' % nam_soil[z]
            if sum(slprop[z][0:nsl[z+1]-1])<1:
                raise ValueError, '\nERROR!\nThe sum of the soil layers proportion of %s is <1!\nCorrect your soil data input!\nEXISTING MARMITES\n' % nam_soil[z]

        return nsl[1:len(nsl)], nam_soil, st, slprop, Sm, Sfc, Sr, Si, Ks


    def inputObs(self, inputObs_fn, inputDate):
        '''
        observations cells for soil moisture and heads (will also compute SATFLOW)
        '''

        # read coordinates and SATFLOW parameters
        filenameIN = os.path.join(self.MM_ws,inputObs_fn)
        if os.path.exists(filenameIN):
            fin = open(filenameIN, 'r')
            lines = fin.readlines()
            fin.close()
        else:
            raise BaseException, "The file %s doesn't exist!!!" % filenameIN

        # define a dictionnary of observatins,  format is: Name (key) x y i j hi h0 RC STO
        obs = {}
        for i in range(1,len(lines)):
            line = lines[i].split()
            name = line[0]
            x = float(line[1])
            y = float(line[2])
            try:
                hi  = float(line[3])
                h0  = float(line[4])
                RC  = float(line[5])
                STO = float(line[6])
            except:
                hi = h0 = RC = STO =  self.hnoflo
            # verify if coordinates are inside MODFLOW grid
            if x < self.xllcorner or \
               x > (self.xllcorner+self.ncol*self.cellsizeMF) or \
               y < self.yllcorner or \
               y > (self.yllcorner+self.nrow*self.cellsizeMF):
                   raise BaseException, 'The coordinates of the observation point %s are not inside the MODFLOW grid' % name
            # compute the cordinates in the MODFLOW grid
            #TODO use the PEST utilities for space extrapolation
            i = self.nrow - np.ceil((y-self.yllcorner)/self.cellsizeMF)
            j = np.ceil((x-self.xllcorner)/self.cellsizeMF) - 1
            obs[name] = {'x':x,'y':y,'i': i, 'j': j, 'hi':hi, 'h0':h0, 'RC':RC, 'STO':STO}

        outpathname=[]
        #  read obs time series
        for o in range(len(obs.keys())):
            outpathname.append(os.path.join(self.MM_ws,'output'+obs.keys()[o]+'.txt'))
        obs_h=[]
        obs_sm=[]
        for o in range(len(obs.keys())):
            obsh_fn=os.path.join(self.MM_ws,'inputObsHEADS_' + obs.keys()[o] +'.txt')
            obssm_fn=os.path.join(self.MM_ws,'inputObsSM_' + obs.keys()[o] + '.txt')
            obs_h.append(self.verifObs(inputDate, obsh_fn))
            obs_sm.append(self.verifObs(inputDate, obssm_fn))

        return obs, outpathname, obs_h, obs_sm


    def verifObs(self, inputDate, filename):
        '''
        Import and process data and parameters
        '''

        try:
            #_________________Open file___________________________#
            obsOutput=np.zeros([len(inputDate)], dtype=float)
            if os.path.exists(filename):
                obsData=np.loadtxt(filename, dtype = str)
                obsDate = pylab.datestr2num(obsData[:,0])
                obsValue = obsData[:,1].astype(float)
                if obsDate[0]<inputDate[0]:
                    print '\nWARNING, Obs data starting before RF and PET,\n these obs data will not be plotted correctly'
                if len(inputDate)<len(obsDate):
                    print '\nWARNING, there is more Obs data than RF and PET data,\n these obs data will not be plotted correctly'
                j=0
                for i in range(len(inputDate)):
                    if j<len(obsDate):
                        if inputDate[i]==obsDate[j]:
                            obsOutput[i]=obsValue[j]
                            j=j+1
                        else:
                            obsOutput[i]=self.hnoflo
                    else:
                        obsOutput[i]=self.hnoflo
            else:
                for i in range(len(inputDate)):
                    obsOutput[i]=self.hnoflo

            obsOutput = np.asarray(obsOutput)

        except (ValueError, TypeError, KeyboardInterrupt, IOError), e:
            obsOutput = 0
            print e

        return obsOutput

    def outputEAgrd(self, outFile_fn, outFolder = []):

        if outFolder == []:
            outFolder = self.MM_ws

        outFile=open(os.path.join(outFolder, outFile_fn), 'w')

        outFile=self.writeHeaderESRIraster(outFile)

        return outFile

    def writeHeaderESRIraster(self, file_asc):
        '''
        Write ESRI ASCII raster header
        '''
        file_asc.write('ncols  ' + str(self.ncol)+'\n' +
                       'nrows  '  + str(self.nrow)+'\n' +
                       'xllcorner  '+ str(self.xllcorner)+'\n' +
                        'yllcorner  '+ str(self.yllcorner)+'\n' +
                        'cellsize   '+ str(self.cellsizeMF)+'\n' +
                        'NODATA_value  '+ str(self.hnoflo)+'\n')
        return file_asc

    def ExportResults(self, i, j, inputDate, _nslmax, results, index, results_S, index_S, h_satflow, heads_MF, obs_h, obs_S, outFileExport, outPESTheads, outPESTsm, obsname):
        """
        Export the processed data in a txt file
        INPUTS:      output fluxes time series and date
        OUTPUT:      output.txt
        """
        for t in range(len(inputDate)):
            year='%4d'%pylab.num2date(inputDate[t]).year
            month='%02d'%pylab.num2date(inputDate[t]).month
            day='%02d'%pylab.num2date(inputDate[t]).day
            date=(day+"/"+month+"/"+year)
            if obs_h[t]!=self.hnoflo:
                obs_h_tmp = str(obs_h[t])
                outPESTheads.write(obsname.ljust(10,' ')+ date.ljust(14,' ')+ '00:00:00        '+ str(heads_MF[t])+ '    \n')
            else:
                obs_h_tmp = str(self.hnoflo)
            # TODO export for all soil layers (currently only the SM computed for the first layer is exported)
            if obs_S[t]!=self.hnoflo:
                obs_S_tmp = str(obs_S[t])
                outPESTsm.write((obsname+'SM').ljust(10,' ')+ date.ljust(14,' ')+ '00:00:00        '+ str(results_S[i,j,index_S.get('iS'),0,t]) + '    \n')
            else:
                obs_S_tmp = str(self.hnoflo)
            # header='Date,RF,E0,PET,PE,RFe,Inter,'+Eu_str+Tu_str+'Eg,Tg,Es,'+S_str+'SUST,SUSTcount,Qs,Qscount,'+Rp_str+'R,hSATFLOW,hMF,hmeas,Smeas,SSATpart,FLOODcount,MB\n'
            Sout = ''
            Rpout=''
            Euout=''
            Tuout=''
            for l in range(_nslmax):
                Sout = Sout + str(results_S[i,j,index_S.get('iS'),l,t]) + ','
                Rpout = Rpout + str(results_S[i,j,index_S.get('iRp'),l,t]) + ','
                Euout = Euout + str(results_S[i,j,index_S.get('iEu'),l,t]) + ','
                Tuout = Tuout + str(results_S[i,j,index_S.get('iTu'),l,t]) + ','
            out_date = pylab.num2date(inputDate[t]).isoformat()[:10]
            out1 = '%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,' % (results[i,j,index.get('iRF'),t], results[i,j,index.get('iE0'),t],results[i,j,index.get('iPET'),t],results[i,j,index.get('iPE'),t],results[i,j,index.get('iRFe'),t],results[i,j,index.get('iINTER'),t])
            out1a = '%.6f,%.6f,' % (results[i,j,index.get('iEg'),t], results[i,j,index.get('iTg'),t])
            out2 = '%.6f,' % (results[i,j,index.get('iEs'),t])
            out3 = '%.6f,%4d,%.6f,%4d,' % (results[i,j,index.get('iSUST'),t],results[i,j,index.get('iPOND'),t],results[i,j,index.get('iQs'),t],results[i,j,index.get('iRunoff'),t])
            out4 = '%.6f,%.6f,%.6f,%.6f,%.6f,%4d,%4d,%6f' % (results[i,j,index.get('iR'),t], h_satflow[t],heads_MF[t],obs_h[t], obs_S[t],results[i,j,index.get('iSATpart'),t],results[i,j,index.get('iFLOOD'),t],results[i,j,index.get('iMB'),t])
            out_line =  out_date, ',', out1, Euout, Tuout, out1a, out2, Sout, out3, Rpout, out4, '\n'
            for l in out_line:
                outFileExport.write(l)

# EOF