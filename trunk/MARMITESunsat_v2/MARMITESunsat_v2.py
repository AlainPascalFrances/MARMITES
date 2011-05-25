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
                PONDm       Max. ponding capacity
            STATE VARIABLES
                RF           Daily rainfall
                PET         Daily evapotranspiration
    OUTPUTS
            RFe              Daily Excess rainfall
            ETu             Daily evapotranspiration
            S               Daily soil moisture
            Rp              Daily percolation
            POND            Daily ponding
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

    def flux(self, RFe, PET, PE, S, Sini, Rp, Zr_elev, VEGarea, HEADStmp, Dltop, Dlbot, Dl, nsl, Sm, Sfc, Sr, Ks, PONDm, PONDratio, t, Si, Rpi, DRN, DRNi, PONDi, E0, dtwt, st, i, j, n, k_Tu_slp, k_Tu_inter):

        def surfwater(s_tmp,Vm,PONDm, E0, i, j, n):
            '''
            Ponding and surface runoff function
            '''
            if (s_tmp-Vm) > 1E-6:
                sust_tmp = s_tmp-Vm
                if sust_tmp-PONDm >= 1E-6:
                    Ro_tmp = sust_tmp-PONDm
                    sust_tmp = PONDm
                else:
                    Ro_tmp = 0.0
                if sust_tmp - E0 >= 1E-6:
                    Estmp = E0
                    sust_tmp = sust_tmp - E0
                else:
                    Estmp = sust_tmp
                    sust_tmp = 0.0
            else:
                sust_tmp = 0.0
                Ro_tmp = 0.0
                Estmp = 0.0
            return sust_tmp, Ro_tmp, Estmp

        ##################

        def perc(s_tmp,Sm,Sfc,Ks, s_lp1, Sm_sp1, i, j, n):
            '''
            Percolation function
            # SWAT PERCOLATION pag 150 chap 2:3.2
            '''
            SWe = s_tmp-Sfc
            TTperc = (Sm-Sfc)/Ks
            if (s_tmp-Sfc)<=1E-6:
                rp_tmp = 0.0
            elif (s_tmp-Sfc)>1E-6:
                rp_tmp = SWe*(1-(pylab.exp(-1/TTperc)))
                # verify the vol. available in deeper soil layer
                if rp_tmp-(Sm_sp1 - s_lp1)>1E-6:
                    rp_tmp = Sm_sp1 - s_lp1
            return rp_tmp

        ##################

        def evp(s_tmp,pet,Sm,Sr, i, j, n):
            '''
            Actual evapotranspiration function
            '''
            # Percent. of soil saturation
            Se=(s_tmp - Sr)/(Sm - Sr)
            if s_tmp - Sr <= 1E-6:
                evp_tmp= 0.0
            elif s_tmp - Sr > 1E-6 and s_tmp - Sm <= 1E-6:
                if (pet*Se - (s_tmp - Sr))> 1E-6:
                    evp_tmp= s_tmp - Sr
                else:
                    evp_tmp= pet*Se
            elif (s_tmp - Sm) > 1E-6:
                if pet-(Sm - Sr) > 1E-6:
                    evp_tmp= Sm - Sr
                else:
                    evp_tmp= pet
            return evp_tmp

        ##################

        # S from previous time step
        # if first time step, use Si
        for l in range(nsl):
            if t == 0:
                Sini[l] = Si[l]
            else:
                Sini[l] = S[t-1,l]

        DRNcorr = 0.0
        if DRN == 0.0:
            if DRNi>0.0:
                DRNcorr = (Sfc[nsl-1]-Sr[nsl-1])*Dl[nsl-1]
                Sini[nsl-1] = Sfc[nsl-1]*Dl[nsl-1]
                Si[nsl-1] = Sini[nsl-1]
        else:
            # GW seepage
            llst = range(nsl-1)
            llst.reverse()
            DRN += Sini[nsl-1]-Sr[nsl-1]*Dl[nsl-1]
            for l in llst:
                if DRN - (Sm[l]*Dl[l]-Sini[l]) > 1E-6:
                    DRN -= (Sm[l]*Dl[l]-Sini[l])
                    Sini[l] = Sm[l]*Dl[l]
                else:
                    Sini[l] += DRN
                    DRN = 0.0
            Sini[nsl-1] = None

        # percolation from previous time step
        # if first time step, use Rpi
        for l in range(nsl-1):
            if t == 0:
                Rpprev = Rpi[l]
            else:
                Rpprev = Rp[t-1,l]
            if Dl[l+1]>0.0:
                if Rpprev - (Sm[l+1]*Dl[l+1]-Sini[l+1]) > 1E-6:
                    Rpprev -= (Sm[l+1]*Dl[l+1]-Sini[l+1])
                    Sini[l+1] = Sm[l+1]*Dl[l+1]
                    Sini[l] += Rpprev
                else:
                    Sini[l+1] += Rpprev
            else:
                Sini[l] += Rpprev

        # infiltration
        Sini[0] += (RFe + PONDi + DRN)

        # POND and Ro
        surfwatertmp = surfwater(Sini[0],Sm[0]*Dl[0],PONDm, E0*PONDratio, i, j, n)
        PONDtmp = surfwatertmp[0]
        Rotmp = surfwatertmp[1]
        Estmp = surfwatertmp[2]
        Sini[0] -= (PONDtmp+Rotmp+Estmp)

        Rptmp = np.zeros([nsl])
        Rtmp = 0.0
        TutmpZr = np.zeros([nsl,len(Zr_elev)])
        Tutmp = np.zeros([nsl])
        Eutmp = np.zeros([nsl])
        Spc = np.zeros([nsl])
        # soil layers and vadose zone layer
        for l in range(nsl):
            # Rp
            if l < (nsl-1):
                if Dl[l+1]>0.0:
                    Rptmp[l] = perc(Sini[l],Sm[l]*Dl[l],Sfc[l]*Dl[l],Ks[l], Sini[l+1], Sm[l+1]*Dl[l+1], i, j, n)
                    Sini[l] -= Rptmp[l]
            else:
                if Dl[l]>0.0:
                    Rptmp[l] = perc(Sini[l],Sm[l]*Dl[l],Sfc[l]*Dl[l],Ks[l],0.0,1.0E6, i, j, n)
                    Rtmp = Rptmp[l]
                    Sini[l] -= Rptmp[l]
            if Dl[l]>0.0:
                # Eu
                if PE>0.0:
                    Eutmp[l] = evp(Sini[l],PE,Sm[l]*Dl[l],Sr[l]*Dl[l], i, j, n)
                Sini[l] -= Eutmp[l]
                PE -= Eutmp[l]
                # Tu
                for z in range(len(Zr_elev)):
                    if PET[z]>0.0:
                        if Dlbot[l] > Zr_elev[z]:
                            TutmpZr[l,z] = evp(Sini[l],PET[z],Sm[l]*Dl[l],Sr[l]*Dl[l], i, j, n)
                        elif Dltop[l] > Zr_elev[z]:
                            PETc = PET[z]*(Dltop[l]-Zr_elev[z])/Dl[l]
                            TutmpZr[l,z] = evp(Sini[l],PETc,Sm[l]*Dl[l],Sr[l]*Dl[l], i, j, n)
                        PET[z] -= TutmpZr[l,z]
                        Tutmp[l] += TutmpZr[l,z]*VEGarea[z]/100
                Sini[l] -= Tutmp[l]

        for l in range(nsl):
            if Dl[l]>0.0:
                Spc[l] = Sini[l]/Dl[l]
            else:
                Spc[l] = Sm[l]
                Sini[l] = Sm[l]*Dl[l]

        # GW evaporation, equation 17 of Shah et al 2007, see ref in the __init__
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

        # Groundwater transpiration
        TgtmpZr = np.zeros([len(Zr_elev)])
        Tgtmp = 0.0
        for z in range(len(Zr_elev)):
            if HEADStmp > Zr_elev[z]:
                # TODO correct average soil moisture in the case vadose layer is filled with GW
                k_Tu = k_Tu_slp[z]*sum(Spc[0:nsl-1])/(nsl-1) + k_Tu_inter[z]
                TgtmpZr[z] = sum(TutmpZr[0:nsl-1,z]/(nsl-1))*(1/k_Tu-1)
                if TgtmpZr[z]>PET[z]:
                    TgtmpZr[z] = PET[z]
                Tgtmp += TgtmpZr[z]*VEGarea[z]/100

        return Estmp, PONDtmp, Rotmp, Rptmp, Eutmp, Tutmp, Sini, Spc, Si, Rtmp, Egtmp, Tgtmp, DRNcorr
        del Estmp, PONDtmp, Rotmp, Rptmp, Eutmp, Tutmp, Sini, Spc, Si, Rtmp, Egtmp, Tgtmp, DRNcorr

#####################

    def run(self, i, j, n,
                  nsl, st, slprop, Sm, Sfc, Sr, Si, PONDi, Rpi, D, Ks, PONDm, PONDratio,
                  ELEV, HEADS, DRN, DRNi,
                  RF, E0, PETveg, RFeveg, PEsoil, VEGarea, Zr,
                  nstp, hdry,
                  k_Tu_slp, k_Tu_inter):

        # Output initialisation
        Ttotal=len(RF)
        # PET for the vegetation patchwork
        PET_tot=np.zeros([len(PETveg[0])], dtype=float)
        # PE for the remaining bare soil
        PE_tot=np.zeros([Ttotal], dtype=float)
        # RFe for the vegetation patchwork
        RFe_tot=np.zeros([len(RFeveg[0])], dtype=float)
        # vegetation interception (weigthed average patchwork)
        INTER=np.zeros([Ttotal], dtype=np.float)
        # Ponding storage
        POND=np.zeros([Ttotal], dtype=np.float)
        # Ponding storage change
        dPOND=np.zeros([Ttotal], dtype=np.float)
        # Surface runoff (hortonian and saturated overland flow)
        Ro=np.zeros([Ttotal], dtype=np.float)
        #Deep percolation
        Rp=np.zeros([Ttotal,nsl], dtype=np.float)
        #Soil moisture storage
        S=np.zeros([Ttotal,nsl], dtype=np.float)
        #Soil moisture storage in volumetric percent of saturation
        Spc=np.zeros([Ttotal,nsl], dtype=np.float)
        #Soil moisture storage changes
        dS=np.zeros([Ttotal,nsl], dtype=np.float)
        #Actual evaporation from unsaturated zone
        Eu=np.zeros([Ttotal,nsl], dtype=np.float)
        #Actual transpiration from unsaturated zone
        Tu=np.zeros([Ttotal,nsl], dtype=np.float)
        #Actual evapotranspiration from surface
        Es=np.zeros([Ttotal], dtype=np.float)
        # GW seepage into soil zone
        SEEPAGE=np.zeros([Ttotal], dtype=int)
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
        nflux1 = 17
        nflux2 = 6
        results1 = np.zeros([Ttotal,nflux1])
        results2 = np.zeros([Ttotal,nsl,nflux2])

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

        Zr_elev = []
        for z in range(len(Zr)):
            Zr_elev.append(Dtop - Zr[z]*1000)

        # PROCESSING THE WHOLE DATA SET
        for t in range(int(nstp)):    # t: current time step

            # Preprocessing of PET/PE/INTER and RFe
            # for different vegetation
            SOILarea = 100
            for v in range(len(PETveg)):
                if VEGarea[v]<>self.hnoflo:
                    RFe_tot[t] += RFeveg[v,t]*VEGarea[v]/100
                    PET_tot[t] += PETveg[v,t]*VEGarea[v]/100
                    SOILarea -= VEGarea[v]
            RFe_tot[t] += RF[t]*SOILarea/100
            INTER[t] = RF[t] - RFe_tot[t]
            PE_tot[t] = PEsoil[t]*SOILarea/100

            # handle drycell
            if HEADS[t]>hdry-1E3:
                HEADStmp = Dbot - 100000
            else:
                HEADStmp=HEADS[t]*1000
            dtwt = Dtop-HEADStmp

            Sini = np.zeros([nsl], dtype=float)

            # heads below soil bottom
            if DRN[t]==0.0:
                # bottom boundary
                Dlbot[nsl-1] = HEADStmp
                Dl[nsl-1] = Dltop[nsl-1] - HEADStmp
                for l in range(nsl):
                    Si[l] = Si[l]*Dl[l]
                # fluxes
                Estmp, PONDtmp, Rotmp, Rptmp, Eutmp, Tutmp, Stmp, Spctmp, Si, Rtmp, Egtmp, Tgtmp, DRNcorr = self.flux(RFe_tot[t], PETveg[:,t], PE_tot[t], S, Sini, Rp, Zr_elev, VEGarea, HEADStmp, Dltop, Dlbot, Dl, nsl, Sm, Sfc, Sr, Ks, PONDm, PONDratio, t, Si, Rpi, DRN[t], DRNi, PONDi, E0[t], dtwt, st, i, j, n, k_Tu_slp, k_Tu_inter)

            # heads above soil bottom
            else:
                # bottom boundary
                Dlbot[nsl-1] = ELEV
                Dl[nsl-1] = 0.0
                for l in range(nsl):
                    Si[l] = Si[l]*Dl[l]
                # fluxes
                Estmp, PONDtmp, Rotmp, Rptmp, Eutmp, Tutmp, Stmp, Spctmp, Si, Rtmp, Egtmp, Tgtmp, DRNcorr = self.flux(RFe_tot[t], PETveg[:,t], PE_tot[t], S, Sini, Rp, Zr_elev, VEGarea, HEADStmp, Dltop, Dlbot, Dl, nsl, Sm, Sfc, Sr, Ks, PONDm, PONDratio, t, Si, Rpi, DRN[t], DRNi, PONDi, E0[t], dtwt, st, i, j, n, k_Tu_slp, k_Tu_inter)
            # fill the output arrays
            POND[t]=PONDtmp
            Ro[t] = Rotmp
            Es[t] = Estmp
            Eg[t] = Egtmp
            Tg[t] = Tgtmp
            R[t]  = Rtmp
            for l in range(nsl):
                S[t,l]=Stmp[l]
                Spc[t,l]=Spctmp[l]
                Rp[t,l]=Rptmp[l]
                Eu[t,l]=Eutmp[l]
                Tu[t,l]=Tutmp[l]
            Rn[t] = R[t] - Eg[t] -Tg[t] - DRNcorr
            ETg[t] = - Eg[t] -Tg[t] - DRNcorr
            if POND[t]>PONDi:
                dPOND[t] = POND[t] - PONDi
            # compute the water mass balance (MB) in the unsaturated zone
            RpinMB=0.0
            RpoutMB=0.0
            EuMB=0.0
            TuMB=0.0
            for l in range(nsl):
                EuMB += Eu[t,l]
                TuMB += Tu[t,l]
                RpoutMB += Rp[t,l]
                if t==0:
                    dS[t,l] = S[t,l]-Si[l]
                    if l<(nsl-1):
                        RpinMB += Rpi[l]
                else:
                    dS[t,l] = (S[t,l]-S[t-1,l])
                    if l<(nsl-1):
                        RpinMB += Rp[t-1,l]
            MB[t]=RF[t]+PONDi+DRN[t]+RpinMB-INTER[t]-Ro[t]-POND[t]-Es[t]-EuMB-TuMB-RpoutMB-sum(dS[t,:])
            # index = {'iRF':0, 'iPET':1, 'iPE':2, 'iRFe':3, 'iPOND':4, 'iRo':5, 'iSEEPAGE':6, 'iEs':7, 'iMB':8, 'iINTER':9, 'iE0':10, 'iEg':11, 'iTg':12, 'iR':13, 'iRn':14, 'idPOND':15, 'ETg':16}
            results1[t,:] = [RF[t], PET_tot[t], PE_tot[t], RFe_tot[t], POND[t], Ro[t], DRN[t], Es[t], MB[t], INTER[t], E0[t], Eg[t], Tg[t], R[t], Rn[t], dPOND[t], ETg[t]]
            # index_S = {'iEu':0, 'iTu':1,'iS':2, 'iRp':3}
            for l in range(nsl):
                results2[t,l,:] = [Eu[t,l], Tu[t,l], Spc[t,l], Rp[t,l], dS[t,l], S[t,l]]

        return results1, results2

        del RF, PET_tot, PE_tot, RFe_tot, POND, Ro, DRN, Es, MB, INTER, E0, Eg, Tg, R, Rn, dPOND, ETg
        del Eu, Tu, Spc, Rp, dS, S
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
        del R, h1, h_tmp, h

#EOF#