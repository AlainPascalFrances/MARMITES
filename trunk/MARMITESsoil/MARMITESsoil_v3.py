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
__version__ = "0.3"
__date__ = "2012"

import numpy as np

class SOIL:

    """
    SOIL class: compute the water balance in soil (depth-wise, each reservoir
    correspond to a soil horizon (A, B, C and rock)
    ## TO BE UPDATED""
    _______________________________________________________________________________

    INPUTS
            PARAMETERS
                MAXIL       Interception looses
                Sm          Max. Soil moiture storage
                Sfc         Soil moiture storage at field capacity
                Sr          Residual Soil moiture storage
                Ssoil_ini          Inicial soil moisture
                Ks          Saturated hydraulic condutivity
                Ssurf_max       Max. ponding capacity
            STATE VARIABLES
                RF           Daily rainfall
                PT           Daily transpiration
                PE           Daily evaporation
    OUTPUTS
            RFe              Daily Excess rainfall
            ETsoil             Daily evapotranspiration
            Ssoil               Daily soil moisture
            Rp              Daily percolation
            Ssurf            Daily ponding
            Ro              Daily runoff

    Provide daily values for the LINRES module
    ______________________________________________________________________________

    """

#####################

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
                       'loam'            : {'dll':33.0,'y0':0.004,'b':0.028, 'ext_d':265.0},
                       'silty clay'      : {'dll':37.0,'y0':0.007,'b':0.046, 'ext_d':335.0},
                       'clay loam'       : {'dll':33.0,'y0':0.008,'b':0.027, 'ext_d':405.0},
                       'silt loam'       : {'dll':38.0,'y0':0.006,'b':0.019, 'ext_d':420.0},
                       'silt'            : {'dll':31.0,'y0':0.007,'b':0.021, 'ext_d':430.0},
                       'silty clay loam' : {'dll':40.0,'y0':0.007,'b':0.021, 'ext_d':450.0},
                       'clay'            : {'dll':45.0,'y0':0.006,'b':0.019, 'ext_d':620.0},
                       'custom'          : {'dll':100.0,'y0':0.00,'b':0.013, 'ext_d':330.0}
}

#####################

    def flux(self, RFe, PT, PE, E0, Zr_elev, VEGarea, HEADS, TopSoilLay, BotSoilLay, Tl, nsl, Sm, Sfc, Sr, Ks, Ssurf_max, Ssurf_ratio, Ssoil_ini, Rp_ini, Ssurf_ini, EXF, dgwt, st, i, j, n, kTu_min, kTu_n, dt, dti, NVEG, LAIveg):

        def surfwater(s_tmp, Sm, Ssurf_max, E0, i, j, n, dt):
            '''
            Ponding and surface runoff function
            '''
            if s_tmp > Sm:
                Ssurf_tmp = s_tmp-Sm
                if (Ssurf_tmp - Ssurf_max) > 1.0E-7:
                    Ro_tmp = (Ssurf_tmp - Ssurf_max)/dt
                    Ssurf_tmp = Ssurf_max
                else:
                    Ro_tmp = 0.0
                if (Ssurf_tmp - E0) > 1.0E-7:
                    Esurf_tmp = E0
                    Ssurf_tmp = Ssurf_tmp - E0*dt
                else:
                    Esurf_tmp = Ssurf_tmp/dt
                    Ssurf_tmp = 0.0
            else:
                Ssurf_tmp = 0.0
                Ro_tmp = 0.0
                Esurf_tmp = 0.0
            return Ssurf_tmp, Ro_tmp, Esurf_tmp

        ##################

        def perc(s_tmp, Sm, Sfc, Ks, s_lp1, Sm_sp1, i, j, n, dt):
            '''
            Percolation function
            '''
            if (s_tmp - Sfc) <=  1.0E-7:
                rp_tmp = 0.0
            elif (s_tmp - Sfc) > 1.0E-7 and (s_tmp - Sm) < 1.0E-7:
                Sg = (s_tmp-Sfc)/(Sm-Sfc) # gravitational water
                if (Ks*Sg*dt) - (s_tmp-Sfc) > 1.0E-7:
                    rp_tmp = (s_tmp-Sfc)/dt
                else:
                    rp_tmp = Ks*Sg
            elif (s_tmp - Sm) > 1.0E-7:
                if (Ks*dt) - (Sm-Sfc) > 1.0E-7:
                    rp_tmp = (Sm-Sfc)/dt
                else:
                    rp_tmp = Ks
            return rp_tmp

        ##################

        def evp(s_tmp,Sm,Sr, pet, i, j, n, dt):
            '''
            Actual evapotranspiration function
            '''
            if (s_tmp - Sr) <= 1.0E-7:
                evp_tmp = 0.0
            elif ((s_tmp - Sr) > 1.0E-7) and (s_tmp - Sm) < 1.0E-7:
                Se = (s_tmp - Sr)/(Sm - Sr)
                if (pet*Se*dt - (s_tmp - Sr)) > 1.0E-7:
                    evp_tmp = (s_tmp - Sr)/dt
                else:
                    evp_tmp = pet*Se
            elif (s_tmp - Sr) > 1.0E-7:
                if (pet*dt-(Sm - Sr)) > 1.0E-7:
                    evp_tmp = (Sm - Sr)/dt
                else:
                    evp_tmp = pet
            return evp_tmp

        ##################
        if EXF < 0.0:
            print 'WARNING!\nEXF < 0.0, value %.6f corrected to 0.0.' % EXF
            EXF = 0.0
        SAT = np.zeros([nsl], dtype = bool)

         # Ssoil_tmp from previous time step
        # if first time step, use Ssoil_ini
        Ssoil_tmp = np.zeros([nsl])
        Rexf_tmp = np.zeros([nsl])
        for l in range(nsl):
            Ssoil_tmp[l] = Ssoil_ini[l]*1.0

        # INFILTRATION
        # first soil layer
        Ssoil_tmp[0] += RFe + Ssurf_ini
        # other soil layer, percolation from previous time step
        for l in range(1,nsl):
            Ssoil_tmp[l] += Rp_ini[l-1]*dti

        # GW EXF
        Ssoil_tmp[nsl-1] += EXF*dt

        # SOIL EXF
        llst = range(nsl)
        llst.reverse()
        for l in llst[:-1]:
            if Ssoil_tmp[l] > Sm[l]*Tl[l]:
                Rexf_tmp[l] += Ssoil_tmp[l] - Sm[l]*Tl[l]
                Ssoil_tmp[l-1] += Rexf_tmp[l]
                Ssoil_tmp[l] = Sm[l]*Tl[l]

        # SAT = True indicates saturation overland flow
        HEADS_corr = HEADS * 1.0
        dgwt_corr = dgwt * 1.0
        if EXF > 0.0:
            for l in llst:
                if Ssoil_tmp[l] >= (Sm[l]*Tl[l]):
                    HEADS_corr += Tl[l]
                    dgwt_corr -= Tl[l]
                    SAT[l] = True
                else:
                    break

        # Ssurf and Ro
        surfwater_tmp = surfwater(Ssoil_tmp[0],Sm[0]*Tl[0],Ssurf_max, E0*Ssurf_ratio, i, j, n, dt)
        Ssurf_tmp = surfwater_tmp[0]
        Ro_tmp = surfwater_tmp[1]
        Esurf_tmp = surfwater_tmp[2]
        Ssoil_tmp[0] -= (Ssurf_tmp + dt*(Ro_tmp + Esurf_tmp))

        Rexf_tmp /= dt

        Rp_tmp = np.zeros([nsl])
        Tsoil_tmpZr = np.zeros([nsl,len(Zr_elev)])
        Tsoil_tmp = np.zeros([nsl])
        Esoil_tmp = np.zeros([nsl])
        Ssoil_pc_tmp = np.zeros([nsl])
        # Groundwater transpiration
        Tg_tmp_Zr = np.zeros([nsl, len(Zr_elev)])
        Tg_tmp = 0.0
        # soil layers

        for l in range(nsl):
            # Esoil
            if Ssurf_tmp == 0.0:
                if PE > 0.0:
                    Esoil_tmp[l] = evp(Ssoil_tmp[l],Sm[l]*Tl[l], Sr[l]*Tl[l], PE, i, j, n, dt)
                    Ssoil_tmp[l] -= Esoil_tmp[l]*dt
                    PE -= Esoil_tmp[l]
            # Tsoil
            for v in range(NVEG):
                if PT[v] > 1E-7:
                    if LAIveg[v] > 1E-7:
                        if VEGarea[v] > 1E-7:
                            if BotSoilLay[l] > Zr_elev[v]:
                                Tsoil_tmpZr[l,v] = evp(Ssoil_tmp[l],Sm[l]*Tl[l], Sr[l]*Tl[l], PT[v], i, j, n, dt)
                            elif TopSoilLay[l] > Zr_elev[v] :
                                PTc = PT[v]*(TopSoilLay[l]-Zr_elev[v])/Tl[l]
                                Tsoil_tmpZr[l,v] = evp(Ssoil_tmp[l],Sm[l]*Tl[l], Sr[l]*Tl[l], PTc, i, j, n, dt)
                            Tsoil_tmp[l] += Tsoil_tmpZr[l,v]*VEGarea[v]*0.01
                            PT[v] -= Tsoil_tmpZr[l,v]
            Ssoil_tmp[l] -= Tsoil_tmp[l]*dt
            # Rp
            if l < (nsl-1):
                if SAT[l+1] == False:
                    Rp_tmp[l] = perc(Ssoil_tmp[l],Sm[l]*Tl[l],Sfc[l]*Tl[l],Ks[l], Ssoil_tmp[l+1], Sm[l+1]*Tl[l+1], i, j, n, dt)
            elif EXF == 0.0:
                Rp_tmp[l] = perc(Ssoil_tmp[l],Sm[l]*Tl[l],Sfc[l]*Tl[l], Ks[l],0.0,1.0E6, i, j, n, dt)
            Ssoil_tmp[l] -= Rp_tmp[l]*dt

        for l in range(nsl):
            Ssoil_pc_tmp[l] = Ssoil_tmp[l]/Tl[l]

        # GW evaporation, equation 17 of Shah et al 2007, see ref in the __init__
        if PE > 0.0:
            dgwt_corr *= 0.1
            y0    = self.paramEg[st]['y0']
            b     = self.paramEg[st]['b']
            dll   = self.paramEg[st]['dll']
            ext_d = self.paramEg[st]['ext_d']
            if dgwt_corr <= dll:
                Eg_tmp = PE
            elif dgwt_corr < ext_d:
                Eg_tmp = PE*(y0 + np.exp(-b*(dgwt_corr-dll)))
            else:
                Eg_tmp = 0.0
            dgwt_corr *= 10.0
        else:
            Eg_tmp = 0.0

        # Groundwater transpiration
        kTu_max = 0.99999
        for v in range(NVEG):
            if HEADS_corr > Zr_elev[v]:
                for l in range(nsl):
                    if Ssoil_pc_tmp[l] <> self.hnoflo:
                        if Ssoil_pc_tmp[l] >= Sr[l]:
                            if Ssoil_pc_tmp[l] <= Sm[l]:
                                kTu = kTu_min[v]+(kTu_max-kTu_min[v])*np.power(1-np.power(np.abs(Ssoil_pc_tmp[l]-Sm[l])/(Sm[l]-Sr[l]),kTu_n[v]),1/kTu_n[v])
                            else:
                                kTu = kTu_max
                        else:
                            kTu = kTu_min[v]
                        Tg_tmp_Zr[l,v] = Tsoil_tmpZr[l,v]*(1.0/kTu-1.0)
                        if Tg_tmp_Zr[l,v] > PT[v]:
                            Tg_tmp_Zr[l,v] = PT[v]
                        PT[v] -= Tg_tmp_Zr[l,v]
                        Tg_tmp += (Tg_tmp_Zr[l,v]*VEGarea[v]*0.01)

        return Esurf_tmp, Ssurf_tmp, Ro_tmp, Rp_tmp, Esoil_tmp, Tsoil_tmp, Ssoil_tmp, Ssoil_pc_tmp, Eg_tmp, Tg_tmp, HEADS_corr, dgwt_corr, SAT, Rexf_tmp
        del Esurf_tmp, Ssurf_tmp, Ro_tmp, Rp_tmp, Esoil_tmp, Tsoil_tmp, Ssoil_tmp, Ssoil_pc_tmp, Ssoil_ini, Eg_tmp, Tg_tmp, dgwt_corr

#####################

    def run(self, _nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _Ssoil_ini, botm_l0, _Ks,
            gridSOIL, gridSOILthick, TopSoil, gridMETEO,
            index, index_S, gridSsurfhmax, gridSsurfw,
            RF_veg_zoneSP, E0_zonesSP, PT_veg_zonesSP, RFe_veg_zonesSP, PE_zonesSP, gridVEGarea,
            LAI_veg_zonesSP, Zr, kTu_min, kTu_n, NVEG,
            cMF, conv_fact, h5_MF, h5_MM, irr_yn,
            RF_irr_zoneSP = [], PT_irr_zonesSP = [], RFe_irr_zoneSP = [],
            crop_irr_SP = [], gridIRR = [],
            Zr_c = [], kTu_min_c = [], kTu_n_c = [], NCROP = []):

        h_MF = None
        h_MF_mem = 'slow'
        try:
            h_MF = h5_MF['heads4MM'][:,:,:]
            h_MF_mem = 'fast'
        except:
            print '\nRAM memory too small compared to the size of the heads array -> slow computing.'
        if cMF.uzf_yn == 1:
            exf_MF = None
            exf_MF_mem = 'slow'
            try:
                # TODO this below assume that there is no grid refinement
                exf_MF = h5_MF['exf4MM'][:,:,:]*conv_fact/(cMF.delr[0]*cMF.delc[0])
                exf_MF_mem = 'fast'
            except:
                print '\nRAM memory too small compared to the size of the exfiltration array -> slow computing.'
        JD_nper = []
        # initial values of SP
        Ssoil_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol,_nslmax])
        Rp_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol,_nslmax])
        Ssurf_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol])
        Ssurf_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol])
        MM_finf_MF = np.zeros([cMF.nrow,cMF.ncol], dtype=float)
        MM_wel_MF  = np.zeros([cMF.nrow,cMF.ncol], dtype=float)
        for n in range(cMF.nper):
            if n == 0:
                JD_nper.append(cMF.JD[0])
            else:
                JD_nper.append(JD_nper[n-1] + cMF.perlen[n])
            tstart_MM = 0
            for t in range(n):
                tstart_MM += cMF.perlen[t]
            tend_MM = tstart_MM + cMF.perlen[n]
            tstart_MF = 0
            for t in range(n):
                tstart_MF += cMF.nstp[t]
            tend_MF = tstart_MF + cMF.nstp[n]
            MM         = np.zeros([cMF.perlen[n],cMF.nrow,cMF.ncol,len(index)], dtype=float)
            MM_S       = np.zeros([cMF.perlen[n],cMF.nrow,cMF.ncol,_nslmax,len(index_S)], dtype=float)
            if h_MF == None:
                h_MF       = h5_MF['heads4MM'][tstart_MF:tend_MF,:,:]
            if cMF.uzf_yn == 1:
                if exf_MF == None:
                    exf_MF = h5_MF['exf4MM'][tstart_MF:tend_MF,:,:]*conv_fact/(cMF.delr[j]*cMF.delc[i])
            # loop into the grid
            for i in range(cMF.nrow):
                for j in range(cMF.ncol):
                    if cMF.iuzfbnd[i][j] != 0.0:
                        SOILzone_tmp = gridSOIL[i,j]-1
                        METEOzone_tmp = gridMETEO[i,j]-1
                        _layer = cMF.iuzfbnd[i][j] - 1
                        slprop = _slprop[SOILzone_tmp]
                        nsl = _nsl[SOILzone_tmp]
                        # thickness of soil layers
                        Tl = gridSOILthick[i,j]*slprop*1000.0
                        # elevation of top and bottom of soil layers
                        TopSoilLay = np.zeros([nsl], dtype=float)
                        BotSoilLay = np.zeros([nsl], dtype=float)
                        for l in range(nsl):
                            if l==0:
                                TopSoilLay[l] = TopSoil[i,j]
                                BotSoilLay[l] = TopSoil[i,j]-Tl[l]
                            else:
                                TopSoilLay[l] = BotSoilLay[l-1]
                                BotSoilLay[l] = TopSoilLay[l]-Tl[l]
                        Ssoil_ini_tmp = []
                        if n == 0:
                            for l in range(nsl):
                                Ssoil_ini_tmp.append(_Ssoil_ini[SOILzone_tmp][l])
                            Ssurf_ini_tmp  = 0.0
                            dti = 1.0
                        else:
                            Ssoil_ini_tmp    = Ssoil_ini_tmp_array[i,j,:]
                            Ssurf_ini_tmp    = Ssurf_ini_tmp_array[i,j]
                        if irr_yn == 1:
                            IRRfield = gridIRR[i,j]
                        if  irr_yn == 1 and IRRfield > 0:
                            NVEG_tmp = 1
                            IRRfield -= 1
                            VEGarea_tmp = None
                            LAIveg_tmp = np.ones((NVEG_tmp,tend_MF-tstart_MF), dtype = np.float)
                            CROP_tmp = crop_irr_SP[IRRfield,tstart_MF:tend_MF]
                            RF_tmp = RF_irr_zoneSP[METEOzone_tmp,IRRfield,tstart_MF:tend_MF]
                            PT_zonesSP_tmp = np.zeros((NVEG_tmp,tend_MF-tstart_MF), dtype = np.float)
                            RFe_zonesSP_tmp = np.zeros((NVEG_tmp,tend_MF-tstart_MF), dtype = np.float)
                            RFe_zonesSP_tmp[0,:] = RFe_irr_zoneSP[METEOzone_tmp,IRRfield,tstart_MF:tend_MF]
                            PT_zonesSP_tmp[0,:] = PT_irr_zonesSP[METEOzone_tmp,IRRfield,tstart_MF:tend_MF]
                            Zr_tmp = Zr_c
                            kTu_min_tmp = kTu_min_c
                            kTu_n_tmp = kTu_n_c
                        else:
                            NVEG_tmp = NVEG
                            CROP_tmp = None
                            RF_tmp = RF_veg_zoneSP[METEOzone_tmp][tstart_MF:tend_MF]
                            PT_zonesSP_tmp = np.zeros((NVEG_tmp,tend_MF-tstart_MF), dtype = np.float)
                            RFe_zonesSP_tmp = np.zeros((NVEG_tmp,tend_MF-tstart_MF), dtype = np.float)
                            LAIveg_tmp = np.zeros((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float)
                            VEGarea_tmp=np.zeros([NVEG_tmp], dtype = np.float)
                            Zr_tmp = np.ones((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float)
                            kTu_min_tmp = np.ones((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float)
                            kTu_n_tmp = np.ones((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float)
                            for v in range(NVEG_tmp):
                                PT_zonesSP_tmp[v,:] = PT_veg_zonesSP[METEOzone_tmp,v,tstart_MF:tend_MF]
                                RFe_zonesSP_tmp[v,:] = RFe_veg_zonesSP[METEOzone_tmp,v,tstart_MF:tend_MF]
                                LAIveg_tmp[v,:] = LAI_veg_zonesSP[v,tstart_MF:tend_MF]
                                VEGarea_tmp[v] = gridVEGarea[v,i,j]
                                Zr_tmp[v] = Zr[v]
                                kTu_min_tmp[v] = kTu_min[v]
                                kTu_n_tmp[v] = kTu_n[v]
                        PE_zonesSP_tmp = PE_zonesSP[METEOzone_tmp,SOILzone_tmp,tstart_MF:tend_MF]
                        E0_zonesSP_tmp = E0_zonesSP[METEOzone_tmp][tstart_MF:tend_MF]
                        if h_MF_mem == 'slow':
                            h_MF_tmp   = h_MF[:,i,j]
                        elif h_MF_mem == 'fast':
                            h_MF_tmp   = h_MF[tstart_MF:tend_MF,i,j]
                        if cMF.uzf_yn == 1:
                            if exf_MF_mem == 'slow':
                                exf_MF_tmp = exf_MF[:,i,j]
                            elif exf_MF_mem == 'fast':
                                exf_MF_tmp = exf_MF[tstart_MF:tend_MF,i,j]
                        else:
                            exf_MF_tmp = 0.0
                        st         = _st[SOILzone_tmp]
                        Sm         = _Sm[SOILzone_tmp]
                        Sfc        = _Sfc[SOILzone_tmp]
                        Sr         = _Sr[SOILzone_tmp]
                        Ks         = _Ks[SOILzone_tmp]
                        Ssurf_max     = 1000*1.12*gridSsurfhmax[i,j]*gridSsurfw[i,j]/cMF.delr[j]
                        Ssurf_ratio   = 1.12*gridSsurfw[i,j]/cMF.delr[j]
                        EXF        = -exf_MF_tmp
                        RF         = RF_tmp
                        E0         = E0_zonesSP_tmp
                        PT         = PT_zonesSP_tmp
                        RFe        = RFe_zonesSP_tmp
                        PE         = PE_zonesSP_tmp
                        nstp       = cMF.nstp[n]
                        perlen     = cMF.perlen[n]
                        # Output initialisation
                        Ttotal = len(RF)
                        # PT for the vegetation patchwork
                        PT_tot = np.zeros([len(PT[0])], dtype = float)
                        # PE for the remaining bare soil
                        PE_tot = np.zeros([Ttotal], dtype = float)
                        # RFe for the vegetation patchwork
                        RFe_tot = np.zeros([len(RFe[0])], dtype = float)
                        # vegetation interception (weigthed average patchwork)
                        INTER_tot = np.zeros([Ttotal], dtype = np.float)
                        # Ponding storage
                        Ssurf = np.zeros([Ttotal], dtype = np.float)
                        # Ponding storage change
                        dSsurf = np.zeros([Ttotal], dtype = np.float)
                        # Surface runoff (hortonian and saturated overland flow)
                        Ro = np.zeros([Ttotal], dtype = np.float)
                        # Percolation
                        Rp = np.zeros([Ttotal,nsl], dtype = np.float)
                        # Percolation into UZF
                        inf = np.zeros([Ttotal], dtype = np.float)
                        # Exfiltration
                        Rexf = np.zeros([Ttotal,nsl], dtype = np.float)
                        #Soil moisture storage
                        Ssoil = np.zeros([Ttotal,nsl], dtype = np.float)
                        #Soil moisture storage in volumetric percent of saturation
                        Ssoil_pc = np.zeros([Ttotal,nsl], dtype = np.float)
                        #Total soil moisture storage changes
                        Ssoil_pc_tot = np.zeros([Ttotal], dtype = np.float)
                        #Soil moisture storage changes
                        dSsoil = np.zeros([Ttotal,nsl], dtype = np.float)
                        #Total soil moisture storage changes
                        dSsoil_tot = np.zeros([Ttotal], dtype = np.float)
                        #Actual evaporation from soil zone
                        Esoil = np.zeros([Ttotal,nsl], dtype = np.float)
                        #Actual transpiration from soil zone
                        Tsoil = np.zeros([Ttotal,nsl], dtype = np.float)
                        #Total actual evapotranspiration from soil zone
                        ETsoil_tot = np.zeros([Ttotal], dtype = np.float)
                        #Saturation overland flow
                        SAT = np.zeros([Ttotal, nsl], dtype = np.int)
                        #Actual evapotranspiration from surface
                        Esurf = np.zeros([Ttotal], dtype = np.float)
                        # HEADS
                        HEADS_MM = np.zeros([Ttotal], dtype = np.float)
                        # dgwt
                        dgwt = np.zeros([Ttotal], dtype = np.float)
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
                        nflux1 = 22
                        nflux2 = 9
                        MM_tmp = np.zeros([Ttotal,nflux1])
                        MM_S_tmp = np.zeros([Ttotal,nsl,nflux2])

                        dt = float(perlen)/float(nstp)

                        # PROCESSING THE WHOLE DATA SET
                        for t in range(int(nstp)):    # t: current time step
                            # Preprocessing of PT/PE/INTER and RFe
                            # for different vegetation
                            Zr_elev = []
                            SOILarea = 100.0
                            if CROP_tmp <> None:
                                Zr_elev.append(TopSoilLay[0] - Zr_tmp[CROP_tmp[t]-1]*1000.0)
                                if CROP_tmp[t] > 0:
                                    VEGarea_tmp = [100.0]
                                else:
                                    VEGarea_tmp = [0.0]
                                kTu_min_tmp = [kTu_min_tmp[CROP_tmp[t]-1]]
                                kTu_n_tmp = [kTu_n_tmp[CROP_tmp[t]-1]]
                            else:
                                for v in range(NVEG_tmp):
                                    Zr_elev.append(TopSoilLay[0] - Zr_tmp[v]*1000.0)
                                kTu_min_tmp = kTu_min_tmp[:]
                                kTu_n_tmp = kTu_n_tmp[:]
                            for v in range(NVEG_tmp):
                                if LAIveg_tmp[v,t] > 1E-7:
                                    RFe_tot[t] += RFe[v,t]*VEGarea_tmp[v]*0.01
                                    PT_tot[t]  += PT[v,t]*VEGarea_tmp[v]*0.01
                                    SOILarea   -= VEGarea_tmp[v]
                            RFe_tot[t]   += RF[t]*SOILarea*0.01
                            INTER_tot[t]  = RF[t] - RFe_tot[t]
                            PE_tot[t]     = PE[t]*SOILarea*0.01
                            # handle drycell
                            if h_MF_tmp[t] > (cMF.hdry-1E3):
                                HEADS_drycell = botm_l0[i,j]*1000.0
                            else:
                                HEADS_drycell = h_MF_tmp[t]*1000.0
                            # dgwt and uzthick
                            uzthick[t] = BotSoilLay[nsl-1] - HEADS_drycell
                            if EXF[t] <= 0.0: # or BotSoilLay[nsl-1] > HEADS_drycell:
                                dgwt[t] = TopSoilLay[0] - HEADS_drycell
                            elif EXF[t] > 0.0:
                                dgwt[t] = sum(Tl)
                                HEADS_drycell = BotSoilLay[nsl-1]
                            # for the first SP, S_ini is expressed in % and has to be converted in mm
                            if n == 0 and t == 0:
                                Ssoil_ini_tmp = Ssoil_ini_tmp * Tl
                            # fluxes
                            Esurf_tmp, Ssurf_tmp, Ro_tmp, Rp_tmp, Esoil_tmp, Tsoil_tmp, Ssoil_tmp, Ssoil_pc_tmp, Eg_tmp, Tg_tmp, HEADS_tmp, dgwt_tmp, SAT_tmp, Rexf_tmp = self.flux(RFe_tot[t], PT[:,t], PE_tot[t], E0[t], Zr_elev, VEGarea_tmp, HEADS_drycell, TopSoilLay, BotSoilLay, Tl, nsl, Sm, Sfc, Sr, Ks, Ssurf_max, Ssurf_ratio, Ssoil_ini_tmp, Rp_ini_tmp_array[i,j,:], Ssurf_ini_tmp, EXF[t], dgwt[t], st, i, j, n, kTu_min_tmp, kTu_n_tmp, dt, dti, NVEG_tmp, LAIveg_tmp[:,t])
                            # fill the output arrays
                            Ssurf[t]   = Ssurf_tmp
                            Ro[t]   = Ro_tmp
                            Esurf[t]   = Esurf_tmp
                            Eg[t]   = Eg_tmp
                            Tg[t]   = Tg_tmp
                            HEADS_MM[t] = HEADS_tmp
                            dgwt[t] = dgwt_tmp
                            Ssoil[t,:]    = Ssoil_tmp[:]
                            Ssoil_pc[t,:] = Ssoil_pc_tmp[:]
                            Ssoil_pc_tot[t] = sum(Ssoil_pc_tmp[:])/nsl
                            Rp[t,:]    = Rp_tmp[:]
                            inf[t]     = Rp_tmp[-1]
                            Rexf[t,:]  = Rexf_tmp[:]
                            Esoil[t,:]    = Esoil_tmp[:]
                            Tsoil[t,:]    = Tsoil_tmp[:]
                            SAT[t,:]   = SAT_tmp[:]
                            ETg[t] = Eg[t] + Tg[t]
                            dSsurf[t] = (Ssurf_ini_tmp - Ssurf[t])/dt
                            # compute the water mass balance (MB) in the soil zone
                            Rp_in_MB  = 0.0
                            Rp_out_MB = 0.0
                            Esoil_MB    = 0.0
                            Tsoil_MB    = 0.0
                            for l in range(nsl):
                                Esoil_MB     += Esoil[t,l]
                                Tsoil_MB     += Tsoil[t,l]
                                if l < (nsl-1):
                                    Rp_in_MB  += Rp_ini_tmp_array[i,j,l]*dti
                                Rp_out_MB += Rp[t,l]*dt
                                dSsoil[t,l] = (Ssoil_ini_tmp[l] - Ssoil[t,l])/dt
                                dSsoil_tot[t] += dSsoil[t,l]
                            dRp_tot = (Rp_in_MB - Rp_out_MB)/dt
                            ETsoil_tot[t] = Esoil_MB + Tsoil_MB
                            # MASS BALANCE COMPUTING
                            if nsl > 1:
                                # surficial soil layer
                                l = 0
                                MB_l[t,l] = (RFe_tot[t] + dSsoil[t,l] + dSsurf[t] + Rexf[t,l+1]) - (Rp[t,l] + Esoil[t,l] + Tsoil[t,l] + Ro_tmp + Esurf_tmp)
                                # intermediate soil layers
                                llst = range(1,nsl-1)
                                for l in llst:
                                    MB_l[t,l] = (Rp_ini_tmp_array[i,j,l-1]*dti - Rp[t,l]*dt)/dt + dSsoil[t,l] + Rexf[t,l+1] - (Esoil[t,l] + Tsoil[t,l] + Rexf[t,l])
                                # last soil layer
                                l = nsl-1
                                MB_l[t,l] = (Rp_ini_tmp_array[i,j,l-1]*dti - Rp[t,l]*dt)/dt + dSsoil[t,l] + EXF[t] - (Esoil[t,l] + Tsoil[t,l] + Rexf[t,l])
                            else:
                                # only one soil layer
                                l = 0
                                MB_l[t,l] = (RFe_tot[t] + dSsoil[t,l] + dSsurf[t] + EXF[t]) - (Rp[t,l] + Esoil[t,l] + Tsoil[t,l] + Ro_tmp + Esurf_tmp)
                                # last soil layer
                            # total mass balance for the soil
                            MB[t] = RFe_tot[t] + dSsurf[t] + dSsoil_tot[t] + dRp_tot + EXF[t] - (Ro[t] + Esurf[t] + Esoil_MB + Tsoil_MB)
                            # export list
                            # indexes of the HDF5 output arrays
                            # index = {'iRF':0, 'iPT':1, 'iPE':2, 'iRFe':3, 'iSsurf':4, 'iRo':5, 'iEXF':6, 'iEsurf':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSsurf':13, 'iETg':14, 'iETsoil':15, 'iSsoil_pc':16, 'idSsoil':17, 'iinf':18, 'iHEADScorr':19, 'idgwt':20, 'iuzthick':21}
                            MM_tmp[t,:] = [RF[t], PT_tot[t], PE_tot[t], RFe_tot[t], Ssurf[t], Ro[t], EXF[t], Esurf[t], MB[t], INTER_tot[t], E0[t], Eg[t], Tg[t], dSsurf[t], ETg[t], ETsoil_tot[t], Ssoil_pc_tot[t], dSsoil_tot[t], inf[t], HEADS_MM[t]*0.001, -dgwt[t]*0.001, uzthick[t]*0.001]
                            # index_S = {'iEsoil':0, 'iTsoil':1,'iSsoil_pc':2, 'iRp':3, 'iRexf':4, 'idSsoil':5, 'iSsoil':6, 'iSAT':7, 'iMB_l':8}
                            for l in range(nsl):
                                MM_S_tmp[t,l,:] = [Esoil[t,l], Tsoil[t,l], Ssoil_pc[t,l], Rp[t,l], Rexf[t,l], dSsoil[t,l], Ssoil[t,l], SAT[t,l], MB_l[t,l]]
                        if (float(cMF.perlen[n])/float(cMF.nstp[n])) != 1.0:
                            for stp in range(cMF.nstp[n]):
                                ts = float(cMF.perlen[n])/float(cMF.nstp[n])
                                tstart =    stp*ts
                                tend   =   (1+stp)*ts
                                for k in range(len(index)):
                                    MM[tstart:tend,i,j,k] = MM_tmp[stp,k]
                                for k in range(len(index_S)):
                                    for l in range(nsl):
                                        MM_S[tstart:tend,i,j,l,k] = MM_S_tmp[stp,l,k]
                        else:
                            for k in range(len(index)):
                                MM[:,i,j,k] = MM_tmp[:,k]
                            for k in range(len(index_S)):
                                for l in range(nsl):
                                    MM_S[:,i,j,l,k] = MM_S_tmp[:,l,k]
                        MM_finf_MF[i,j] = MM_S_tmp[:,nsl-1,index_S.get('iRp')].sum()/conv_fact
                        MM_wel_MF[i,j] = MM_tmp[:,index.get('iETg')].sum()/conv_fact
                        del MM_tmp, MM_S_tmp
                        # setting initial conditions for the next SP
                        Ssoil_ini_tmp_array[i,j,:]  = MM_S[cMF.nstp[n]-1,i,j,:,index_S.get('iSsoil')]
                        Rp_ini_tmp_array[i,j,:]  = MM_S[cMF.nstp[n]-1,i,j,:,index_S.get('iRp')]
                        Ssurf_ini_tmp_array[i,j]    = MM[cMF.nstp[n]-1,i,j,index.get('iSsurf')]
                    else:
                        if cMF.perlen[n]>1:
                            MM[:,i,j,:] = cMF.hnoflo
                            MM_S[:,i,j,:,:] = cMF.hnoflo
                        else:
                            MM[:,i,j,:] = cMF.hnoflo
                            MM_S[:,i,j,:,:] = cMF.hnoflo
                        MM_finf_MF[i,j] = 0.0
                        MM_wel_MF[i,j] = 0.0
            dti = float(cMF.perlen[n])/float(cMF.nstp[n])
            h5_MM['MM'][tstart_MM:tend_MM,:,:,:]     = MM[:,:,:,:]
            h5_MM['MM_S'][tstart_MM:tend_MM,:,:,:,:] = MM_S[:,:,:,:,:]
            h5_MM['finf'][n,:,:]                     = MM_finf_MF
            h5_MM['ETg'][n,:,:]                      = MM_wel_MF
        del MM, MM_S, MM_finf_MF, MM_wel_MF, Ssoil_ini_tmp_array, Rp_ini_tmp_array, Ssurf_ini_tmp_array, dti
        del h_MF, exf_MF, h_MF_tmp, exf_MF_tmp
        h5_MM.close()

#####################

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
            h_tmp = (h1[t] + h0*1000.0)*0.001
            h[t] = h_tmp

        return h
        del R, h1, h_tmp, h

#EOF#