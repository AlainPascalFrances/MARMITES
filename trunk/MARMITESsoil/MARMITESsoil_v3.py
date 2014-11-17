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

class clsMMsoil:

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
                       'custom'          : {'dll':100.0,'y0':0.000,'b':0.013, 'ext_d':330.0}
                       }
                     #                        'custom'          : {'dll':100.0,'y0':0.00,'b':0.013, 'ext_d':330.0}
#####################

    def flux(self, cMF, perleni, RFe, PT, PE, E0surf_max, Zr_elev, VEGarea, HEADS, TopSoilLay, BotSoilLay, Tl, nsl, Sm, Sfc, Sr, Ks, Ssurf_max, Ssoil_ini, Rp_ini, Ssurf_ini, EXF, dgwt, st, i, j, n, kT_min, kT_max, kT_n, NVEG, LAIveg):

        def surfwater(s_tmp, Sm, Ssurf_max, E0, perlen):
            '''
            Ponding and surface runoff function
            '''
            if (Sm - s_tmp) < 1.0E-7:
                Ssurf_tmp = s_tmp - Sm
                if (Ssurf_tmp - Ssurf_max) > 1.0E-7:
                    Ro_tmp = (Ssurf_tmp - Ssurf_max)/perlen
                    Ssurf_tmp = Ssurf_max
                else:
                    Ro_tmp = 0.0
                if (Ssurf_tmp - E0) > 1.0E-7:
                    Esurf_tmp = E0
                    Ssurf_tmp = Ssurf_tmp - E0*perlen
                else:
                    Esurf_tmp = Ssurf_tmp/perlen
                    Ssurf_tmp = 0.0
            else:
                Ssurf_tmp = 0.0
                Ro_tmp = 0.0
                Esurf_tmp = 0.0
            return Ssurf_tmp, Ro_tmp, Esurf_tmp

        ##################

        def perc(s_tmp, Sm, Sfc, Ks, perlen):
            '''
            Percolation function
            '''
            if (s_tmp - Sfc) < 1.0E-7:
                rp_tmp = 0.0
            elif (s_tmp - Sfc) > 1.0E-7 and (s_tmp - Sm) < 1.0E-7:
                Sg = (s_tmp-Sfc)/(Sm-Sfc) # gravitational water
                if (Ks*Sg*perlen) - (s_tmp-Sfc) > 1.0E-7:
                    rp_tmp = (s_tmp-Sfc)/perlen
                else:
                    rp_tmp = Ks*Sg
            elif (s_tmp - Sm) > 1.0E-7:
                if (Ks) - (Sm-Sfc) > 1.0E-7:
                    rp_tmp = (Sm-Sfc)//perlen
                else:
                    rp_tmp = Ks
            return rp_tmp

        ##################

        def evp(s_tmp,Sm,Sr, pet, perlen):
            '''
            Actual evapotranspiration function
            '''
            if (s_tmp - Sr) < 1.0E-7:
                evp_tmp = 0.0
            elif ((s_tmp - Sr) > 1.0E-7) and (s_tmp - Sm) < 1.0E-7:
                Se = (s_tmp - Sr)/(Sm - Sr)
                if (pet*Se*perlen - (s_tmp - Sr)) > 1.0E-7:
                    evp_tmp = (s_tmp - Sr)/perlen
                else:
                    evp_tmp = pet*Se
            elif (s_tmp - Sr) > 1.0E-7:
                if (pet*perlen-(Sm - Sr)) > 1.0E-7:
                    evp_tmp = (Sm - Sr)/perlen
                else:
                    evp_tmp = pet
            return evp_tmp

        ##################
        if EXF < 0.0:
            print 'WARNING!\nEXFg < 0.0, value %.6f corrected to 0.0.' % EXF
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
            Ssoil_tmp[l] += Rp_ini[l-1]*perleni

        # GW EXF
        perlen = cMF.perlen[n]
        Ssoil_tmp[nsl-1] += EXF*perlen

        # SOIL EXF
        llst = range(nsl)
        llst.reverse()
        for l in llst[:-1]:
            if Ssoil_tmp[l] >= (Sm[l]*Tl[l]):
                Rexf_tmp[l] += Ssoil_tmp[l] - Sm[l]*Tl[l]
                Ssoil_tmp[l-1] += Rexf_tmp[l]
                Ssoil_tmp[l] = Sm[l]*Tl[l]
        Rexf_tmp /= perlen

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
        surfwater_tmp = surfwater(Ssoil_tmp[0], Sm[0]*Tl[0], Ssurf_max, E0surf_max, perlen)
        Ssurf_tmp = surfwater_tmp[0]
        Ro_tmp = surfwater_tmp[1]
        Esurf_tmp = surfwater_tmp[2]
        Ssoil_tmp[0] -= (Ssurf_tmp + perlen*(Ro_tmp + Esurf_tmp))

        Rp_tmp = np.zeros([nsl])
        Tsoil_tmpZr = np.zeros([nsl,len(Zr_elev)])
        Tsoil_tmp = np.zeros([nsl])
        Esoil_tmp = np.zeros([nsl])
        Ssoil_pc_tmp = np.zeros([nsl])

        # soil layers
        for l in range(nsl):
            # Rp
            if l < (nsl-1):
                if SAT[l+1] == False:
                    Rp_tmp[l] = perc(Ssoil_tmp[l], Sm[l]*Tl[l], Sfc[l]*Tl[l], Ks[l], perlen)
            elif EXF == 0.0:
                Rp_tmp[l] = perc(Ssoil_tmp[l], Sm[l]*Tl[l], Sfc[l]*Tl[l], Ks[l], perlen)
            Ssoil_tmp[l] -= Rp_tmp[l]*perlen
            # Esoil
            if Ssurf_tmp == 0.0:
                if PE > 0.0:
                    Esoil_tmp[l] = evp(Ssoil_tmp[l],Sm[l]*Tl[l], Sr[l]*Tl[l], PE, perlen)
                    Ssoil_tmp[l] -= Esoil_tmp[l]*perlen
                    PE -= Esoil_tmp[l]
            # Tsoil
            for v in range(NVEG):
                if PT[v] > 1.0E-7:
                    if LAIveg[v] > 1.0E-7:
                        if VEGarea[v] > 1.0E-7:
                            if BotSoilLay[l] > Zr_elev[v]:
                                Tsoil_tmpZr[l,v] = evp(Ssoil_tmp[l],Sm[l]*Tl[l], Sr[l]*Tl[l], PT[v], perlen)
                            elif TopSoilLay[l] > Zr_elev[v] :
                                PTc = PT[v]*(TopSoilLay[l]-Zr_elev[v])/Tl[l]
                                Tsoil_tmpZr[l,v] = evp(Ssoil_tmp[l],Sm[l]*Tl[l], Sr[l]*Tl[l], PTc, perlen)
                            Tsoil_tmp[l] += Tsoil_tmpZr[l,v]*VEGarea[v]*0.01
                            PT[v] -= Tsoil_tmpZr[l,v]
                            if PT[v] < 0.0:
                               PT[v] = 0.0 
            Ssoil_tmp[l] -= Tsoil_tmp[l]*perlen

        for l in range(nsl):
            Ssoil_pc_tmp[l] = Ssoil_tmp[l]/Tl[l]
            
        sy_tmp = cMF.cPROCESS.float2array(cMF.sy_actual)[cMF.outcropL[i,j]-1,i,j]
        # GW evaporation Eg, equation 17 of Shah et al 2007, see ref in the __init__
        if cMF.wel_yn == 1:
            if Ssurf_tmp == 0.0:
                if PE > 0.0:
                    dgwt_corr *= 0.1
                    y0    = self.paramEg[st]['y0']
                    b     = self.paramEg[st]['b']
                    dll   = self.paramEg[st]['dll']
                    ext_d = self.paramEg[st]['ext_d']
                    if dgwt_corr <= dll:
                        Eg_tmp = PE*1.0
                    elif dgwt_corr < ext_d:
                        Eg_tmp = PE*(y0 + np.exp(-b*(dgwt_corr-dll)))
                        dgwt_corr += Eg_tmp/sy_tmp
                        HEADS_corr -= Eg_tmp/sy_tmp
                    else:
                        Eg_tmp = 0.0
                    dgwt_corr *= 10.0
                else:
                    Eg_tmp = 0.0
            else:
                Eg_tmp = 0.0
        else:
            Eg_tmp = 0.0                

        # Groundwater transpiration Tg
        Tg_tmp_Zr = np.zeros([nsl, len(Zr_elev)])
        Tg_tmp = 0.0
        if cMF.wel_yn == 1:
            order = np.asarray(Zr_elev).flatten().argsort()
#            print 'Zr_elev', Zr_elev, np.array(Zr_elev)[order]
#            print 'NVEG', range(NVEG), np.array(range(NVEG))[order]
#            print 'kT_max', kT_max, np.array(kT_max)[order]
            for jj, (Zr_elev_, v, kT_min_, kT_max_, kT_n_) in enumerate(zip(np.array(Zr_elev)[order],np.array(range(NVEG))[order], np.array(kT_min)[order], np.array(kT_max)[order], np.array(kT_n)[order])):
                #print '----'
                if HEADS_corr > Zr_elev_:
                    for l in range(nsl):
                        #print 'veg%d, soil layer %d' %(v,l)
                        if Ssoil_pc_tmp[l] != cMF.hnoflo:
                            if Ssoil_pc_tmp[l] >= Sr[l]:
                                if Ssoil_pc_tmp[l] <= Sm[l]:
                                    kT = kT_min_+(kT_max_-kT_min_)*np.power(1-np.power(np.abs(Ssoil_pc_tmp[l]-Sm[l])/(Sm[l]-Sr[l]),kT_n_),1/kT_n_)
                                else:
                                    kT = kT_max_
                            else:
                                kT = kT_min_
                            Tg_tmp_Zr[l,v] = PT[v]*kT
                            PT[v] -= Tg_tmp_Zr[l,v]
                            Tg_tmp1 = (Tg_tmp_Zr[l,v]*VEGarea[v]*0.01)
                            if Tg_tmp1 > 1E-6 and (HEADS_corr-Tg_tmp1/sy_tmp) < Zr_elev_:
                                #print 'WARNING at cell [i=%d,j=%d,L=%d] with NVEG = %d\nTg = %.2f, HEADS_corr = %.2f, dgwt_corr= %.2f, Zr_elev = %.2f, Sy = %.3f' % (i+1,j+1,cMF.outcropL[i,j], v, Tg_tmp1, HEADS_corr, dgwt_corr, Zr_elev_, sy_tmp)
                                Tg_tmp1 = (HEADS_corr - Zr_elev[v])*sy_tmp
                                #print 'New Tg = %.2f mm' % Tg_tmp1
#                            else:                        
#                                print 'Tg veg%d = %.2f mm' %(v,Tg_tmp1)
                            Tg_tmp += Tg_tmp1
                            dgwt_corr += Tg_tmp1/sy_tmp
                            HEADS_corr -= Tg_tmp1/sy_tmp                                           
#                else:
#                    print 'veg%d: root too short!' % v
        else:
            Tg_tmp = 0.0
        
        del sy_tmp
    
        return Esurf_tmp, Ssurf_tmp, Ro_tmp, Rp_tmp, Esoil_tmp, Tsoil_tmp, Ssoil_tmp, Ssoil_pc_tmp, Eg_tmp, Tg_tmp, HEADS_corr, dgwt_corr, SAT, Rexf_tmp

#####################

    def runMMsoil(self, _nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _Ssoil_ini, botm_l0, _Ks,
            gridSOIL, gridSOILthick, TopSoil, gridMETEO,
            index, index_S, gridSsurfhmax, gridSsurfw,
            RF_veg_zoneSP, E0_zonesSP, PT_veg_zonesSP, RFe_veg_zonesSP, PE_zonesSP, gridVEGarea,
            LAI_veg_zonesSP, Zr, kT_min, kT_max, kT_n, NVEG,
            cMF, conv_fact, h5_MF, h5_MM, irr_yn,
            RF_irr_zoneSP = [], PT_irr_zonesSP = [], RFe_irr_zoneSP = [],
            crop_irr_SP = [], gridIRR = [],
            Zr_c = [], kT_min_c = [], kT_max_c = [], kT_n_c = [], NCROP = []):

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
        tstart_MM = 0
        tstart_MF = 0
        # initial values of SP
        Ssoil_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol,_nslmax])
        Rp_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol,_nslmax])
        Ssurf_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol])
        Ssurf_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol])
        MM_finf_MF = np.zeros([cMF.nrow,cMF.ncol], dtype=float)
        MM_wel_MF  = np.zeros([cMF.nrow,cMF.ncol], dtype=float)
        for n in range(cMF.nper):
            if n > 0:
                tstart_MM += cMF.perlen[n-1]
                tstart_MF += cMF.nstp[n-1]
            tend_MM = tstart_MM + cMF.perlen[n]
            tend_MF = tstart_MF + cMF.nstp[n]
            MM         = np.zeros([cMF.perlen[n],cMF.nrow,cMF.ncol,len(index)], dtype=float)
            MM_S       = np.zeros([cMF.perlen[n],cMF.nrow,cMF.ncol,_nslmax,len(index_S)], dtype=float)
            if h_MF == None:
                h_MF       = h5_MF['heads4MM'][tstart_MF:tend_MF,:,:]
            if cMF.uzf_yn == 1:
                if exf_MF == None:
                    exf_MF = h5_MF['exf4MM'][tstart_MF:tend_MF,:,:]*conv_fact/(cMF.delr[0]*cMF.delc[0])
            # loop into the grid
            for i in range(cMF.nrow):
                for j in range(cMF.ncol):
                    if cMF.outcropL[i,j] > 0:
                        SOILzone_tmp = gridSOIL[i,j]-1
                        METEOzone_tmp = gridMETEO[i,j]-1
                        slprop = _slprop[SOILzone_tmp]
                        nsl = _nsl[SOILzone_tmp]
                        # thickness of soil layers
                        Tl = gridSOILthick[i,j]*slprop*1000.0
                        # elevation of top and bottom of soil layers
                        TopSoilLay = np.zeros([nsl], dtype=float)
                        BotSoilLay = np.zeros([nsl], dtype=float)
                        for l in range(nsl):
                            if l==0:
                                TopSoilLay[0] = TopSoil[i,j]
                                BotSoilLay[0] = TopSoil[i,j]-Tl[0]
                            else:
                                TopSoilLay[l] = BotSoilLay[l-1]
                                BotSoilLay[l] = TopSoilLay[l]-Tl[l]
                        Ssoil_ini_tmp = []
                        if n == 0:
                            for l in range(nsl):
                                Ssoil_ini_tmp.append(_Ssoil_ini[SOILzone_tmp][l])
                            Ssurf_ini_tmp  = 0.0
                            perleni = 1.0
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
                            kT_min_tmp = kT_min_c
                            kT_max_tmp = kT_max_c
                            kT_n_tmp = kT_n_c
                        else:
                            NVEG_tmp = NVEG
                            CROP_tmp = None
                            RF_tmp = RF_veg_zoneSP[METEOzone_tmp][tstart_MF:tend_MF]
                            PT_zonesSP_tmp = np.zeros((NVEG_tmp,tend_MF-tstart_MF), dtype = np.float)
                            RFe_zonesSP_tmp = np.zeros((NVEG_tmp,tend_MF-tstart_MF), dtype = np.float)
                            LAIveg_tmp = np.zeros((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float)
                            VEGarea_tmp = np.zeros([NVEG_tmp], dtype = np.float)
                            Zr_tmp = np.ones((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float)
                            kT_min_tmp = np.ones((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float)
                            kT_max_tmp = np.ones((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float)
                            kT_n_tmp = np.ones((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float)
                            for v in range(NVEG_tmp):
                                PT_zonesSP_tmp[v,:]  = PT_veg_zonesSP[METEOzone_tmp,v,tstart_MF:tend_MF]
                                RFe_zonesSP_tmp[v,:] = RFe_veg_zonesSP[METEOzone_tmp,v,tstart_MF:tend_MF]
                                LAIveg_tmp[v,:]      = LAI_veg_zonesSP[v,tstart_MF:tend_MF]
                                VEGarea_tmp[v]       = gridVEGarea[v,i,j]
                                Zr_tmp[v]            = Zr[v]
                                kT_min_tmp[v]        = kT_min[v]
                                kT_max_tmp[v]        = kT_max[v]
                                kT_n_tmp[v]          = kT_n[v]
                        PE_zonesSP_tmp = PE_zonesSP[METEOzone_tmp,SOILzone_tmp,tstart_MF:tend_MF]
                        E0_zonesSP_tmp = E0_zonesSP[METEOzone_tmp][tstart_MF:tend_MF]
                        if h_MF_mem == 'slow':
                            h_MF_tmp   = h_MF[:,i,j]
                        elif h_MF_mem == 'fast':
                            h_MF_tmp   = h_MF[tstart_MF:tend_MF,i,j]
                        if cMF.uzf_yn == 1:
                            if exf_MF_mem == 'slow':
                                exf_MF_tmp = -exf_MF[:,i,j]
                            elif exf_MF_mem == 'fast':
                                exf_MF_tmp = -exf_MF[tstart_MF:tend_MF,i,j]
                        else:
                            exf_MF_tmp = 0.0
                        st         = _st[SOILzone_tmp]
                        Sm         = _Sm[SOILzone_tmp]
                        Sfc        = _Sfc[SOILzone_tmp]
                        Sr         = _Sr[SOILzone_tmp]
                        Ks         = _Ks[SOILzone_tmp]
                        Ssurf_max  = np.power(cMF.delr[j],3)*gridSsurfhmax[i,j]*gridSsurfw[i,j]*1.126847784/np.power(100.0,2)/10.0   #1000*1.12*gridSsurfhmax[i,j]*gridSsurfw[i,j]/cMF.delr[j]
                        E0surf_max = cMF.delr[j]*gridSsurfw[i,j]*1.126847784*E0_zonesSP_tmp/np.power(100.0,2)  #1.12*gridSsurfw[i,j]/cMF.delr[j]
                        # Output initialisation
                        # PT for the vegetation patchwork
                        PT_tot = np.zeros([len(PT_zonesSP_tmp[0])], dtype = float)
                        # RFe for the vegetation patchwork
                        RFe_tot = np.zeros([len(RFe_zonesSP_tmp[0])], dtype = float)
                        #Soil moisture storage changes
                        dSsoil = np.zeros([nsl], dtype = np.float)
                        # MASS BALANCE each soil layer
                        MB_l = np.zeros([nsl], dtype = np.float)
                        # output arrays
                        nflux1 = 22
                        nflux2 = 9
                        MM_tmp = np.zeros([nflux1])
                        MM_S_tmp = np.zeros([nsl,nflux2])

                        # PROCESSING THE WHOLE DATA SET
                        # Preprocessing of PT/PE/INTER and RFe
                        # for different vegetation
                        Zr_elev = []
                        SOILarea = 100.0
                        if CROP_tmp != None:
                            Zr_elev.append(TopSoilLay[0] - Zr_tmp[CROP_tmp-1]*1000.0)
                            if CROP_tmp > 0:
                                VEGarea_tmp = [100.0]
                            else:
                                VEGarea_tmp = [0.0]
                            kT_min_tmp = [kT_min_tmp[CROP_tmp-1]]
                            kT_max_tmp = [kT_max_tmp[CROP_tmp-1]]
                            kT_n_tmp = [kT_n_tmp[CROP_tmp-1]]
                        else:
                            for v in range(NVEG_tmp):
                                Zr_elev.append(TopSoilLay[0] - Zr_tmp[v]*1000.0)
                            kT_min_tmp = kT_min_tmp[:]
                            kT_max_tmp = kT_max_tmp[:]
                            kT_n_tmp = kT_n_tmp[:]
                        for v in range(NVEG_tmp):
                            if LAIveg_tmp[v] > 1.0E-7:
                                RFe_tot += RFe_zonesSP_tmp[v]*VEGarea_tmp[v]*0.01
                                PT_tot  += PT_zonesSP_tmp[v]*VEGarea_tmp[v]*0.01
                                SOILarea -= VEGarea_tmp[v]
                        RFe_tot   += RF_tmp*SOILarea*0.01
                        INTER_tot  = RF_tmp - RFe_tot
                        PE_tot     = PE_zonesSP_tmp*SOILarea*0.01
                        # handle drycell
                        if np.abs(h_MF_tmp - cMF.hdry) < 1.0E-5:
                            HEADS_drycell = botm_l0[i,j]*1000.0
                        else:
                            HEADS_drycell = h_MF_tmp*1000.0
                        # dgwt and uzthick
                        if exf_MF_tmp <= 0.0: # or BotSoilLay[nsl-1] > HEADS_drycell:
                            dgwt = TopSoilLay[0] - HEADS_drycell
                        elif exf_MF_tmp > 0.0:
                            dgwt = sum(Tl)
                            HEADS_drycell = BotSoilLay[nsl-1]
                        uzthick = BotSoilLay[nsl-1] - HEADS_drycell
                        # for the first SP, S_ini is expressed in % and has to be converted in mm
                        if n == 0:
                            Ssoil_ini_tmp = Ssoil_ini_tmp * Tl
                        # fluxes
                        Esurf_tmp, Ssurf_tmp, Ro_tmp, Rp_tmp, Esoil_tmp, Tsoil_tmp, Ssoil_tmp, Ssoil_pc_tmp, Eg_tmp, Tg_tmp, HEADS_MM, dgwt_tmp, SAT_tmp, Rexf_tmp = self.flux(cMF, perleni, RFe_tot, PT_zonesSP_tmp[:], PE_zonesSP_tmp*SOILarea*0.01, E0surf_max, Zr_elev, VEGarea_tmp, HEADS_drycell, TopSoilLay, BotSoilLay, Tl, nsl, Sm, Sfc, Sr, Ks, Ssurf_max, Ssoil_ini_tmp, Rp_ini_tmp_array[i,j,:], Ssurf_ini_tmp, exf_MF_tmp, dgwt, st, i, j, n, kT_min_tmp, kT_max_tmp, kT_n_tmp, NVEG_tmp, LAIveg_tmp[:])
                        Ssoil_pc_tot = sum(Ssoil_pc_tmp[:])/nsl
                        inf     = Rp_tmp[-1]
                        ETg = Eg_tmp + Tg_tmp
                        dSsurf = (Ssurf_tmp - Ssurf_ini_tmp)/cMF.perlen[n]
                        # compute the water mass balance (MB) in the soil zone
                        Rp_in_MB  = 0.0
                        Rp_out_MB = 0.0
                        Esoil_MB    = 0.0
                        Tsoil_MB    = 0.0
                        dSsoil_tot = 0.0
                        for l in range(nsl):
                            Esoil_MB     += Esoil_tmp[l]
                            Tsoil_MB     += Tsoil_tmp[l]
                            if l < (nsl-1):
                                Rp_in_MB  += Rp_ini_tmp_array[i,j,l]*perleni
                            Rp_out_MB += Rp_tmp[l]*cMF.perlen[n]
                            dSsoil[l] = (Ssoil_tmp[l]-Ssoil_ini_tmp[l])/cMF.perlen[n]
                            dSsoil_tot += dSsoil[l]
                        dRp_tot = (Rp_in_MB - Rp_out_MB)/cMF.perlen[n]
                        ETsoil_tot = Esoil_MB + Tsoil_MB
                        # MASS BALANCE COMPUTING
                        if nsl > 1:
                            # surficial soil layer
                            l = 0
                            MB_l[l] = (RFe_tot + Rexf_tmp[l+1]) - (Rp_tmp[l] + Esoil_tmp[l] + Tsoil_tmp[l] + Ro_tmp + Esurf_tmp + dSsoil[l] + dSsurf)
                            # intermediate soil layers
                            llst = range(1,nsl-1)
                            for l in llst:
                                MB_l[l] = (Rp_ini_tmp_array[i,j,l-1]*perleni - Rp_tmp[l]*cMF.perlen[n])/cMF.perlen[n] + Rexf_tmp[l+1] - (Esoil_tmp[l] + Tsoil_tmp[l] + Rexf_tmp[l]+ dSsoil[l] )
                            # last soil layer
                            l = nsl-1
                            MB_l[l] = (Rp_ini_tmp_array[i,j,l-1]*perleni - Rp_tmp[l]*cMF.perlen[n])/cMF.perlen[n] + exf_MF_tmp - (Esoil_tmp[l] + Tsoil_tmp[l] + Rexf_tmp[l] + dSsoil[l] )
                        else:
                            # only one soil layer
                            l = 0
                            MB_l[l] = (RFe_tot + exf_MF_tmp) - (Rp_tmp[l] + Esoil_tmp[l] + Tsoil_tmp[l] + Ro_tmp + Esurf_tmp + dSsurf + dSsoil[l])
                            # last soil layer
                        # total mass balance for the soil
                        MB = RFe_tot + dRp_tot + exf_MF_tmp - (Ro_tmp + Esurf_tmp + Esoil_MB + Tsoil_MB + dSsoil_tot + dSsurf)
                        # export list
                        # indexes of the HDF5 output arrays
                        # index_MM = {'iRF':0, 'iPT':1, 'iPE':2, 'iRFe':3, 'iSsurf':4, 'iRo':5, 'iEXFg':6, 'iEsurf':7, 'iMB':8, 'iI':9, 'iE0':10, 'iEg':11, 'iTg':12, 'idSsurf':13, 'iETg':14, 'iETsoil':15, 'iSsoil_pc':16, 'idSsoil':17, 'iinf':18, 'ihcorr':19, 'idgwt':20, 'iuzthick':21}
                        MM_tmp[:] = [RF_tmp, PT_tot, PE_tot, RFe_tot, Ssurf_tmp, Ro_tmp, exf_MF_tmp, Esurf_tmp, MB, INTER_tot, E0_zonesSP_tmp, Eg_tmp, Tg_tmp, dSsurf, ETg, ETsoil_tot, Ssoil_pc_tot, dSsoil_tot, inf, HEADS_MM*0.001, -dgwt_tmp*0.001, uzthick*0.001]
                        # index_MM_soil = {'iEsoil_l':0, 'iTsoil_l':1,'iSsoil_pc_l':2, 'iRp_l':3, 'iEXFg_l':4, 'idSsoil_l':5, 'iSsoil_l':6, 'iSAT_l':7, 'iMB_l':8}
                        for l in range(nsl):
                            MM_S_tmp[l,:] = [Esoil_tmp[l], Tsoil_tmp[l], Ssoil_pc_tmp[l], Rp_tmp[l], Rexf_tmp[l], dSsoil[l], Ssoil_tmp[l], SAT_tmp[l], MB_l[l]]

                        # fill MF and MM output arrays
                        for k in range(len(index)):
                            MM[:,i,j,k] = MM_tmp[k]
                        for k in range(len(index_S)):
                            for l in range(nsl):
                                MM_S[:,i,j,l,k] = MM_S_tmp[l,k]
                        # # Volumetric recharge rate UZF1 package
                        MM_finf_MF[i,j] = MM_S_tmp[nsl-1,index_S.get('iRsoil_l')]/conv_fact
                        # Volumetric recharge rate WEL package
                        MM_wel_MF[i,j] = MM_tmp[index.get('iETg')]/conv_fact
                        del MM_tmp, MM_S_tmp
                        # setting initial conditions for the next SP
                        Ssoil_ini_tmp_array[i,j,:]  = MM_S[cMF.perlen[n]-1,i,j,:,index_S.get('iSsoil_l')]
                        Rp_ini_tmp_array[i,j,:]     = MM_S[cMF.perlen[n]-1,i,j,:,index_S.get('iRsoil_l')]
                        Ssurf_ini_tmp_array[i,j]    = MM[cMF.perlen[n]-1,i,j,index.get('iSsurf')]
                    else:
                        if cMF.perlen[n]>1:
                            MM[:,i,j,:] = cMF.hnoflo
                            MM_S[:,i,j,:,:] = cMF.hnoflo
                        else:
                            MM[:,i,j,:] = cMF.hnoflo
                            MM_S[:,i,j,:,:] = cMF.hnoflo
                        MM_finf_MF[i,j] = 0.0
                        MM_wel_MF[i,j] = 0.0
            #dti = float(cMF.perlen[n])
            h5_MM['MM'][tstart_MM:tend_MM,:,:,:]     = MM[:,:,:,:]
            h5_MM['MM_S'][tstart_MM:tend_MM,:,:,:,:] = MM_S[:,:,:,:,:]
            h5_MM['finf'][n,:,:]                     = MM_finf_MF
            h5_MM['ETg'][n,:,:]                      = MM_wel_MF
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
                Rg           Daily gross recharge
    OUTPUTS
            h               Daily water level
     _________________________
    """

    def runSATFLOW(self, Rg, hi, h0, RC, STO):
        h1 = np.zeros([len(Rg)], dtype=np.float)
        h_tmp = (hi*1000.0 + Rg[0]/STO - hi*1000.0/RC)
        h1[0] = h_tmp
        for t in range(1,len(Rg)):
            h_tmp = h1[t-1] + Rg[t]/STO -h1[t-1]/RC
            h1[t] = h_tmp

        h = np.zeros([len(Rg)], dtype=np.float)
        for t in range(0,len(Rg)):
            h_tmp = (h1[t] + h0*1000.0)*0.001
            h[t] = h_tmp

        return h
        del Rg, h1, h_tmp, h

##################
if __name__ == "__main__":
    print '\nWARNING!\nStart MARMITES-MODFLOW models using the script startMARMITES_v3.py\n'

#EOF#