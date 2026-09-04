# -*- coding: utf-8 -*-
"""
MARMITES is a distributed depth-wise lumped-parameter model for
spatio-temporal assessment of water fluxes in the unsaturated zone.
MARMITES is a french word to
design a big cooking pot used by sorcerers for all kinds of experiments!
The main objective of MARMITES development is to partition the rainfall
in the several fluxes from the unsaturated and saturated zone (source
reservoir). It applies the concepts enunciated by Lubcynski (2010):                        ## TO BE UPDATED## 20100823
ET = Es + Ei + ETu + ETg
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
import sys
from decimal import Decimal, ROUND_HALF_EVEN

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
                P           Daily rainfall
                PT           Daily transpiration
                PE           Daily evaporation
    OUTPUTS
            Pe              Effective rainfall
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
        # table 3, values for bare soil
        # dll [cm], y0 [], b [cm^-1], ext_d [cm]
        self.paramEg = {'sand'            : {'dll':16.0,'y0':0.000,'b':0.171, 'ext_d':50.0},
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
                        'sandy loam field_enrico': {'dll':100.0,'y0':0.000,'b':0.013, 'ext_d':475.0},
                        'sandy loam field': {'dll':115.3,'y0':0.023,'b':0.013, 'ext_d':1000.0}
                       }
                     #                        'custom'          : {'dll':100.0,'y0':0.00,'b':0.013, 'ext_d':330.0}
#####################

    def flux(self, cMF, perleni, Pe, PT, PE, Eosurf_max, Zr_elev, VEGarea, HEADSini, TopSoilLay, BotSoilLay, Tl, nsl, Sm, Sfc, Sr, Ks, Ssurf_max, Ssoil_ini, Ssurf_ini, EXF_ini, dgwt, st, i, j, n, kTg_min, kTg_max, kT_f, kT_s, NVEG, LAIveg):

        ##################

        def perc(s_tmp, Sm, Sfc, Ks, perlen):
            '''
            Percolation function
            '''
            global rp_tmp
            if s_tmp <= Sfc:
                rp_tmp = 0.0
            elif s_tmp  > Sfc: # and (s_tmp - Sm) <= 1.0E-5:
                Sg = (s_tmp-Sfc)/(Sm-Sfc) # gravitational water
                if (Ks*Sg*perlen) > (s_tmp-Sfc):
                    rp_tmp = (s_tmp-Sfc)/perlen
                else:
                    rp_tmp = Ks*Sg
            #elif (s_tmp - Sm) > 1.0E-5:
            #    if (Ks) - (Sm-Sfc) > 1.0E-5:
            #        rp_tmp = (Sm-Sfc)/perlen
            #    else:
            #        rp_tmp = Ks
            return rp_tmp

        ##################

        def evp(s_tmp,Sm,Sr, pet, perlen):
            '''
            Actual evapotranspiration function
            '''
            #global evp_tmp
            pet = Decimal(float(pet)).quantize(Decimal('0.0001'), rounding=ROUND_HALF_EVEN)
            if s_tmp <= Sr:
                evp_tmp = 0.0
            elif s_tmp > Sr: # and (s_tmp - Sm) < 1.0E-5:
                Se = (s_tmp - Sr)/(Sm - Sr)
                if pet*Se*perlen > (s_tmp - Sr):
                    evp_tmp = (s_tmp - Sr)/perlen
                else:
                    evp_tmp = pet*Se
            #elif (s_tmp - Sr) > 1.0E-5:
            #    if (pet*perlen-(Sm - Sr)) > 1.0E-5:
            #        evp_tmp = (Sm - Sr)/perlen
            #    else:
            #        evp_tmp = pet
            return evp_tmp

        ##################

        if EXF_ini < 0.0:
            print('WARNING!\nEXFg < 0.0, value %.6f corrected to 0.0.' % EXF_ini)
            EXF_ini = 0.0

        # INITIALIZATION
        SAT = list(np.zeros([nsl], dtype = bool))
        perlen = cMF.perlen[n]
        # Surface
        Ssurf_tmp = Decimal(0.0).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN)
        Ssurf_tmp  += Decimal(float(Pe)).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN) + Decimal(Ssurf_ini).quantize(Decimal('.001'), rounding=ROUND_HALF_EVEN)
        # Soil
        Ssoil_tmp = []
        Rexf_tmp =  []
        for l in range(nsl):
            Rexf_tmp.append(Decimal(0.0).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN))
        #Rexf_tmp1 = list(np.zeros([nsl]))
        for l in range(nsl):
            Ssoil_tmp.append(Decimal(Ssoil_ini[l]).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN))
        # PERCOLATION
#        for l in range(1,nsl):
#            Ssoil_tmp[l] += Decimal(Rp_ini[l-1]*perleni).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN)
        # EXFILTRATION
        Ssoil_tmp[-1] += Decimal(float(EXF_ini*perleni)).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN)

        # SOIL EXF
        llst = list(range(nsl))
        llst.reverse()
        if EXF_ini > 0.0:
            for l in llst:
                if Ssoil_tmp[l] >= Sm[l]*Tl[l]:
                    Rexf_tmp[l] = Ssoil_tmp[l] - Sm[l]*Tl[l]
                    Ssoil_tmp[l] = Sm[l]*Tl[l]
                if l!=0:
                    Ssoil_tmp[l-1] += Rexf_tmp[l]
            Rexf_tmp /= perlen
            # SURFACE
            Ssurf_tmp += Rexf_tmp[0]*perlen                   
        #del Rexf_tmp1
            
        #I
        if Ssurf_tmp > (Sm[0]*Tl[0] - Ssoil_tmp[0]):
            I = Sm[0]*Tl[0] - Ssoil_tmp[0]
            Ssurf_tmp -= Sm[0]*Tl[0] - Ssoil_tmp[0]
        else:
            I = Ssurf_tmp
            Ssurf_tmp = 0.0
        Ssoil_tmp[0] += I
        I /= perlen

        # SAT = True indicates saturation overland flow
        HEADSini_corr = HEADSini * 1.0
        dgwt_corr = dgwt
        if EXF_ini > 0.0:
            for l in llst:
                if Ssoil_tmp[l] == Sm[l]*Tl[l] :
                    HEADSini_corr += float(Tl[l])
                    dgwt_corr -= Tl[l]
                    SAT[l] = True
                else:
                    break
            
        # SURFACE storage, Ro, Eo
        if Ssurf_tmp > Ssurf_max:
            Ro_tmp = (Ssurf_tmp - Ssurf_max)/perlen
            Ssurf_tmp = Ssurf_max
        else:
            Ro_tmp = 0.0
        if Ssurf_tmp > Eosurf_max:
            Eow_tmp = Eosurf_max
            Ssurf_tmp -= Eosurf_max*perlen
        else:
            Eow_tmp = Ssurf_tmp/perlen
            Ssurf_tmp = 0.0

        Rp_tmp = list(np.zeros([nsl]))
        Tsoil_tmpZr = list(np.zeros([nsl,NVEG]))
        Tsoil_tmp = list(np.zeros([nsl]))
        Esoil_tmp = list(np.zeros([nsl]))
        Ssoil_pc_tmp = list(np.zeros([nsl]))

        # soil layers
        for l in range(nsl):
            # Tsoil
            for v in range(NVEG):
                if PT[v] > 0.0:
                    if LAIveg[v] > 0.0:
                        if VEGarea[v] > 0.0:
                            if BotSoilLay[l] > Zr_elev[v]:
                                Tsoil_tmpZr[l][v] = evp(Ssoil_tmp[l], Sm[l]*Tl[l], Sr[l]*Tl[l], PT[v], perlen)
                            elif TopSoilLay[l] > Zr_elev[v] :
                                PTc = PT[v]*(TopSoilLay[l]-np.float32(Zr_elev[v]))/np.float32(Tl[l])
                                Tsoil_tmpZr[l][v] = evp(Ssoil_tmp[l], Sm[l]*Tl[l], Sr[l]*Tl[l], PTc, perlen)
                            Tsoil_tmp[l] += Tsoil_tmpZr[l][v]*VEGarea[v]*0.01
                            PT[v] -= Tsoil_tmpZr[l][v]
                            if PT[v] < 0.0:
                               PT[v] = 0.0 
            Ssoil_tmp[l] -= Decimal(Tsoil_tmp[l]*perlen).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
            #Tsoil_tmp[l] = Decimal(Tsoil_tmp[l]).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
            # Esoil
            if Ssurf_tmp == 0.0:
                if PE > 0.0:
                    Esoil_tmp[l] = evp(Ssoil_tmp[l],Sm[l]*Tl[l], Sr[l]*Tl[l], PE, perlen)
                    Ssoil_tmp[l] -= Decimal(Esoil_tmp[l]*perlen).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
                    PE = Decimal(float(PE)- float(Esoil_tmp[l])).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
            # Rp
            if l < (nsl-1):
                if SAT[l+1] == False:
                    Rp_tmp[l] = perc(Ssoil_tmp[l], Sm[l]*Tl[l], Sfc[l]*Tl[l], Ks[l], perlen)
                    if Rp_tmp[l]*perlen > (Sm[l+1]*Tl[l+1] - Ssoil_tmp[l+1]):
                        Rp_tmp[l] = (Sm[l+1]*Tl[l+1] - Ssoil_tmp[l+1])/perlen
                    Ssoil_tmp[l+1] += Decimal(Rp_tmp[l] * perlen).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
            elif EXF_ini == 0.0:
                Rp_tmp[l] = perc(Ssoil_tmp[l], Sm[l]*Tl[l], Sfc[l]*Tl[l], Ks[l], perlen)
            Ssoil_tmp[l] -= Decimal(Rp_tmp[l]*perlen).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)

        for l in range(nsl):
            Ssoil_pc_tmp[l] = Ssoil_tmp[l]/Tl[l]
            
        sy_tmp = cMF.cPROCESS.float2array(cMF.sy_actual)[cMF.outcropL[i,j]-1,i,j]
        sy_tmp = Decimal(sy_tmp).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
        dgwt_corr_tmp = Decimal(float(dgwt_corr)*0.1).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
        HEADSini_corr_tmp = Decimal(float(HEADSini_corr)*0.1).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
        # GW evaporation Eg, equation 17 of Shah et al 2007, see ref in the __init__
        if cMF.wel_yn == 1:
            if Ssurf_tmp == 0.0:
                if PE > 0.0:
                    y0    = Decimal(self.paramEg[st]['y0']).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
                    b     = Decimal(self.paramEg[st]['b']).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
                    dll   = Decimal(self.paramEg[st]['dll']).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
                    ext_d = Decimal(self.paramEg[st]['ext_d']).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
                    #print 'dgwt_corr_tmp %.2f, dll %.2f, ext_d %.2f' % (dgwt_corr_tmp, dll, ext_d)
                    if dgwt_corr_tmp <= dll:
                        Eg_tmp = PE
                    elif dgwt_corr_tmp < ext_d:
                        Eg_tmp = PE*(y0 + np.exp(-b*(dgwt_corr_tmp-dll)))
                    else:
                        Eg_tmp = Decimal(0.0).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
                    if Eg_tmp>0.0:
                        if (dgwt_corr_tmp + Decimal(0.1).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)*Eg_tmp/sy_tmp) > ext_d:
                            Eg_tmp = Decimal(10.0).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)*(ext_d - dgwt_corr_tmp)*sy_tmp
                            dgwt_corr_tmp = ext_d
                            HEADSini_corr_tmp -= ext_d
                        else:
                            dgwt_corr_tmp += Decimal(0.1).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)*Eg_tmp/sy_tmp
                            HEADSini_corr_tmp -= Decimal(0.1).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)*Eg_tmp/sy_tmp
                    #print 'new dgwt_corr_tmp %.2f, Eg_tmp %.2f, PE %.2f' %(dgwt_corr_tmp, Eg_tmp, PE)
                    #print '------'                    
                else:
                    Eg_tmp = Decimal(0.0).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
            else:
                Eg_tmp = Decimal(0.0).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
        else:
            Eg_tmp = Decimal(0.0).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN)
        Eg_tmp = float(Eg_tmp)
        dgwt_corr_tmp = 10.0*float(dgwt_corr_tmp)
        HEADSini_corr_tmp = 10.0*float(HEADSini_corr_tmp)

        # Groundwater transpiration Tg
        Tg_tmp_Zr = np.zeros([nsl, len(Zr_elev)])
        Tg_tmp = 0.0
        if cMF.wel_yn == 1:
            order = np.asarray(Zr_elev).flatten().argsort()
            #print('Zr_elev', Zr_elev, np.array(Zr_elev)[order])
            #print('NVEG', range(NVEG), np.array(range(NVEG))[order])
            #print('kTg_max', kTg_max, np.array(kTg_max)[order])
            for jj, (Zr_elev_, v, kTg_min_, kTg_max_, kT_f_, kT_s_) in enumerate(zip(np.array(Zr_elev)[order],np.array(list(range(NVEG)))[order], np.array(kTg_min)[order], np.array(kTg_max)[order], np.array(kT_f)[order], np.array(kT_s)[order])):
                #print '----'
                if HEADSini_corr_tmp > Zr_elev_:
                    for l in range(nsl):
                        #print 'veg%d, soil layer %d' %(v,l)
                        if np.isclose(float(Ssoil_pc_tmp[l]), cMF.hnoflo):
                            pass;
                        else:
                            if Ssoil_pc_tmp[l] > Sr[l]:
                                if Ssoil_pc_tmp[l] < Sm[l]:
                                    Ssoil_norm_tmp = (float(Ssoil_pc_tmp[l]) - float(Sr[l]))/(float(Sm[l])-float(Sr[l]))
                                    kTg = float(kTg_max_) - float(kTg_max_-kTg_min_) / (1.0 + np.exp((Ssoil_norm_tmp-float(kT_f_))/float(kT_s_)))
                                elif Ssoil_pc_tmp[l] == Sm[l]:
                                    kTg = float(kTg_max_)
                                elif Ssoil_pc_tmp[l] > Sm[l]:
                                    print("WARNING!\nComputing of Tg: soil moisture higher than porosity!\nSoil moisture = %.4f, phi = %.4f" % (Ssoil_pc_tmp[l], Sm[l]))
                                    kTg = float(kTg_max_)
                            elif Ssoil_pc_tmp[l] == Sr[l]:
                                kTg = float(kTg_min_)
                            elif Ssoil_pc_tmp[l] < Sr[l]:
                                print("WARNING!\nComputing of Tg: soil moisture lower than wilting point!\nSoil moisture = %.4f, WP = %.4f" % (Ssoil_pc_tmp[l], Sr[l]))
                                kTg = float(kTg_min_)
                            Tg_tmp_Zr[l,v] = PT[v]*kTg
                            PT[v] -= Tg_tmp_Zr[l,v]
                            Tg_tmp1 = (Tg_tmp_Zr[l,v]*VEGarea[v]*0.01)
                            if Decimal(Tg_tmp1).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN) > 0.0 and Decimal(float(HEADSini_corr_tmp)-Tg_tmp1/float(sy_tmp)).quantize(Decimal('0.00001'), rounding=ROUND_HALF_EVEN) < Zr_elev_:
                                #print 'WARNING at cell [i=%d,j=%d,L=%d] with NVEG = %d\nTg = %.2f, HEADSini_corr_tmp = %.2f, dgwt_corr_tmp= %.2f, Zr_elev = %.2f, Sy = %.3f' % (i+1,j+1,cMF.outcropL[i,j], v, Tg_tmp1, HEADSini_corr_tmp, dgwt_corr_tmp, Zr_elev_, sy_tmp)
                                Tg_tmp1 = (HEADSini_corr_tmp - float(Zr_elev[v]))*float(sy_tmp)
                                #print 'New Tg = %.2f mm' % Tg_tmp1
#                            else:                        
#                                print 'Tg veg%d = %.2f mm' %(v,Tg_tmp1)
                            Tg_tmp += Tg_tmp1
                            dgwt_corr_tmp += Tg_tmp1/float(sy_tmp)
                            HEADSini_corr_tmp -= Tg_tmp1/float(sy_tmp)
#                else:
#                    print 'veg%d: root too short!' % v
        else:
            Tg_tmp = 0.0

        del sy_tmp, HEADSini_corr_tmp, dgwt_corr_tmp
    
        return Eow_tmp, Ssurf_tmp, Ro_tmp, Rp_tmp, Esoil_tmp, Tsoil_tmp, Ssoil_tmp, Ssoil_pc_tmp, Eg_tmp, Tg_tmp, HEADSini_corr, dgwt_corr, SAT, Rexf_tmp, I

#####################

    def runMMsoil(self, _nsl, _nslmax, _st, _Sm, _Sfc, _Sr, _slprop, _Ssoil_ini, botm_l0, _Ks,
            gridSOIL, gridSOILthick, TopSoil, gridMETEO,
            index, index_S, gridSsurfhmax, gridSsurfw,
            P_veg_zoneSP, Eo_zonesSP, PT_veg_zonesSP, Pe_veg_zonesSP, PE_zonesSP, gridVEGarea,
            LAI_veg_zonesSP, Zr, kTg_min, kTg_max, kT_f, kT_s, NVEG,
            cMF, conv_fact, h5_MF, h5_MM, irr_yn,
            P_irr_zoneSP = [], PT_irr_zonesSP = [], Pe_irr_zoneSP = [],
            crop_irr_SP = [], gridIRR = [],
            Zr_c = [], kTg_min_c = [], kTg_max_c = [], kT_f_c= [], kT_s_c= [],
            verbose = 0, report = None, report_fn = None, stdout = None):

        global dgwt, IRRfield
        h_MF_ini = None
        h_MF_ini_mem = 'slow'
        try:
            h_MF_ini = h5_MF['heads4MM'][:,:,:]
            h_MF_ini_mem = 'fast'
        except:
            print('\nRAM memory too small compared to the size of the heads array -> slow computing.')
        if cMF.uzf_yn == 1:
            exf_MF_ini = None
            exf_MF_ini_mem = 'slow'
            try:
                # TODO this below assume that there is no grid refinement
                exf_MF_ini = h5_MF['exf4MM'][:,:,:]*conv_fact/(cMF.delr[0]*cMF.delc[0])
                exf_MF_ini_mem = 'fast'
            except:
                print('\nRAM memory too small compared to the size of the exfiltration array -> slow computing.')
        tstart_MM = 0
        tstart_MF = 0
        # initial values of SP
        Ssoil_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol,_nslmax])
        #Rp_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol,_nslmax])
        Ssurf_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol])
        Ssurf_ini_tmp_array = np.zeros([cMF.nrow,cMF.ncol])
        MM_perc_MF = np.zeros([cMF.nrow,cMF.ncol], dtype=np.float32)
        MM_wel_MF  = np.zeros([cMF.nrow,cMF.ncol], dtype=np.float32)
        for n in range(cMF.nper):
            if n > 0:
                tstart_MM += cMF.perlen[n-1]
                tstart_MF += cMF.nstp[n-1]
            tend_MM = tstart_MM + cMF.perlen[n]
            tend_MF = tstart_MF + cMF.nstp[n]
            MM         = np.zeros([cMF.perlen[n],cMF.nrow,cMF.ncol,len(index)], dtype=np.float32)
            MM_S       = np.zeros([cMF.perlen[n],cMF.nrow,cMF.ncol,_nslmax,len(index_S)], dtype=np.float32)
            if h_MF_ini_mem == 'slow':
                h_MF_ini = h5_MF['heads4MM'][tstart_MF:tend_MF,:,:]
            if cMF.uzf_yn == 1:
                if exf_MF_ini_mem == 'slow':
                    exf_MF_ini = h5_MF['exf4MM'][tstart_MF:tend_MF,:,:]*conv_fact/(cMF.delr[0]*cMF.delc[0])
            # loop into the grid
            for i in range(cMF.nrow):
                for j in range(cMF.ncol):
                    if cMF.outcropL[i,j] > 0:
                        SOILzone_tmp = gridSOIL[i,j]-1
                        METEOzone_tmp = gridMETEO[i,j]-1
                        slprop = _slprop[SOILzone_tmp]
                        nsl = _nsl[SOILzone_tmp]
                        # thickness of soil layers
                        Tl_ = gridSOILthick[i,j]*slprop*1000.0
                        Tl = []
                        for l in list(range(nsl)):
                            Tl.append(Decimal(Tl_[l]).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN))
                        del Tl_
                        # elevation of top and bottom of soil layers
                        TopSoilLay = np.zeros([nsl], dtype=np.float32)
                        BotSoilLay = np.zeros([nsl], dtype=np.float32)
                        for l in range(nsl):
                            if l==0:
                                TopSoilLay[0] = Decimal(TopSoil[i,j]).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN)
                                BotSoilLay[0] = Decimal(TopSoil[i,j]).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN) - Tl[0]
                            else:
                                TopSoilLay[l] = Decimal(float(BotSoilLay[l-1])).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN)
                                BotSoilLay[l] = Decimal(float(TopSoilLay[l])).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN) - Tl[l]
                        Ssoil_ini_tmp = []
                        if n == 0:
                            for l in range(nsl):
                                Ssoil_ini_tmp.append(Decimal(_Ssoil_ini[SOILzone_tmp][l]).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN))
                            Ssurf_ini_tmp  = 0.0
                            perleni = 1.0
                        else:
                            Ssoil_ini_tmp = []
                            for l in range(nsl):
                                Ssoil_ini_tmp.append(Decimal(float(Ssoil_ini_tmp_array[i,j,l])).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN))
                            Ssurf_ini_tmp    = Decimal(Ssurf_ini_tmp_array[i,j]).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN)
                        if irr_yn == 1:
                            IRRfield = gridIRR[i,j]
                        if  irr_yn == 1 and IRRfield > 0:
                            NVEG_tmp = 1
                            IRRfield -= 1
                            VEGarea_tmp = None
                            LAIveg_tmp = np.ones((NVEG_tmp,tend_MF-tstart_MF), dtype = np.float32)
                            CROP_tmp = crop_irr_SP[IRRfield,tstart_MF:tend_MF]
                            P_tmp = P_irr_zoneSP[METEOzone_tmp,IRRfield,tstart_MF:tend_MF]
                            PT_zonesSP_tmp = np.zeros((NVEG_tmp,tend_MF-tstart_MF), dtype = np.float32)
                            Pe_zonesSP_tmp = np.zeros((NVEG_tmp,tend_MF-tstart_MF), dtype = np.float32)
                            Pe_zonesSP_tmp[0,:] = Pe_irr_zoneSP[METEOzone_tmp,IRRfield,tstart_MF:tend_MF]
                            PT_zonesSP_tmp[0,:] = PT_irr_zonesSP[METEOzone_tmp,IRRfield,tstart_MF:tend_MF]
                            Zr_tmp = Zr_c
                            kTg_min_tmp = kTg_min_c
                            kTg_max_tmp = kTg_max_c
                            kT_f_tmp = kT_f_c
                            kT_s_tmp = kT_s_c
                        else:
                            NVEG_tmp = NVEG
                            CROP_tmp = None
                            P_tmp = P_veg_zoneSP[METEOzone_tmp][tstart_MF:tend_MF]
                            PT_zonesSP_tmp = np.zeros((NVEG_tmp,tend_MF-tstart_MF), dtype = np.float32)
                            Pe_zonesSP_tmp = np.zeros((NVEG_tmp,tend_MF-tstart_MF), dtype = np.float32)
                            LAIveg_tmp = np.zeros((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float32)
                            VEGarea_tmp = np.zeros([NVEG_tmp], dtype = np.float32)
                            Zr_tmp = np.ones((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float32)
                            kTg_min_tmp = np.ones((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float32)
                            kTg_max_tmp = np.ones((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float32)
                            kT_f_tmp = np.ones((NVEG_tmp,tend_MF - tstart_MF), dtype = np.float32)
                            kT_s_tmp = np.ones((NVEG_tmp, tend_MF - tstart_MF), dtype=np.float32)
                            for v in range(NVEG_tmp):
                                PT_zonesSP_tmp[v,:]  = PT_veg_zonesSP[METEOzone_tmp,v,tstart_MF:tend_MF]
                                Pe_zonesSP_tmp[v,:] = Pe_veg_zonesSP[METEOzone_tmp,v,tstart_MF:tend_MF]
                                LAIveg_tmp[v,:]      = LAI_veg_zonesSP[v,tstart_MF:tend_MF]
                                VEGarea_tmp[v]       = gridVEGarea[v,i,j]
                                Zr_tmp[v]            = Zr[v]
                                kTg_min_tmp[v]        = kTg_min[v]
                                kTg_max_tmp[v]        = kTg_max[v]
                                kT_f_tmp[v]          = kT_f[v]
                                kT_s_tmp[v]          = kT_s[v]
                        PE_zonesSP_tmp = PE_zonesSP[METEOzone_tmp,SOILzone_tmp,tstart_MF:tend_MF]
                        Eo_zonesSP_tmp = Eo_zonesSP[METEOzone_tmp][tstart_MF:tend_MF]
                        if h_MF_ini_mem == 'slow':
                            h_MF_ini_tmp   = h_MF_ini[:,i,j]
                        elif h_MF_ini_mem == 'fast':
                            h_MF_ini_tmp   = h_MF_ini[tstart_MF:tend_MF,i,j]
                        if cMF.uzf_yn == 1:
                            if exf_MF_ini_mem == 'slow':
                                exf_MF_ini_tmp = -exf_MF_ini[:,i,j]
                            elif exf_MF_ini_mem == 'fast':
                                exf_MF_ini_tmp = -exf_MF_ini[tstart_MF:tend_MF,i,j]
                        else:
                            exf_MF_ini_tmp = 0.0
                        st         = _st[SOILzone_tmp]
                        Sm         = _Sm[SOILzone_tmp]
                        Sfc        = _Sfc[SOILzone_tmp]
                        Sr         = _Sr[SOILzone_tmp]
                        Ks         = _Ks[SOILzone_tmp]
                        shapeFactor = 1.126847784
                        Ssurf_max  = Decimal(float(np.power(cMF.delr[j],3)*gridSsurfhmax[i,j]*gridSsurfw[i,j]*shapeFactor/np.power(100.0,2)/10.0)).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN)   #1000*1.12*gridSsurfhmax[i,j]*gridSsurfw[i,j]/cMF.delr[j]
                        Eosurf_max = Decimal(float(cMF.delr[j]*gridSsurfw[i,j]*shapeFactor*Eo_zonesSP_tmp/np.power(100.0,2))).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN)  #1.12*gridSsurfw[i,j]/cMF.delr[j]
                        # Output initialisation
                        # PT for the vegetation patchwork
                        PT_tot = np.zeros([len(PT_zonesSP_tmp[0])], dtype = np.float32)
                        # Pe for the vegetation patchwork
                        Pe_tot = np.zeros([len(Pe_zonesSP_tmp[0])], dtype = np.float32)
                        #Soil moisture storage changes
                        dSsoil = np.zeros([nsl], dtype = np.float32)
                        # MASS BALANCE each soil layer
                        MB_l = np.zeros([nsl], dtype = np.float32)

                        # PROCESSING THE WHOLE DATA SET
                        # Preprocessing of PT/PE/INTER and Pe
                        # for different vegetation
                        Zr_elev = []
                        SOILarea = 100.0
                        if CROP_tmp != None:
                            Zr_elev.append(Decimal(float(TopSoilLay[0])).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN) - Decimal(float(Zr_tmp[CROP_tmp-1])*1000.0).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN))
                            if CROP_tmp > 0:
                                VEGarea_tmp = [100.0]
                            else:
                                VEGarea_tmp = [0.0]
                            kTg_min_tmp = [kTg_min_tmp[CROP_tmp-1]]
                            kTg_max_tmp = [kTg_max_tmp[CROP_tmp-1]]
                            kT_f_tmp = [kT_f_tmp[CROP_tmp-1]]
                            kT_s_tmp = [kT_s_tmp[CROP_tmp - 1]]
                        else:
                            for v in range(NVEG_tmp):
                                Zr_elev.append(Decimal(float(TopSoilLay[0])).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN) - Decimal(Zr_tmp[v][0] * 1000.0).quantize(Decimal('.00001'), rounding=ROUND_HALF_EVEN))
                            kTg_min_tmp = kTg_min_tmp[:]
                            kTg_max_tmp = kTg_max_tmp[:]
                            kT_f_tmp = kT_f_tmp[:]
                            kT_s_tmp = kT_s_tmp[:]
                        for v in range(NVEG_tmp):
                            if LAIveg_tmp[v] > 1.0E-5:
                                Pe_tot += Pe_zonesSP_tmp[v]*VEGarea_tmp[v]*0.01
                                PT_tot  += PT_zonesSP_tmp[v]*VEGarea_tmp[v]*0.01
                                SOILarea -= VEGarea_tmp[v]
                        Pe_tot   += P_tmp*SOILarea*0.01
                        INTER_tot  = P_tmp - Pe_tot
                        PE_tot     = PE_zonesSP_tmp*SOILarea*0.01
                        # handle drycell
                        if np.abs(h_MF_ini_tmp - cMF.hdry) < 1.0E-5:
                            HEADSini_drycell = botm_l0[i,j]*1000.0
                        else:
                            HEADSini_drycell = h_MF_ini_tmp*1000.0
                        # dgwt and uzthick
                        if exf_MF_ini_tmp <= 0.0: # or BotSoilLay[nsl-1] > HEADSini_drycell:
                            dgwt = TopSoilLay[0] - HEADSini_drycell
                        elif exf_MF_ini_tmp > 0.0:
                            dgwt = sum(Tl)
                            HEADSini_drycell = BotSoilLay[nsl-1]
                        uzthick = BotSoilLay[nsl-1] - HEADSini_drycell
                        # for the first SP, S_ini is expressed in % and has to be converted in mm
                        if n == 0:
                            for k in range(len(Tl)):
                                Ssoil_ini_tmp[k] = Ssoil_ini_tmp[k] * Tl[k]
                        # MAIN SUB-ROUTINE fluxes
                        Eow_tmp, Ssurf_tmp, Ro_tmp, Rp_tmp, Esoil_tmp, Tsoil_tmp, Ssoil_tmp, Ssoil_pc_tmp, Eg_tmp, Tg_tmp, HEADSini_MM, dgwt_tmp, SAT_tmp, Rexf_tmp, I = self.flux(cMF, perleni, Pe_tot, PT_zonesSP_tmp[:], PE_zonesSP_tmp*SOILarea*0.01, Eosurf_max, Zr_elev, VEGarea_tmp, HEADSini_drycell, TopSoilLay, BotSoilLay, Tl, nsl, Sm, Sfc, Sr, Ks, Ssurf_max, Ssoil_ini_tmp, Ssurf_ini_tmp, exf_MF_ini_tmp, dgwt, st, i, j, n, kTg_min_tmp, kTg_max_tmp, kT_f_tmp, kT_s_tmp, NVEG_tmp, LAIveg_tmp[:])
                        Ssoil_pc_tot = sum(Ssoil_pc_tmp[:])/nsl
                        perc     = Rp_tmp[-1]
                        ETg = Eg_tmp + Tg_tmp
                        dSsurf = (float(Ssurf_tmp) - float(Ssurf_ini_tmp))/cMF.perlen[n]
                        # compute the water mass balance (MB) in the soil zone
                        #Rp_in_MB   = 0.0
                        #Rp_out_MB  = 0.0
                        Esoil_MB   = 0.0
                        Tsoil_MB   = 0.0
                        dSsoil_tot = 0.0
                        for l in range(nsl):
                            Esoil_MB     += float(Esoil_tmp[l])
                            Tsoil_MB     += Tsoil_tmp[l]
                            #if l > 0: #< (nsl-1):
                            #    Rp_in_MB  += Rp_ini_tmp_array[i,j,l-1]*perleni
                            #Rp_out_MB += float(Rp_tmp[l])*cMF.perlen[n]
                            dSsoil[l] = float((Ssoil_tmp[l]-Ssoil_ini_tmp[l])/cMF.perlen[n])
                            dSsoil_tot += dSsoil[l]
                        #dRp_tot = (Rp_in_MB - Rp_out_MB)/cMF.perlen[n]
                        ETsoil_tot = Esoil_MB + Tsoil_MB
                       
                       # MASS BALANCE COMPUTING
                        MBsurf = Pe_tot + float(Rexf_tmp[0]) - (float(Eow_tmp) + float(Ro_tmp) + float(I) + float(dSsurf))
                        if nsl > 1:
                            # surficial soil layer
                            l = 0
                            In = float(I) + float(Rexf_tmp[l+1])
                            Out = float(Rp_tmp[l]) + float(Rexf_tmp[l]) + float(Esoil_tmp[l]) + Tsoil_tmp[l] + dSsoil[l]
                            MB_l[l] = In - Out
                            del In, Out
                            # intermediate soil layers
                            llst = list(range(1,nsl-1))
                            for l in llst:
                                # Rexf[l] already added to Ssoil ???
                                In =  float(Rp_tmp[l-1]) + float(Rexf_tmp[l+1]) #
                                Out = float(Rp_tmp[l]) + float(Rexf_tmp[l]) + float(Esoil_tmp[l]) + Tsoil_tmp[l] + dSsoil[l]
                                MB_l[l] = In - Out
                                del In, Out
                            # last soil layer
                            l = nsl-1
                            In = float(Rp_tmp[l-1]) + float(exf_MF_ini_tmp/cMF.perlen[n])  #*perleni/cMF.perlen[n-1])
                            Out = float(Rp_tmp[l]) + float(Rexf_tmp[l]) + float(Esoil_tmp[l]) + Tsoil_tmp[l] + dSsoil[l]   #float(sum(Rexf_tmp))
                            MB_l[l] = In - Out
                            del In, Out
                        else:
                            # only one soil layer
                            l = 0
                            In = float(I) + float(exf_MF_ini_tmp/cMF.perlen[n])
                            Out = float(Rp_tmp[l]) + float(Rexf_tmp[l]) + float(Esoil_tmp[l]) + float(Tsoil_tmp[l]) +  dSsoil[l]
                            MB_l[l] = In - Out
                            del In, Out
                            # last soil layer
                        # total mass balance for the soil
                        In = float(I) + float(exf_MF_ini_tmp/cMF.perlen[n])
                        Out = float(Rexf_tmp[0]) + Esoil_MB + Tsoil_MB + dSsoil_tot + float(Rp_tmp[-1])
                        MB = In - Out
                        del In, Out
                      
                        # export list
                        # indexes of the HDF5 output arrays
                        # index_MM = {'iP':0, 'iPT':1, 'iPE':2, 'iPe':3, 'iSsurf':4, 'iRo':5, 'iEXFg':6, 'iEow':7, 'iMB':8, 'iEi':9, 'iEo':10, 'iEg':11, 'iTg':12, 'idSsurf':13, 'iETg':14, 'iETsoil':15, 'iSsoil_pc':16, 'idSsoil':17, 'iperc':18, 'ihcorr':19, 'idgwt':20, 'iuzthick':21, 'iI':22, 'iMBsurf':23}
                        MM_tmp = [P_tmp, PT_tot, PE_tot, Pe_tot, np.float32(Ssurf_tmp), np.float32(Ro_tmp), exf_MF_ini_tmp, np.float32(Eow_tmp), np.float32(MB), np.float32(INTER_tot), Eo_zonesSP_tmp, np.float32(Eg_tmp), np.float32(Tg_tmp), np.float32(dSsurf), np.float32(ETg), np.float32(ETsoil_tot), np.float32(Ssoil_pc_tot), np.float32(dSsoil_tot), np.float32(perc), np.float32(HEADSini_MM)*0.001, np.float32(-dgwt_tmp)*0.001, uzthick*0.001, np.float32(I), MBsurf]
                        # index_MM_soil = {'iEsoil':0, 'iTsoil':1,'iSsoil_pc':2, 'iRsoil':3, 'iExf':4, 'idSsoil':5, 'iSsoil':6, 'iSAT':7, 'iMB':8}
                        MM_S_tmp = np.zeros([nsl,len(index_S)], dtype = np.float32)
                        for l in range(nsl):
                            MM_S_tmp[l,:] = [np.float32(Esoil_tmp[l]), np.float32(Tsoil_tmp[l]), np.float32(Ssoil_pc_tmp[l]), Rp_tmp[l], np.float32(Rexf_tmp[l]), dSsoil[l], np.float32(Ssoil_tmp[l]), SAT_tmp[l], MB_l[l]]

                        # fill MF and MM output arrays
                        for k in range(len(index)):
                            MM[:,i,j,k] = MM_tmp[k]
                            #h5_MM['MM'][tstart_MM:tend_MM,i,j,k] = MM_tmp[k]
                        for k in range(len(index_S)):
                            for l in range(nsl):
                                MM_S[:,i,j,l,k] = MM_S_tmp[l,k]
                                #h5_MM['MM_S'][tstart_MM:tend_MM,i,j,l,k] = MM_S_tmp[l,k]
                        # # Volumetric recharge rate UZF1 package
                        MM_perc_MF[i,j] = MM_S_tmp[nsl-1,index_S.get(b'iRsoil')]/conv_fact
                        #h5_MM['perc'][n,i,j] = MM_S_tmp[nsl - 1, index_S.get(b'iRsoil')] / conv_fact
                        # Volumetric recharge rate WEL package
                        MM_wel_MF[i,j] = MM_tmp[index.get(b'iETg')]/conv_fact
                        #h5_MM['ETg'][n,i,j] = MM_tmp[index.get(b'iETg')] / conv_fact
                        del MM_tmp, MM_S_tmp
                        # setting initial conditions for the next SP
                        Ssoil_ini_tmp_array[i,j,:]  = MM_S[cMF.perlen[n]-1,i,j,:,index_S.get(b'iSsoil')]
                        #Rp_ini_tmp_array[i,j,:]     = MM_S[cMF.perlen[n]-1,i,j,:,index_S.get(b'iRsoil')]
                        Ssurf_ini_tmp_array[i,j]    = MM[cMF.perlen[n]-1,i,j,index.get(b'iSsurf')]
                        #Ssoil_ini_tmp_array[i,j,:]  = h5_MM['MM_S'][cMF.perlen[n]-1,i,j,:,index_S.get(b'iSsoil')]
                        #Rp_ini_tmp_array[i,j,:]     = h5_MM['MM_S'][cMF.perlen[n]-1,i,j,:,index_S.get(b'iRsoil')]
                        #Ssurf_ini_tmp_array[i,j]    = h5_MM['MM'][cMF.perlen[n]-1,i,j,index.get(b'iSsurf')]
                    else:
                        if cMF.perlen[n]>1:
                            MM[:,i,j,:] = cMF.hnoflo
                            MM_S[:,i,j,:,:] = cMF.hnoflo
                            #h5_MM['MM'][tstart_MM:tend_MM,i,j,:] = cMF.hnoflo
                            #h5_MM['MM_S'][tstart_MM:tend_MM,i,j,:,:] = cMF.hnoflo
                        else:
                            MM[:,i,j,:] = cMF.hnoflo
                            MM_S[:,i,j,:,:] = cMF.hnoflo
                            #h5_MM['MM'][tstart_MM:tend_MM,i,j,:] = cMF.hnoflo
                            #h5_MM['MM_S'][tstart_MM:tend_MM,i,j,:,:] = cMF.hnoflo
                        MM_perc_MF[i,j] = 0.0
                        MM_wel_MF[i,j] = 0.0
                        #h5_MM['perc'][n, i, j] = 0.0
                        #h5_MM['ETg'][n, i, j] = 0.0
                        #dti = float(cMF.perlen[n])
            h5_MM['MM'][tstart_MM:tend_MM,:,:,:]     = MM[:,:,:,:]
            h5_MM['MM_S'][tstart_MM:tend_MM,:,:,:,:] = MM_S[:,:,:,:,:]
            h5_MM['perc'][n,:,:]                     = MM_perc_MF
            h5_MM['ETg'][n,:,:]                      = MM_wel_MF
            if (n/100.0).is_integer() and n>0:
                print("Processed data up to stress period %d of %d" % (n,cMF.nper))
                if verbose == 0:
                    sys.stdout = stdout
                    report.close()
                    stdout = sys.stdout
                    report = open(report_fn, 'a')
                    sys.stdout = report
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
        h1 = np.zeros([len(Rg)], dtype=np.float32)
        h_tmp = (hi*1000.0 + Rg[0]/STO - hi*1000.0/RC)
        h1[0] = h_tmp
        for t in range(1,len(Rg)):
            h_tmp = h1[t-1] + Rg[t]/STO -h1[t-1]/RC
            h1[t] = h_tmp

        h = np.zeros([len(Rg)], dtype=np.float32)
        for t in range(0,len(Rg)):
            h_tmp = (h1[t] + h0*1000.0)*0.001
            h[t] = h_tmp

        return h
        del Rg, h1, h_tmp, h

##################
if __name__ == "__main__":
    print('\nWARNING!\nStart MARMITES-MODFLOW models using the script startMARMITES_v3.py\n')

#EOF#