# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:        PET_PM_FAO56
# Purpose:
#
# Author:      frances08512
#
# Created:     25-11-2010
# Copyright:   (c) frances08512 2010
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__author__ = "Alain P. Francés <frances08512@itc.nl>"
__version__ = "0.3"
__date__ = "2012"

# Make a help entry for this library
def PET_PM_FAO56():
    """ Script to compute potential evapo(transpi)ration PET based on hourly measurements of
        meteorological data.
        The main equation is based on equation 3 of FAO 56:
        lambda.ET=[DELTA(Rns-Rnl-G)+(ro_a.Cp/r_a).(e_s-e_a)]/[DELTA+gama.(1+r_s/r_a)]
        The FAO 56 terminology and units are followed in this script.

        PET is defined here following Gieske 2003 (pag 67):
            PET: the maximum possible evaporation according to prevailing atmospheric
            conditions and vegetative properties. The land surface in question (can
            be any part of the landscape that contains a certain fraction of vegetation)
            should be well supplied by water such that soil moisture forms no limitation
            in the stomatal aperture. [...] the biophysical properties of a potentially
            evaporating vegetation are spatially and temporally variables.

        References:

        Allen, R.G., Pereira, L.S., Raes, D. and Smith, M., 1998.
        Crop evapotranspiration : guidelines for computing crop water requirements.
        FAO irrigation and drainage paper, 56. FAO, Rome, 300 pp.

        Gieske, A.S.M., 2003.
        Operational solutions of actual evapotranspiration.
        In: Understanding water in a dry environment : hydrological processes in arid and semi arid zones
        Ed. by I. Simmers. Rotterdam : Balkema, 2003.
        ISBN 9058096181 (IAH International Contributions to Hydrogeology ; 23) pp. 65-114.

        Dingman, S.L., 2002.
        Physical hydrology.
        Prentice Hall, Upper Saddle River, 646 pp.

        surface resistance bare soil

        van de Griend, A.A. and Owe, M., 1994.
        Bare soil surface resistance to evaporation by vapor diffusion under semiarid conditions.
        Water Resources Research, 30(2): 181-188.

        Liu www.hydrol-earth-syst-sci.net/11/769/2007/
        """

    print 'A Python script to compute potential evaporation'
    print 'from hourly meteorological data of different land covers'
    print '(vegetation, open water, bare soil)'
    print 'using Penman-Monteith equation.'
    print 'Type PET.PM() for help'
    print 'Author: ',__author__
    print 'Version: ',__version__
    print 'Date: ',__date__
    return

import numpy as np
import tempfile
import matplotlib as mpl

'''
    ================================================================
    PM function
    ================================================================
'''

def process(cUTIL, datenum = np.array([]), datenum_d = np.array([]), J = np.array([]),  time = np.array([])\
              ,pathMMsurf = tempfile.gettempdir()\
              ,RF = np.array([]), IRR = None, Ta = np.array([]), RHa = np.array([]) \
              ,Pa = np.array([]), u_z_m = np.array([]), Rs = np.array([]) \
              ,phi = 0.0, Lm = 0.0 , Z = 0.0, Lz = 0.0, FC = 0.0, z_m=2.0, z_h=2.0 \
              ,NVEG = 1, VegType = np.array([]), S_w_v = np.array([]), C_leaf_star_v = np.array([]) \
              ,alfa_v = np.array([]), f_s_v = np.array([]), LAI_v = np.array([]), h_v = np.array([])\
              ,NSOIL = 0, SoilType = [], por = np.array([]), fc = np.array([])\
              ,alfa_s = np.array([])
              ,alfa_w = 0.06\
              ,NFIELD = 0, alfa_f = np.array([]), f_s_f = np.array([]), LAI_f = np.array([])\
              ,C_leaf_star_f = np.array([]), h_f = np.array([]), S_w_f = np.array([])
              ):

    '''
    Function to calculate hourly Penman-Monteith evapo(transpi)ration (mm/day).

    Input:
    -------------------------------
        METEO TIME SERIES
        datenum: date in numerical format [python date format]
        time: time [string in format 'hh:mm']
        RF: hourly rainfall [mm]
        Ta: air temperature measured at 2m [ºC]
        RHa: relative air humidity [%] measured at heigth z_h [m]
        Pa: air pressure [kPa]
        u_z_m: windspeed measured at heigth z_m [m.s-1]
        Rs: incoming solar (=shortwave) radiation [MJ.m-2.hour-1]

        METEO PARAMETERS
        phi: latitude of the meteo station [degres, >0 hemisphere N]
        Lm: longitude of the meteo station [degres west from Greenwich]
        Z: altitude of the station above sea level [m]
        Lz: longitude of the center of the local time zone [degres west of Greenwich]
        z_m: heigth of u_z_m measurement [m]
        z_h: heigth of humidity measurement [m]

# TO BE UPDATED, see startMARMITESsurface.py
        VEGETATION PARAMETERS
        VegType: name of the vegetation type [string, 1 word, no space allowed]
        h: heigth of plant [m]
        C_leaf_star: maximum leaf conductance [m.s-1]
        LAI_d: leaf area index dry season [m2.m-2]
        LAI_w: leaf area index wet season [m2.m-2]
        S_d: canopy capacity dry season [mm]
        S_w: canopy capacity wet season [mm]
        f_s_vd: shelter factor dry season []
        f_s_vw: shelter factor wet season []
        alfa_vd: vegetation albedo in dry season
        alfa_vw: vegetation albedo in wet season
        J_vd: starting julian day of the dry season [int 1-365]
        J_vw: starting julian day of the wet season [int 1-365]

# TO BE UPDATED, see startMARMITESsurface.py
        SOIL PARAMETERS
        NSOIL: number of soil types
        to repeat NSOIL times in a same line
        por: surface (1cm) soil porosity [m3.m-3]
        fc: surface (1cm) soil field capacity [m3.m-3]
        alfa_sd: soil albedo in dry season
        alfa_sw: soil albedo in wet season
        J_sd: starting julian day of the dry season [int 1-365]
        J_sw: starting julian day of the wet season [int 1-365]
        TRANS_sdw: transition period between dry and wet season [days]
        TRANS_swd: transition period between wet and dry season [days]

        WATER PARAMETERS
        alfa_w: water albedo

    Output:
    ---------
        PT_PM_VEG     : hourly PT [mm/day] of NVEG type(s) of vegetation covers
        PE_PM_SOIL    : hourly PE [mm/day] of NSOIL type(s) of soil types
        E0            : hourly PE [mm/day] of open water
        NOTE THAT OUTPUTS ARE IN SOLAR TIME

    '''

###########################################

    # Set constants
    sigma = 2.043E-10 # Stefan Boltzmann constant MJ.m-2.h-1 FAO56 pag 74
    lambdav = 2.45    # latent heat of evaporation [MJ.kg-1] FAO56 pag 31
                      # Gieske 2003 pag 74 Eq33/Dingman 2002
                      # lambda=2.501-2.361E-3*t, with t temperature evaporative surface (ºC)
                      # see script Lambda_function_t.py
    Gsc = 0.082       # solar constant [MJ.m-2.min-1] FAO56 pag 47 Eq28
    eps = 0.622       # ratio molecular weigth of vapour/dry air FAO56 p26 BOX6
    R = 0.287         # specific gas [kJ.kg-1.K-1]    FAO56 p26 box6
    Cp = 1.013E-3     # specific heat at cte pressure [MJ.kg-1.ªC-1] FAO56 p26 box6
    k = 0.41          # karman's cte   []  FAO 56 Eq4

    # compute DELTA - SLOPE OF SATURATION VAPOUR PRESSURE CURVE
    # [kPa.ªC-1]
    # FAO56 pag 37 Eq13
    DELTA = 4098*(0.6108*np.exp(17.27*Ta/(Ta+237.3)))/pow(Ta+237.3,2)
    print "\nDELTA computed!"

    # compute dr - inverse distance to the sun
    # FAO56 pag47 Eq23
    dr = 1.0+0.033*np.cos(2.0*np.pi*J/365.0)

    # compute delta - solar declination
    # [rad]
    # FAO56 pag47 Eq24
    delta = 0.409*np.sin(2*np.pi*J/365.0-1.39)

    # compute Sc - seasonnal correction of solar time
    # [hour]
    # FAO56 pag47 Eq32
    Sc = []
    b = 2.0*np.pi*(J-81.0)/364.0    # Eq 34
    Sc = 0.1645*np.sin(2*b) - 0.1255*np.cos(b) - 0.025*np.sin(b)

    # compute w - solar time angle at the midpoint of the period (time)
    # [rad]
    # FAO56 pag48 Eq31
    w = []
    for j in range(len(J)):
        t = -0.5 + FC + float(time[j].split(':')[0]) + \
                  float(time[j].split(':')[1])/60.0
        w.append((np.pi/12)*((t+0.06667*(Lz-Lm)+Sc[j])-12))
    w = np.asarray(w)

    # compute w1 - solar time angle at the beginning of the period (time)
    # [rad]
    # FAO56 pag47 Eq29
    tl = 1  # hourly data
    w1 = (w - np.pi*tl/24.0)

    # compute w2 - solar time angle at the end of the period (time + 1h)
    # [rad]
    # FAO56 pag47 Eq30
    w2 = w + np.pi*tl/24.0

    # compute ws - sunset hour angle
    # [rad]
    # FAO56 pag47 Eq25
    ws = np.arccos(-np.tan(phi*np.pi/180.0)*np.tan(delta))

    # compute Ra - extraterrestrial radiation
    # [MJ.m-2.hour-1]
    # FAO56 pag47 Eq28
    Ra = []
    #Ra_Watts = []
    for j in range(len(J)):
        if (w1[j] < -ws[j] or w2[j] > -ws[j]):
            if w1[j] < -ws[j] : w1[j] = -ws[j]
            if w2[j] < -ws[j] : w2[j] = -ws[j]
        if (w1[j] < ws[j] or w2[j] > ws[j]):
            if w1[j] >  ws[j] : w1[j] =  ws[j]
            if w2[j] >  ws[j] : w2[j] =  ws[j]
        Ra.append((12*60/np.pi)*Gsc*dr[j]* \
              ((w2[j]-w1[j])*np.sin(phi*np.pi/180)*np.sin(delta[j]) + \
               np.cos(phi*np.pi/180)*np.cos(delta[j])*(np.sin(w2[j])-np.sin(w1[j]))))
        #Ra_Watts.append(Ra[j]*24/0.08864)
    Ra = np.asarray(Ra)

    # compute Rs0 - clear-sky solar (shortwave) radiation
    # [MJ.m-2.hour-1]
    # FAO56 pag51 Eq37
    Rs0 = (0.75+2E-5*Z)*Ra
#    Rs0_Watts = Rs0*24.0/0.08864

    # correcting Rs measurement values
    badvalues = 0
    Rs_corr = []
    for j in range(len(J)):
        if Rs[j]>Rs0[j]:
            Rs_corr.append(Rs0[j])
            badvalues = badvalues + 1
        elif Rs[j]<0.0:
            Rs_corr.append(0.0)
        else:
            Rs_corr.append(Rs[j])
    if badvalues > 0:
        print "\n" + '%.1f' %(100.0*badvalues/len(J)) + " % (" + str(badvalues) + '/' + str(len(J)) +') shortwave radiation measurements that were higher than the clear-sky one were corrected.'
        if badvalues/len(J)>0.15:
            print '\nCalibrate your sensors!!!'
    Rs_corr = np.asarray(Rs_corr)

    # compute Rns - NET SHORTWAVE RADIATION
    # [MJ.m-2.hour-1]
    # FAO56 pag51 Eq37
    # for each type of vegetation, crop and soil (albedo dependent)
    Rns_VEG = np.zeros((NVEG, len(J)), dtype = float)
    for v in range(NVEG):
        Rns_VEG[v] = (1 - alfa_v[v])*Rs_corr
    if IRR <> None:
        Rns_FIELD = np.zeros((NFIELD, len(J)), dtype = float)
        for f in range(NFIELD):
            Rns_FIELD[f] = (1 - alfa_f[f])*Rs_corr
    Rns_SOIL = np.zeros((NSOIL, len(J)), dtype = float)
    for s in range(NSOIL):
        Rns_SOIL[s] = (1 - alfa_s[s])*Rs_corr
    Rns_WATER = (1-alfa_w)*Rs_corr
    print "\nRns computed!"

    # compute e0_Ta - saturation vapour pressure at actual air temperature
    # [kPa]
    # FAO56 pag36 Eq11
    e0_Ta = 0.6108*np.exp(17.27*Ta/(Ta+237.3))

    # es - MEAN SATURATION VAPOUR PRESSURE
    # [kPa]
    # FAO56 pag74 eq53
    # es = e0_Ta
    print "\ne_s (e0_Ta) computed!"

    # correcting RHa measurement values
    badhighvalues = 0
    badlowvalues = 0
    for j in range(len(J)):
        if RHa[j]>100.0:
            RHa[j] = 100.0
            badhighvalues = badhighvalues + 1
        elif RHa[j]<25.0:
            RHa[j] = 25.0  #observed value in calibrated sensors in Sardon
            badlowvalues = badlowvalues + 1
    if badhighvalues > 0:
        print "\n" + '%.1f' % (100*badhighvalues/len(J)) + " % (" + str(badhighvalues) + "/" + str(len(J)) + ") relative humidity measurements > 100% were fixed to 100%."
        if badvalues/len(J)>0.15:
            print '\nCalibrate your sensors!!!'
    if badlowvalues > 0:
        print "\n" + '%.1f' % (100*badlowvalues/len(J)) + " % (" + str(badlowvalues) + "/" + str(len(J)) + ") relative humidity measurements < 25% were fixed to 25%."
        if badvalues/len(J)>0.15:
            print '\nCalibrate your sensors!!!'

    # compute e_a - ACTUAL VAPOUR PRESSURE
    # [kPa]
    # FAO56 pag74 Eq54
    e_a = e0_Ta*RHa/100
    print "\ne_a computed!"

    # compute Rnl - NET LONGWAVE RADIATION
    # [MJ.m-2.hour-1]
    # FAO56 pag51 Eq37 and pag74 of hourly computing
    Rnl = []
    r = []
    i = 0
    for j in range(len(J)):
        if (ws[j] - 0.52) <= w[j] <= (ws[j] - 0.10):  # FAO56: (ws[j] - 0.79) <= w[j] <= (ws[j] - 0.52)
            i = 1
            if Rs0[j] > 0:
                if Rs_corr[j]/Rs0[j] > 0.3:
                    r_sunset = Rs_corr[j]/Rs0[j]
                else:
                    r_sunset = 0.3
            else:
                r_sunset = 0.75  #see FAO56 pag75
        if ((ws[j] - 0.10) < w[j] or w[j] <= (-ws[j]+ 0.10)):
            if i>0:
                r.append(r_sunset)
            else:
                r.append(0.75) #see FAO56 pag75
        else:
            r.append(Rs_corr[j]/Rs0[j])
        Rnl.append(sigma*pow(Ta[j] + 273.16,4)*(0.34-0.14*np.sqrt(e_a[j]))*(1.35*r[j]-0.35))
    Rnl = np.asarray(Rnl)
    r = np.asarray(r)
##        if Rnl[j]<0:
##            r=0.8
##            Rnl[j] = sigma*pow(Ta[j] + 273.16,4)*(0.34-0.14*np.sqrt(e_a[j]))*(1.35*r-0.35)

    print "\nRnl computed!"

    # compute G - SOIL HEAT FLUX
    # [MJ.m-2.hour-1]
    # FAO56 pag55 Eq45 and 46
    G_VEG = np.zeros((NVEG, len(J)), dtype = float)
    for v in range(NVEG):
        G_v = []
        Rn = Rns_VEG[v] - Rnl
        for j in range(len(J)):
            if w[j]<-ws[j] or w[j]>ws[j]:
                G_v.append(0.5*Rn[j])
            else:
                G_v.append(0.1*Rn[j])
        G_VEG[v] = np.asarray(G_v)
        del Rn, G_v
    if IRR <> None:
        G_FIELD = np.zeros((NFIELD, len(J)), dtype = float)
        for f in range(NFIELD):
            G_f = []
            Rn = Rns_FIELD[f] - Rnl
            for j in range(len(J)):
                if w[j]<-ws[j] or w[j]>ws[j]:
                    G_f.append(0.5*Rn[j])
                else:
                    G_f.append(0.1*Rn[j])
            G_FIELD[f] = np.asarray(G_f)
            del Rn, G_f
    if NSOIL>0:
        G_SOIL = np.zeros((NSOIL, len(J)), dtype = float)
        for s in range(NSOIL):
            G_s = []
            Rn = Rns_SOIL[s] - Rnl
            for j in range(len(J)):
                if w[j]<-ws[j] or w[j]>ws[j]:
                    G_s.append(0.5*Rn[j])
                else:
                    G_s.append(0.1*Rn[j])
            G_SOIL[s] = np.asarray(G_s)
            del Rn, G_s
    G_w = []
    Rn = Rns_WATER - Rnl
    for j in range(len(J)):
        if w[j]<-ws[j] or w[j]>ws[j]:
            G_w.append(0.5*Rn[j])
        else:
            G_w.append(0.1*Rn[j])
    G_WATER = np.asarray(G_w)
    del Rn, G_w
    print "\nG computed!"
##    G_VEG_Watts = []
##    for j in range(len(J)):
##        G_VEG_Watts.append(G_VEG[0][j]*24/0.08864)

    # ro_a - MEAN AIR DENSITY AT CTE PRESSURE
    # [kg.m-3]
    # FAO56 pag26 box6
    ro_a = Pa/(R*1.01*(Ta+273.16))
    print "\nro_a computed!"

    # gama - PSYCHROMETRIC CONSTANT
    # [kPa.ªC-1]
    # FAO56 pag31 eq8
    gama = Cp*Pa/(eps*lambdav)
    print "\ngama computed!"

    # r_s - SURFACE RESISTANCE
    # [s.m-1]
    # VEG: Dingman pag 208 (canopy conductance) (equivalent to FAO56 pag21 Eq5)
    r_s_VEG = []
    f_k = []
    DELTArho_v = []
    f_rho = []
    f_T = []
    if NVEG>0:
        for v in range(NVEG):
            r_s_VEG.append([])
            for j in range(len(J)):
                if Rs_corr[j]<0:
                    Rs_tmp = 0
                elif Rs_corr[j]>86.5/24:
                    Rs_tmp = 86.5
                else:
                    Rs_tmp = Rs[j]*24
                f_k.append(12.78*Rs_tmp/(11.57*Rs_tmp+104.4))
                DELTArho_v.append(2.17*(e0_Ta[j]-e_a[j])/(Ta[j]+273.16))
                if DELTArho_v[j]<0:
                    f_rho.append(1)
                elif DELTArho_v[j]>0.01152:
                    f_rho.append(0.233)
                else:
                    f_rho.append(1 - 66.6*DELTArho_v[j])
                if Ta[j]<0:
                    f_T.append(0.0)
                elif Ta[j]>40.0:
                    f_T.append(0.0)
                else:
                    f_T.append(Ta[j]*pow(40-Ta[j],1.18)/691)
                if v == 0:
                    C_leaf = C_leaf_star_v[v]
                else:
                    C_leaf = C_leaf_star_v[v] * f_k[j] * f_rho[j] * f_T[j]
                f_temp = f_s_v[v][j]*LAI_v[v][j]*C_leaf
                if f_temp == 0.0:
                    r_s_VEG[v].append(1.0E6)
                else:
                    r_s_VEG[v].append(1.0/f_temp)
        r_s_VEG = np.asarray(r_s_VEG)
    del f_k, DELTArho_v, f_rho, f_T
    # FIELD/CROP: Dingman pag 208 (canopy conductance) (equivalent to FAO56 pag21 Eq5)
    if IRR <> None:
        r_s_FIELD = []
        f_k = []
        DELTArho_f = []
        f_rho = []
        f_T = []
        for f in range(NFIELD):
            r_s_FIELD.append([])
            for j in range(len(J)):
                if Rs_corr[j]<0.0:
                    Rs_tmp = 0.0
                elif Rs_corr[j]>86.5/24.0:
                    Rs_tmp = 86.5
                else:
                    Rs_tmp = Rs[j]*24.0
                f_k.append(12.78*Rs_tmp/(11.57*Rs_tmp+104.4))
                DELTArho_f.append(2.17*(e0_Ta[j]-e_a[j])/(Ta[j]+273.16))
                if DELTArho_f[j]<0.0:
                    f_rho.append(1.0)
                elif DELTArho_f[j]>0.01152:
                    f_rho.append(0.233)
                else:
                    f_rho.append(1.0 - 66.6*DELTArho_f[j])
                if Ta[j]<0:
                    f_T.append(0.0)
                elif Ta[j]>40.0:
                    f_T.append(0.0)
                else:
                    f_T.append(Ta[j]*pow(40-Ta[j],1.18)/691.0)
                C_leaf = C_leaf_star_f[f][j] * f_k[j] * f_rho[j] * f_T[j]
                f_temp = f_s_f[f][j]*LAI_f[f][j]*C_leaf
                if f_temp == 0.0:
                    r_s_FIELD[f].append(1.0E6)
                else:
                    r_s_FIELD[f].append(1.0/f_temp)
        del f_k, DELTArho_f, f_rho, f_T
        r_s_FIELD = np.asarray(r_s_FIELD)
     # SOIL: equation 20 of van de Griend and Owe, 1994
    r_s_SOIL = []
    if NSOIL>0:
        for s in range(NSOIL):
            r_s_SOIL.append(10.0*np.exp(0.3563*100.0*(fc[s]-por[s])))
        print "\nr_s computed!"
    r_s_SOIL = np.asarray(r_s_SOIL)

    # correction windspeed measurement and scaling at h+2m
    # [m.s-1]
    # FAO56 pag56 eq47
#    u_z_m = np.where(u_z_m <= 0.0, 1E-9, u_z_m)
    u_2 = np.where(z_m <> 2.0,u_z_m*4.87/(np.log(67.8*z_m-5.42)),u_z_m)
    u_2 = np.where(u_2 <= 0.0, 1E-9, u_2)
    u_hplus2_v = np.zeros((NVEG, len(J)), dtype = float)
    u_hplus2_v = u_2 * (np.log(67.8*(h_v+2.0)-5.42)) / 4.87
    # FIELD/CROP
    if IRR <> None:
        u_hplus2_f = np.zeros((NFIELD, len(J)), dtype = float)
        u_hplus2_f = u_2*(np.log(67.8*(h_f+2.0)-5.42)) / 4.87

    # r_a - AERODYNAMIC RESISTANCE
    # [s.m-1]
    r_a_VEG = np.zeros((NVEG, len(J)), dtype = float)
    if NVEG>0:
        for v in range(NVEG):
            if v == 0:
                # FAO56 pag20 eq4- (d - zero displacement plane, z_0m - roughness length momentum transfer, z_0h - roughness length heat and vapour transfer, [m], FAO56 pag21 BOX4
                r_a_VEG[0] = np.log((2-(2*h_v[0]/3))/(0.123*h_v[0]))*np.log((2-(2*h_v[0]/3))/(0.0123*h_v[0]))/(pow(k,2)*u_2)
            else:
                # DINGMAN pag 296
                r_a_VEG[v] = pow(np.log((h_v[v]+2-(0.7*h_v[v]))/(0.1*h_v[v])),2)/(pow(k,2)*u_hplus2_v[v])
    # FIELD/CROP
    if IRR <> None:
        # DINGMAN pag 296
        r_a_FIELD = np.where(LAI_f > 0.0, pow(np.log((h_f+2-(0.7*h_f))/(0.1*h_f)),2)/(pow(k,2)*u_hplus2_f),0.0)
    # r_a of SOIL
    # Liu www.hydrol-earth-syst-sci.net/11/769/2007/
    # only function of ws, it is assumed that roughness are the same for any type of soil
    if NSOIL > 0:
        r_a_SOIL = np.log((2.0)/0.0058)*np.log(2.0/0.0058)/(pow(k,2)*u_2)
##    # BavelHillel_PotActualEvaporationBareSoilSurface_1976
##    # SOIL ra = ([ln(2.0/Z0)]^2)/(0.16*U2) with Z0 = 0.01
##    r_a_SOIL = []
##    Z0 = 0.01
##    if NSOIL > 0:
##        for j in range(len(J)):
##            r_a_SOIL.append(pow(np.log(2.0/Z0),2)/(0.16*u_2[j]))
    print "\nr_a computed!"

    # equation variables of outputs
    outputVAR = [DELTA,gama, Rs, Rs_corr, Rs0, Rnl, r]

    # PET, PE and E0 computing
    # WARNING: note the conditions above of >0.0 and <E0 to correct the values
    # and avoid instabilities in the further computing of MARMITESunsat
    # One should analyse the graph of the variables to ensure that PET, PE and E0
    # are properly computed

    # E0 - open water
    # computed by Penman equation in Gieske 2003 pag 94 eq63 and 64
    E0 = []
    value = (DELTA*(Rns_WATER-Rnl-G_WATER) + gama*0.26*(1+0.54*u_2)*(e0_Ta-e_a))/(lambdav*(DELTA + gama[j]))
    E0 = np.where(value >= 0.0, value, 0.0)
    del value
    print "\nE0 of open water computed!"
    # PT/PE - Penman-Montheith
    # mm.hour-1
    # FAO56 pag19 eq3
    # VEG
    PT_PM_VEG = []
    Erf_VEG = []  # Evaporation during RF
    if NVEG>0:
        for v in range(NVEG):
            PT_PM_VEG_v = []
            Erf_v = []
            for j in range(len(J)):
                if r_a_VEG[v][j] > 0.0:
                    value = (DELTA[j]*(Rns_VEG[v][j]-Rnl[j]-G_VEG[v][j])+3600*ro_a[j]*Cp*(e0_Ta[j]-e_a[j])/r_a_VEG[v][j])/  \
                            (lambdav*(DELTA[j] + gama[j]*(1+r_s_VEG[v][j]/r_a_VEG[v][j])))
                else:
                    value = 0.0
                if value >=0.0:
                    PT_PM_VEG_v.append(value)
                else:
                    PT_PM_VEG_v.append(0.0)
                if r_a_VEG[v][j] > 0.0:
                    value = (DELTA[j]*(Rns_VEG[v][j]-Rnl[j]-G_VEG[v][j])+3600*ro_a[j]*Cp*(e0_Ta[j]-e_a[j])/r_a_VEG[v][j])/  \
                            (lambdav*(DELTA[j] + gama[j]))
                else:
                    value = 0.0
                if value >= 0.0:
                    Erf_v.append(value)
                else:
                    Erf_v.append(0.0)
            PT_PM_VEG.append(PT_PM_VEG_v)
            Erf_VEG.append(Erf_v)
        PT_PM_VEG = np.asarray(PT_PM_VEG)
        Erf_VEG = np.asarray(Erf_VEG)
        print "\nPT of VEGETATION computed!"
        del PT_PM_VEG_v, Erf_v, value
    # FIELD
    if IRR <> None:
        PT_PM_FIELD = []
        Erf_FIELD = []  # Evaporation during RF
        if NFIELD>0:
            for f in range(NFIELD):
                PT_PM_FIELD_f = []
                Erf_f = []
                for j in range(len(J)):
                    if r_a_FIELD[f][j] > 0.0:
                        value = (DELTA[j]*(Rns_FIELD[f][j]-Rnl[j]-G_FIELD[f][j])+3600*ro_a[j]*Cp*(e0_Ta[j]-e_a[j])/r_a_FIELD[f][j])/  \
                                (lambdav*(DELTA[j] + gama[j]*(1+r_s_FIELD[f][j]/r_a_FIELD[f][j])))
                    else:
                        value = 0.0
                    if value >=0.0:
                        PT_PM_FIELD_f.append(value)
                    else:
                        PT_PM_FIELD_f.append(0.0)
                    if r_a_FIELD[f][j] > 0.0:
                        value = (DELTA[j]*(Rns_FIELD[f][j]-Rnl[j]-G_FIELD[f][j])+3600*ro_a[j]*Cp*(e0_Ta[j]-e_a[j])/r_a_FIELD[f][j])/  \
                                (lambdav*(DELTA[j] + gama[j]))
                    else:
                        value = 0.0
                    if value >= 0.0:
                        Erf_f.append(value)
                    else:
                        Erf_f.append(0.0)
                PT_PM_FIELD.append(PT_PM_FIELD_f)
                Erf_FIELD.append(Erf_f)
            PT_PM_FIELD = np.asarray(PT_PM_FIELD)
            Erf_FIELD = np.asarray(Erf_FIELD)
            print "\nPT of FIELDS/CROPS computed!"
            del PT_PM_FIELD_f, Erf_f, value
    # for SOIL
    PE_PM_SOIL = []
    if NSOIL > 0:
        for s in range(NSOIL):
            value = (DELTA*(Rns_SOIL[s]-Rnl-G_SOIL[s])+3600*ro_a*Cp*(e0_Ta-e_a)/r_a_SOIL)/  \
                        (lambdav*(DELTA + gama*(1+r_s_SOIL[s]/r_a_SOIL)))
            PE_PM_SOIL.append(np.where(value >= 0.0, np.where(value <= E0, value, E0), 0.0))
        PE_PM_SOIL = np.asarray(PE_PM_SOIL)
        print "\nPE of SOIL computed!"
        del value

    # #### DAILY SUM ##############################################
    if IRR <> None:
        RF_irr = np.zeros([NFIELD, len(RF)], dtype = float)
        for f in range(NFIELD):
            try:
                RF_irr[f] = RF + IRR[f]
            except:
                cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nIrrigation time serie incompatible with rainfall time serie!")
            RF_irr_d =  np.zeros([NFIELD, len(datenum_d)], dtype=float)
            RFint_irr = np.zeros([NFIELD, len(datenum_d)], dtype=float)
            RF_irr_duration =  np.zeros([NFIELD, len(datenum_d)], dtype=float)
    J_day = np.zeros([len(datenum_d)], dtype=float)
    PT_PM_VEG_d = np.zeros([NVEG,len(datenum_d)], dtype=float)
    LAI_veg_d = np.zeros([NVEG,len(datenum_d)], dtype=float)
    if IRR <> None:
        PT_PM_FIELD_d = np.zeros([NFIELD,len(datenum_d)], dtype=float)
    if NSOIL > 0:
        PE_PM_SOIL_d = np.zeros([NSOIL,len(datenum_d)], dtype=float)
    E0_d =  np.zeros([len(datenum_d)], dtype=float)
    RF_veg_d =  np.zeros([len(datenum_d)], dtype=float)
    RFint_veg = np.zeros([len(datenum_d)], dtype=float)
    RF_veg_duration =  np.zeros([len(datenum_d)], dtype=float)
    t_d = 0
    n1 = 0
    n1_d = []
    actual_day = mpl.dates.num2date(datenum[0]).isoformat()[:10]
    J_day[0] = J[0]
    for t in range(len(datenum)):
        if actual_day == mpl.dates.num2date(datenum[t]).isoformat()[:10]:
            n1 = n1 + 1
            for v in range(NVEG):
                PT_PM_VEG_d[v,t_d] = PT_PM_VEG_d[v,t_d] + PT_PM_VEG[v,t]
                LAI_veg_d[v,t_d] = (LAI_veg_d[v,t_d] + LAI_v[v,t])/2.0
            for s in range(NSOIL):
                PE_PM_SOIL_d[s,t_d] = PE_PM_SOIL_d[s,t_d]  + PE_PM_SOIL[s,t]
            E0_d[t_d] = (E0_d[t_d] + E0[t])
            if RF[t]>0:
                RF_veg_d[t_d] = RF_veg_d[t_d] + RF[t]
                RF_veg_duration[t_d] = RF_veg_duration[t_d] + 1.0
            if IRR <> None:
                for f in range(NFIELD):
                    if RF_irr[f][t]>0:
                        RF_irr_d[f][t_d] = RF_irr_d[f][t_d] + RF_irr[f][t]
                        RF_irr_duration[f][t_d] = RF_irr_duration[f][t_d] + 1.0
                    PT_PM_FIELD_d[f,t_d] = PT_PM_FIELD_d[f,t_d] + PT_PM_FIELD[f,t]
        else:
            if RF_veg_duration[t_d]>0:
                RFint_veg[t_d] = RF_veg_d[t_d]/RF_veg_duration[t_d]
            if IRR <> None:
                for f in range(NFIELD):
                    if RF_irr_duration[f][t_d]>0:
                        RFint_irr[f][t_d] = RF_irr_d[f][t_d]/RF_irr_duration[f,t_d]
            n1_d.append(n1+1)
            n1 = 0
            t_d = t_d + 1
            actual_day = mpl.dates.num2date(datenum[t]).isoformat()[:10]
            J_day[t_d] = J[t]
            for v in range(NVEG):
                PT_PM_VEG_d[v,t_d] = PT_PM_VEG[v,t]
                LAI_veg_d[v,t_d] = LAI_v[v,t]
            for s in range(NSOIL):
                PE_PM_SOIL_d[s,t_d] = PE_PM_SOIL[s,t]
            E0_d[t_d] = E0[t]
            if RF[t]>0:
                RF_veg_d[t_d] = RF_veg_d[t_d] + RF[t]
                RF_veg_duration[t_d] = RF_veg_duration[t_d] + 1.0
            if IRR <> None:
                for f in range(NFIELD):
                    RF_irr_d[f][t_d] = RF_irr_d[f][t_d] + RF_irr[f][t]
                    RF_irr_duration[f][t_d] = RF_irr_duration[f][t_d] + 1.0
                    PT_PM_FIELD_d[f,t_d] = PT_PM_FIELD[f][t]
    if n1<>0:
        if RF_veg_duration[t_d]>0:
            RFint_veg[t_d] = RF_veg_d[t_d]/RF_veg_duration[t_d]
        if IRR <> None:
            for f in range(NFIELD):
                RFint_irr[f][t_d] = RF_irr_d[f][t_d]/RF_irr_duration[f][t_d]
        n1_d.append(n1+1)

    #  #####  COMPUTING RF/INTERCEPTION ##############################################
    print "\nComputing INTERCEPTION..."
    """Gash sparse model
        Gash 1979
        Gash, Lloyd and Lachaud 1995
        Valente et al 1997
        Finch 2001
    """
    c = 1   # I is computed for each landcover assuming that evaporation is 1D vertical with no horizontal interaction or advection, see Gash 1995 chapter 2.3 pag 82
    RFtreshold = 0.4 # INDICATE HERE A RF TRESHOLD LIKE INDICATED IN PEREIRA 2009 or GASH 1979

    # select hourly RF average (and corresponding evaporation) above treshold necessary to saturate the canopy, following method of Gash 1979 pag 49
    def RFsat(datenum_, RF_, RFtreshold_):
        RF_sat_ = []
        for t in range(len(datenum_)):
            if RF_[t] > RFtreshold_:
                RF_sat_.append(RF_[t])
        if len(RF_sat_)>0:
            avRF_sat_ = sum(RF_sat_)/len(RF_sat_)
        else:
            print '\nWARNING!!! There is no RF during the input period!'
            avRF_sat_ = 0.0
        return avRF_sat_, len(RF_sat_)

    def Erf_sat(datenum_, RF_, Erf_, avRF_sat_, S_, c_, RFtreshold_):
        Erf_sat_ = []
        for t in range(len(datenum_)):
            if RF_[t] > RFtreshold_:
                Erf_sat_.append(Erf_[t])
        if len(Erf_sat_)>0:
            avErf_sat_ = sum(Erf_sat_)/len(Erf_sat_)
            Pgl_ = -avRF_sat_*(S_/c_)*np.log(1-avErf_sat_/avRF_sat_)/avErf_sat_
        else:
            print '\nWARNING!!! There is no RF during the input period!'
            avErf_sat_ = 0.0
            Pgl_ = 0.0
        return Pgl_, avErf_sat_, len(Erf_sat_)

    # INTERCEPTION VEG
    print '\n-----------\nInterception by VEGETATION'
    I_veg_d = np.zeros([NVEG,len(datenum_d)], dtype=float)
    RFe_veg_d = np.zeros([NVEG,len(datenum_d)], dtype=float)
    avRFsat_veg, len_RFsat_veg = RFsat(datenum_ = datenum, RF_ = RF, RFtreshold_ = RFtreshold)
    print 'RF average = ' + '%.4f' %avRFsat_veg + ' mm/h'
    print 'Values computed for ' + str(len_RFsat_veg) + ' events of hourly RF higher than ' + str(RFtreshold) + ' mm'
    for v in range(NVEG):
        Pgl, avErf_sat, len_Erf_sat = Erf_sat(datenum_ = datenum, RF_ = RF, Erf_ = Erf_VEG[v], avRF_sat_ = avRFsat_veg, S_ = S_w_v[v], c_ = c, RFtreshold_ = RFtreshold)
        print '\nVegetation type ' + str(v)  + ' (' + VegType[v] + ')'
        print 'Evaporation average during RF = ' + '%.4f' % (avErf_sat) + ' mm/h (' + str(len_Erf_sat)   + ' values)'
        print 'RF treshold to saturate the canopy (Pgl) = ' + '%.2f' %Pgl
        if avRFsat_veg > 0.0:
            I_veg_d[v] = np.where(RF_veg_d > Pgl, c*Pgl + (c*avErf_sat/avRFsat_veg)*(RF_veg_d-Pgl), c*RF_veg_d)  #- c*S_w_v[v]
            RFe_veg_d[v] = np.where(RF_veg_d > Pgl, RF_veg_d-I_veg_d[v], 0.0)

    # INTERCEPTION IRR
    if IRR <> None:
        print '\n-----------\nInterception by CROPS in FIELDS'
        I_irr_d = np.zeros([NFIELD,len(datenum_d)], dtype=float)
        RFe_irr_d = np.zeros([NFIELD,len(datenum_d)], dtype=float)
        for f in range(NFIELD):
            print '\nFIELD #%i' % (f+1)
            avRFsat_irr, len_RFsat_irr = RFsat(datenum_ = datenum, RF_ = RF_irr[f], RFtreshold_ = RFtreshold)
            print 'RF average = ' + '%.4f' %avRFsat_irr + ' mm/h'
            print 'Values computed for ' + str(len_RFsat_irr) + ' events of hourly RF higher than ' + str(RFtreshold) + ' mm'
            if S_w_f[f].sum() > 0.0:
                Pgl, avErf_sat, len_Erf_sat = Erf_sat(datenum_ = datenum, RF_ = RF, Erf_ = Erf_FIELD[f], avRF_sat_ = avRFsat_irr, S_ = S_w_f[f], c_ = c, RFtreshold_ = RFtreshold)
                print 'Evaporation average during RF = ' + '%.4f' % (avErf_sat) + ' mm/h (' + str(len_Erf_sat)   + ' values)'
                print 'Averaged RF treshold to saturate the canopy (Pgl) = ' + '%.2f' % (sum(Pgl)/sum(RF_irr[f]>0.0))
                if avRFsat_irr > 0.0:
                    I_irr_d[f] = np.where(RF_irr_d[f] > Pgl, c*Pgl + (c*avErf_sat/avRFsat_irr)*(RF_irr_d[f]-Pgl), c*RF_irr_d[f])  #- c*S_w_f[f]
                    RFe_irr_d[f] = np.where(RF_irr_d[f] > Pgl, RF_irr_d[f]-I_irr_d[f], 0.0)
            else:
                print 'No crop'
                RFe_irr_d[f] = RF_irr_d[f]
    print '-----------'

    # RETURN ARRAY VALUES TO MAIN PROGRAM
    if IRR == None:
        PT_PM_FIELD = PT_PM_FIELD_d= Erf_FIELD = RF_irr_d = RFint_irr = RF_irr_duration = RFe_irr_d = I_irr_d = []
    return J, J_day, outputVAR, PT_PM_VEG, Erf_VEG, PE_PM_SOIL, E0, PT_PM_VEG_d, PE_PM_SOIL_d, E0_d, RF_veg_d, RFint_veg, RF_veg_duration, n1_d, RFe_veg_d, I_veg_d, LAI_veg_d, PT_PM_FIELD, PT_PM_FIELD_d, Erf_FIELD, RF_irr_d, RFint_irr, RF_irr_duration, RFe_irr_d, I_irr_d

# EOF #