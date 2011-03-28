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
__version__ = "0.0"
__date__ = "November 2010"

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
    print 'from hourly meteorological data for different land covers'
    print '(vegetation, open water, bare soil)'
    print 'using Penman-Monteith equation.'
    print 'Type PET.PM() for help'
    print 'Author: ',__author__
    print 'Version: ',__version__
    print 'Date: ',__date__
    return

import numpy as np
import pylab
from sys import path

'''
    ================================================================
    PM function
    ================================================================
'''

def process(datenum = np.array([]), datenum_d = np.array([])
              ,RF = np.array([]), Ta = np.array([]), RHa = np.array([]) \
              ,Pa = np.array([]), u_z_m = np.array([]), Rs = np.array([]) \
              ,phi = 0.0, Lm = 0.0 , Z = 0.0, Lz = 0.0, FC = 0.0, z_m=2.0, z_h=2.0 \
              ,NVEG = 1, VegType =np.array([]), h = np.array([]) \
              ,S = np.array([]), C_leaf_star = np.array([]) \
              ,LAI_d = np.array([]), LAI_w = np.array([])  \
              ,f_s_d = np.array([]), f_s_w = np.array([])  \
              ,alfa_vd = np.array([]), alfa_vw = np.array([]), J_vd = 91, J_vw= 305, TRANS_vdw = 20\
              ,NSOIL = 1, SoilType = [], Sy = np.array([])\
              ,alfa_sd = np.array([]), alfa_sw = np.array([]), J_sd = 166, J_sw = 274, TRANS_sdw = 20
              ,alfa_w = 0.06):
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
        f_s_d: shelter factor dry season []
        f_s_w: shelter factor wet season []
        alfa_vd: vegetation albedo in dry season
        alfa_vw: vegetation albedo in wet season
        J_vd: starting julian day of the dry season [int 1-365]
        J_vw: starting julian day of the wet season [int 1-365]

# TO BE UPDATED, see startMARMITESsurface.py
        SOIL PARAMETERS
        NSOIL: number of soil types
        to repeat NSOIL times in a same line
        Sy: surface (1cm) soil specific yield [m3.m-3]
        alfa_sd: soil albedo in dry season
        alfa_sw: soil albedo in wet season
        jd_sd: starting julian day of the dry season [int 1-365]
        jd_sw: starting julian day of the wet season [int 1-365]

        WATER PARAMETERS
        alfa_w: water albedo

    Output:
    ---------
        PET_PM_VEG    : hourly PET [mm/day] for NVEG type(s) of vegetation covers
        PE_PM_SOIL    : hourly PE [mm/day] for NSOIL type(s) of soil types
        E0            : hourly PE [mm/day] for open water
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

    # compute J - julian day
    YYYY = []
    MM = []
    DD = []
    HH = []
    MN = []
    J = []
    time = []
    for t in range(len(datenum)):
        YYYY.append(float('%4d'%pylab.num2date(datenum[t]).year))
        MM.append(float('%02d'%pylab.num2date(datenum[t]).month))
        DD.append(float('%02d'%pylab.num2date(datenum[t]).day))
        HH.append('%02d'%pylab.num2date(datenum[t]).hour)
        MN.append('%02d'%pylab.num2date(datenum[t]).minute)
        J.append(DD[t] - 32 + int(275*MM[t]/9) + 2 * int(3/(MM[t] + 1)) + int(MM[t]/100-np.mod(YYYY[t],4)/4+0.975) )
        time.append(HH[t]+':' + MN[t])

    # compute DELTA - SLOPE OF SATURATION VAPOUR PRESSURE CURVE
    # [kPa.ªC-1]
    # FAO56 pag 37 Eq13
    DELTA = []
    for j in range(len(J)):
        DELTA.append(4098*(0.6108*np.exp(17.27*Ta[j]/(Ta[j]+237.3)))/pow(Ta[j]+237.3,2))
    print "\nDELTA computed!"

    # compute dr - inverse distance to the sun
    # FAO56 pag47 Eq23
    dr = []
    for j in range(len(J)):
        dr.append(1+0.033*pylab.cos(2*np.pi*J[j]/365))

    # compute delta - solar declination
    # [rad]
    # FAO56 pag47 Eq24
    delta = []
    for j in range(len(J)):
        delta.append(0.409*pylab.sin(2*np.pi*J[j]/365-1.39))

    # compute Sc - seasonnal correction for solar time
    # [hour]
    # FAO56 pag47 Eq32
    Sc = []
    for j in range(len(J)):
        b = 2*np.pi*(J[j]-81)/364    # Eq 34
        Sc.append(0.1645*pylab.sin(2*b) - 0.1255*pylab.cos(b) - 0.025*pylab.sin(b))

    # compute w - solar time angle at the midpoint of the period (time)
    # [rad]
    # FAO56 pag48 Eq31
    w = []
    for j in range(len(J)):
        t = -0.5 + FC + float(time[j].split(':')[0]) + \
                  float(time[j].split(':')[1])/60.0
        w.append((np.pi/12)*((t+0.06667*(Lz-Lm)+Sc[j])-12))

    # compute w1 - solar time angle at the beginning of the period (time)
    # [rad]
    # FAO56 pag47 Eq29
    w1 = []
    tl = 1  # hourly data
    for j in range(len(J)):
        w1.append(w[j] - np.pi*tl/24)

    # compute w2 - solar time angle at the end of the period (time + 1h)
    # [rad]
    # FAO56 pag47 Eq30
    w2 = []
    for j in range(len(J)):
        w2.append(w[j] + np.pi*tl/24)

    # compute ws - sunset hour angle
    # [rad]
    # FAO56 pag47 Eq25
    ws = []
    for j in range(len(J)):
        ws.append(pylab.arccos(-pylab.tan(phi*np.pi/180)*pylab.tan(delta[j])))

    # compute Ra - extraterrestrial radiation
    # [MJ.m-2.hour-1]
    # FAO56 pag47 Eq28
    Ra = []
    Ra_Watts = []
    for j in range(len(J)):
        if (w1[j] < -ws[j] or w2[j] > -ws[j]):
            if w1[j] < -ws[j] : w1[j] = -ws[j]
            if w2[j] < -ws[j] : w2[j] = -ws[j]
        if (w1[j] < ws[j] or w2[j] > ws[j]):
            if w1[j] >  ws[j] : w1[j] =  ws[j]
            if w2[j] >  ws[j] : w2[j] =  ws[j]
        Ra.append((12*60/np.pi)*Gsc*dr[j]* \
              ((w2[j]-w1[j])*pylab.sin(phi*np.pi/180)*pylab.sin(delta[j]) + \
               pylab.cos(phi*np.pi/180)*pylab.cos(delta[j])*(pylab.sin(w2[j])-pylab.sin(w1[j]))))
        Ra_Watts.append(Ra[j]*24/0.08864)

    # compute Rs0 - clear-sky solar (shortwave) radiation
    # [MJ.m-2.hour-1]
    # FAO56 pag51 Eq37
    Rs0 = []
    Rs0_Watts = []
    for j in range(len(J)):
        Rs0.append((0.75+2E-5*Z)*Ra[j])
        Rs0_Watts.append(Rs0[j]*24/0.08864)

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
        print "\n" + '%.1f' %(100*badvalues/len(J)) + " % (" + str(badvalues) + '/' + str(len(J)) +') shortwave radiation measurements that were higher than the clear-sky one were corrected.'
        if badvalues/len(J)>0.15:
            print '\nCalibrate your sensors!!!'

    # compute Rns - NET SHORTWAVE RADIATION
    # [MJ.m-2.hour-1]
    # FAO56 pag51 Eq37
    # for each type of vegetation and soil (albedo dependent)
    Rns_VEG = []
    if NVEG > 0:
        for v in range(NVEG):
            Rns_v = []
            alfa = []
            for j in range(len(J)):
                if J[j] < J_vd[v] or J[j] > J_vw[v]: # wet period
                    alfa.append(alfa_vw[v])
                else:
                    if J[j] < J_vd[v] + TRANS_vdw[v]: # transition wet to dry
                        alfa.append(alfa[j-1]-(alfa_vw[v]-alfa_vd[v])/(24*TRANS_vdw[v]+24))
                    elif J[j] > J_vw[v] - TRANS_vdw[v]: # transition dry to wet
                        alfa.append(alfa[j-1]+(alfa_vw[v]-alfa_vd[v])/(24*TRANS_vdw[v]+24))
                    else:                              # dry period
                        alfa.append(alfa_vd[v])
                Rns_v.append((1-alfa[j])*Rs_corr[j])
            Rns_VEG.append(Rns_v)
            del Rns_v, alfa
    Rns_SOIL = []
    if NSOIL > 0:
        for s in range(NSOIL):
            Rns_s = []
            alfa = []
            for j in range(len(J)):
                if J[j] < J_sd[s] or J[j] > J_sw[s]:
                    alfa.append(alfa_sw[s])
                else:
                    if J[j] < J_sd[s] + TRANS_sdw[s]:
                        alfa.append(alfa[j-1]-(alfa_sw[s]-alfa_sd[s])/(24*TRANS_sdw[s]+24))
                    elif J[j] > J_sw[s] - TRANS_sdw[s]:
                        alfa.append(alfa[j-1]+(alfa_sw[s]-alfa_sd[s])/(24*TRANS_sdw[s]+24))
                    else:
                        alfa.append(alfa_sd[s])
                Rns_s.append((1-alfa[j])*Rs_corr[j])
            Rns_SOIL.append(Rns_s)
            del Rns_s, alfa
    Rns_WATER = []
    for j in range(len(J)):
        Rns_WATER.append((1-alfa_w)*Rs_corr[j])
    print "\nRns computed!"

    # compute e0_Ta - saturation vapour pressure at actual air temperature
    # [kPa]
    # FAO56 pag36 Eq11
    e0_Ta = []
    for i in range(len(J)):
        e0_Ta.append(0.6108*np.exp(17.27*Ta[i]/(Ta[i]+237.3)))

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
        print "\n" + '%.1f' % (100*badlowvalues/len(J)) + " % (" + str(badlowvalues) + "/" + str(len(J)) + ") relative humidity measurements > 100% were fixed to 100%."
        if badvalues/len(J)>0.15:
            print '\nCalibrate your sensors!!!'

    # compute e_a - ACTUAL VAPOUR PRESSURE
    # [kPa]
    # FAO56 pag74 Eq54
    e_a = []
    for i in range(len(J)):
        e_a.append(e0_Ta[i]*RHa[i]/100)
    print "\ne_a computed!"

    # compute Rnl - NET LONGWAVE RADIATION
    # [MJ.m-2.hour-1]
    # FAO56 pag51 Eq37 and pag74 for hourly computing
    Rnl = []
    r = []
    i = 0
    for j in range(len(J)):
        if (ws[j] - 0.52) <= w[j] <= (ws[j] - 0.10):  # FAO56: (ws[j] - 0.79) <= w[j] <= (ws[j] - 0.52)
            i = 1
            if Rs0[j] > 0:
                r_sunset = Rs_corr[j]/Rs0[j]
            else:
                r_sunset = 0.75  #see FAO56 pag75
        if ((ws[j] - 0.10) < w[j] or w[j] <= (-ws[j]+ 0.10)):
            if i>0:
                r.append(r_sunset)
            else:
                r.append(0.75) #see FAO56 pag75
        else:
            r.append(Rs_corr[j]/Rs0[j])
        Rnl.append(sigma*pow(Ta[j] + 273.16,4)*(0.34-0.14*pylab.sqrt(e_a[j]))*(1.35*r[j]-0.35))
##        if Rnl[j]<0:
##            r=0.8
##            Rnl[j] = sigma*pow(Ta[j] + 273.16,4)*(0.34-0.14*pylab.sqrt(e_a[j]))*(1.35*r-0.35)
    print "\nRnl computed!"

    # compute G - SOIL HEAT FLUX
    # [MJ.m-2.hour-1]
    # FAO56 pag55 Eq45 and 46
    G_VEG = []
    if NVEG>0:
        for v in range(NVEG):
            G_v = []
            for j in range(len(J)):
                Rn = Rns_VEG[v][j] - Rnl[j]
                if w[j]<-ws[j] or w[j]>ws[j]:
                    G_v.append(0.5*Rn)
                else:
                    G_v.append(0.1*Rn)
            G_VEG.append(G_v)
    G_SOIL = []
    if NSOIL>0:
        for s in range(NSOIL):
            G_s = []
            for j in range(len(J)):
                Rn = Rns_SOIL[s][j] - Rnl[j]
                if w[j]<-ws[j] or w[j]>ws[j]:
                    G_s.append(0.5*Rn)
                else:
                    G_s.append(0.1*Rn)
            G_SOIL.append(G_s)
    G_WATER = []
    for j in range(len(J)):
        Rn = Rns_WATER[j] - Rnl[j]
        if w[j]<-ws[j] or w[j]>ws[j]:
            G_WATER.append(0.5*Rn)
        else:
            G_WATER.append(0.1*Rn)
    print "\nG computed!"
    G_VEG_Watts = []
    for j in range(len(J)):
        G_VEG_Watts.append(G_VEG[0][j]*24/0.08864)

    # ro_a - MEAN AIR DENSITY AT CTE PRESSURE
    # [kg.m-3]
    # FAO56 pag26 box6
    ro_a = []
    for j in range(len(J)):
        ro_a.append(Pa[j]/(R*1.01*(Ta[j]+273.16)))
    print "\nro_a computed!"

    # gama - PSYCHROMETRIC CONSTANT
    # [kPa.ªC-1]
    # FAO56 pag31 eq8
    gama = []
    for j in range(len(J)):
        gama.append(Cp*Pa[j]/(eps*lambdav))
    print "\ngama computed!"

    # r_s - SURFACE RESISTANCE
    # [s.m-1]
    # VEG: Dingman pag 208 (canopy conductance) (equivalent to FAO56 pag21 Eq5)
    r_s_VEG = []
    f_k = []
    DELTArho_v = []
    f_rho = []
    f_T = []
    for v in range(NVEG):
        r_s_VEG.append([])
    if NVEG>0:
        for v in range(NVEG):
            f_s_tmp = []
            LAI_tmp =[]
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
                    C_leaf = C_leaf_star[v]
                else:
                    C_leaf = C_leaf_star[v] * f_k[j] * f_rho[j] * f_T[j]
                if J[j] < J_vd[v] or J[j] > J_vw[v]:    #wet period
                    f_s_tmp.append(f_s_w[v])
                    LAI_tmp.append(LAI_w[v])
                else:
                    if J[j] < J_vd[v] + TRANS_vdw[v]:# transition wet to dry
                        f_s_tmp.append(f_s_tmp[j-1]-(f_s_w[v]-f_s_d[v])/(24*TRANS_vdw[v]+24))
                        LAI_tmp.append(LAI_tmp[j-1]-(LAI_w[v]-LAI_d[v])/(24*TRANS_vdw[v]+24))
                    elif J[j] > J_vw[v] - TRANS_vdw[v]:  # transition dry to wet
                        f_s_tmp.append(f_s_tmp[j-1]+(f_s_w[v]-f_s_d[v])/(24*TRANS_vdw[v]+24))
                        LAI_tmp.append(LAI_tmp[j-1]+(LAI_w[v]-LAI_d[v])/(24*TRANS_vdw[v]+24))
                    else:                               #dry period
                        f_s_tmp.append(f_s_d[v])
                        LAI_tmp.append(LAI_d[v])
                f_temp = f_s_tmp[j]*LAI_tmp[j]*C_leaf
                if f_temp == 0.0:
                    r_s_VEG[v].append(1.0E6)
                else:
                    r_s_VEG[v].append(1.0/f_temp)
     # SOIL: van de Griend and Owe, 1994
    r_s_SOIL = []
    if NSOIL>0:
        for s in range(NSOIL):
             r_s_SOIL.append(10*np.exp(0.3563*Sy[s]))
        print "\nr_s computed!"


    # correction windspeed measurement and scaling at h+2m
    # [m.s-1]
    # FAO56 pag56 eq47
    u_2 = []
    u_hplus2 = []
    for v in range(NVEG):
        u_hplus2.append([])
    for j in range(len(J)):
        if u_z_m[j]<0.0:
            u_z_m[j] = 0.0
        if z_m <> 2.0:
            u_2.append(u_z_m[j] * 4.87 / (pylab.log(67.8*z_m-5.42)))
        else:
            u_2.append(u_z_m[j])
        for v in range(NVEG):
            u_hplus2[v].append(u_2[j] * (pylab.log(67.8*(h[v]+2.0)-5.42)) / 4.87)

    # r_a - AERODYNAMIC RESISTANCE
    # [s.m-1]
    r_a_VEG = []
    if NVEG>0:
        for v in range(NVEG):
            r_a_j=[]
            for j in range(len(J)):
                if u_z_m[j] == 0.0:
                    r_a_j.append(1.0E6)
                else:
                    if v == 0:      # FAO56 pag20 eq4- (d - zero displacement plane, z_0m - roughness length momentum transfer, z_0h - roughness length heat and vapour transfer, [m], FAO56 pag21 BOX4
                        r_a_j.append(pylab.log((2-(2*h[v]/3))/(0.123*h[v]))*pylab.log((2-(2*h[v]/3))/(0.0123*h[v]))/(pow(k,2)*u_2[j]))
                    else:           # DINGMAN pag 296
                        r_a_j.append(pow(pylab.log((h[v]+2-(0.7*h[v]))/(0.1*h[v])),2)/(pow(k,2)*u_hplus2[v][j]))
            r_a_VEG.append(r_a_j)
##    if NVEG>0:
##        for v in range(NVEG):
##            if T[v] == 1:
##                # AWSET manual pag24, r_a=94/u_10 for tree (Thompson et al, 1981)
##                # From Pereira, Gash, David, Valente - Evaporation of intercepted rainfall from isolated evergreen oak trees:
##                #                                 Do the crowns behave as wet bulbs?
##                #                                 Agricultural and forest meteorology 149 - 2009 (667-679)
##                # 1/(0.16*pow(u_TreeCrown,0.441))
##                r_a_j_tree=[]
##                for j in range(len(J)):
###                    if u_hplus2[v][j] == 0.0:
##                    if u_TreeCrown[v][j] == 0.0:
##                        r_a_j_tree.append(1.0E6)
##                    else:
###                        r_a_j_tree.append(1/(0.16*pow(u_TreeCrown[v][j],0.441)))
##                        r_a_j_tree.append(94/u_10[j])
##                r_a_VEG.append(r_a_j_tree)
##            else:
##                if z_m < h[v]:
##                    print "Wind speed and relative humidity measurements have to be done at a heigth higher than the crop heigth!   \
##                         \nNot possible to compute aerodynamic parameters for vegetation type #" + str(NVEG[v]) + \
##                         ", see FAO56 pag 20 Eq4."
##                    r_a_VEG.append(np.zeros([len(J)]))
##                else:
##                    r_a_j=[]
##                    for j in range(len(J)):
##                        if u_z_m[j] == 0.0:
##                            r_a_j.append(1.0E6)
##                        else:
##                            r_a_j.append(pylab.log((z_m-d[v])/z_0m[v])*pylab.log((z_h-d[v])/z_0h[v])/(pow(k,2)*u_z_m[j]))
##                    r_a_VEG.append(r_a_j)
    # r_a for SOIL
    # Liu www.hydrol-earth-syst-sci.net/11/769/2007/
    r_a_SOIL = []
    # only function of ws, it is assumed that roughness are the same for any type of soil
    if NSOIL > 0:
        for j in range(len(J)):
            if u_z_m[j] == 0:
                r_a_SOIL.append(1E6)
            else:
                r_a_SOIL.append(pylab.log((z_m)/0.0058)*pylab.log(z_h/0.0058)/(pow(k,2)*u_z_m[j]))
    ##    # BavelHillel_PotActualEvaporationBareSoilSurface_1976
    ##    # SOIL ra = ([ln(2.0/Z0)]^2)/(0.16*U2) with Z0 = 0.01
    ##    r_a_SOIL = []
    ##    Z0 = 0.01
    ##    for j in range(NSOIL):
    ##        r_a_i=[]
    ##        for i in range(len(J)):
    ##            r_a_i.append(pow(pylab.log(2.0/Z0),2)/(0.16*u_2[i]))
    ##        r_a_SOIL.append(r_a_i)
    print "\nr_a computed!"


    # equation variables for outputs
    outputVAR = [DELTA,gama, Rs, Rs_corr, Rs0, Rnl, r]

    # PE(T) - Penman-Montheith
    # mm.hour-1
    # FAO56 pag19 eq3
    # for VEG
    PET_PM_VEG = []
    Erf = []  # Evaporation during RF
    if NVEG>0:
        for v in range(NVEG):
            PET_PM_VEG_v = []
            Erf_v = []
            for j in range(len(J)):
                PET_PM_VEG_v.append(  (DELTA[j]*(Rns_VEG[v][j]-Rnl[j]-G_VEG[v][j])+3600*ro_a[j]*Cp*(e0_Ta[j]-e_a[j])/r_a_VEG[v][j])/  \
                            (lambdav*(DELTA[j] + gama[j]*(1+r_s_VEG[v][j]/r_a_VEG[v][j])))                  \
                            )
                Erf_v.append(  (DELTA[j]*(Rns_VEG[v][j]-Rnl[j]-G_VEG[v][j])+3600*ro_a[j]*Cp*(e0_Ta[j]-e_a[j])/r_a_VEG[v][j])/  \
                            (lambdav*(DELTA[j] + gama[j]))                 \
                            )
            PET_PM_VEG.append(PET_PM_VEG_v)
            Erf.append(Erf_v)
        print "\nPET for VEGETATION computed!"
    # for SOIL
    PE_PM_SOIL = []
    if NSOIL > 0:
        for s in range(NSOIL):
            PE_PM_SOIL_s = []
            for j in range(len(J)):
                PE_PM_SOIL_s.append(  (DELTA[j]*(Rns_SOIL[s][j]-Rnl[j]-G_SOIL[s][j])+3600*ro_a[j]*Cp*(e0_Ta[j]-e_a[j])/r_a_SOIL[j])/  \
                            (lambdav*(DELTA[j] + gama[j]*(1+r_s_SOIL[s]/r_a_SOIL[j])))                  \
                            )
            PE_PM_SOIL.append(PE_PM_SOIL_s)
        print "\nPE for SOIL computed!"
    # E0 - open water
    # computed by Penman equation in Gieske 2003 pag 94 eq63 and 64
    E0 = []
    for j in range(len(J)):
        E0.append(  (DELTA[j]*(Rns_WATER[j]-Rnl[j]-G_WATER[j])  +         \
                     gama[j]*0.26*(1+0.54*u_2[j])*(e0_Ta[j]-e_a[j])  )/   \
                     (lambdav*(DELTA[j] + gama[j])))
    print "\nE0 for open water computed!"


    # #### DAILY SUM ##############################################
    J_d = []
    PET_PM_VEG_d = np.zeros([NVEG,len(datenum_d)], dtype=float)
    PE_PM_SOIL_d = np.zeros([NSOIL,len(datenum_d)], dtype=float)
    E0_d =  np.zeros([len(datenum_d)], dtype=float)
    RF_d =  np.zeros([len(datenum_d)], dtype=float)
    RFint = np.zeros([len(datenum_d)], dtype=float)
    RF_duration =  np.zeros([len(datenum_d)], dtype=float)
    t_d = 0
    n1 = 0
    n1_d = []
    actual_day = pylab.num2date(datenum[0]).isoformat()[:10]
    J_d.append(J[0])
    for t in range(len(datenum)):
        if actual_day == pylab.num2date(datenum[t]).isoformat()[:10]:
            n1 = n1 + 1
            for v in range(NVEG):
                PET_PM_VEG_d[v][t_d] = PET_PM_VEG_d[v][t_d] + PET_PM_VEG[v][t]
            for s in range(NSOIL):
                PE_PM_SOIL_d[s][t_d] = PE_PM_SOIL_d[s][t_d]  + PE_PM_SOIL[s][t]
            E0_d[t_d] = (E0_d[t_d] + E0[t])
            if RF[t]>0:
                RF_d[t_d] = RF_d[t_d] + RF[t]
                RF_duration[t_d] = RF_duration[t_d] + 1.0
        else:
            if RF_duration[t_d]>0:
                RFint[t_d] = RF_d[t_d]/RF_duration[t_d]
            n1_d.append(n1+1)
            n1 = 0
            t_d = t_d + 1
            actual_day = pylab.num2date(datenum[t]).isoformat()[:10]
            J_d.append(J[t])
            for v in range(NVEG):
                PET_PM_VEG_d[v][t_d] = PET_PM_VEG[v][t]
            for s in range(NSOIL):
                PE_PM_SOIL_d[s][t_d] = PE_PM_SOIL[s][t]
            E0_d[t_d] = E0[t]
            if RF[t]>0:
                RF_d[t_d] = RF_d[t_d] + RF[t]
                RF_duration[t_d] = RF_duration[t_d] + 1.0
    if n1<>0:
        if RF_duration[t_d]>0:
            RFint[t_d] = RF_d[t_d]/RF_duration[t_d]
        n1_d.append(n1+1)


    #  #####  COMPUTING RF/INTERCEPTION ##############################################
    print '\n##############'
    print "Computing INTERCEPTION..."

    """Gash sparse model
        Gash 1979
        Gash, Lloyd and Lachaud 1995
        Valente et al 1997
        Finch 2001
    """

    c=1   # I is computed for each landcover assuming that evaporation is 1D vertical with no horizontal interaction or advection, see Gash 1995 chapter 2.3 pag 82
    I_d = []
    RFe_d = []
    Pgl = []

    # select hourly RF average (and corresponding evaporation) above treshold necessary to saturate the canopy, following method of Gash 1979 pag 49
    RF_sat =  []
    Erf_sat = []
    for v in range(NVEG):
        Erf_sat.append([])
        I_d.append([])
        RFe_d.append([])
        Pgl.append([])
    datenum_RF_sat = []
    RFtreshold = 0.4 # INDICATE HERE A RF TRESHOLD LIKE INDICATED IN PEREIRA 2009 or GASH 1979
    for t in range(len(datenum)):
        if RF[t]>RFtreshold:
            datenum_RF_sat.append(datenum)
            RF_sat.append(RF[t])
            for v in range(NVEG):
                Erf_sat[v].append(Erf[v][t])
    if len(RF_sat)>0:
        avRF_sat = sum(RF_sat)/len(RF_sat)
        print '\nRF average = ' + '%.4f' %avRF_sat + ' mm/h'
        avErf_sat=[]
        for v in range(NVEG):
            avErf_sat.append(sum(Erf_sat[v])/len(Erf_sat[v]))
            Pgl[v] = -avRF_sat*(S[v]/c)*pylab.log(1-avErf_sat[v]/avRF_sat)/avErf_sat[v]
            print '\nVegetation type ' + str(v)  + ' (' + VegType[v] + ')'
            print 'Erf average during RF = ' + '%.4f' %avErf_sat[v] + ' mm/h (' + str(len(Erf_sat[v]))   + ' values)'
            print 'RF treshold to saturate the canopy (Pgl) = ' + '%.2f' %Pgl[v]
        print '\nValues computed for ' + str(len(RF_sat)) + ' events of hourly RF higher than ' + str(RFtreshold) + ' mm'
    else:
        print '\nWARNING!!! There is no RF during the input period!'
        avErf_sat = 0

    # INTERCEPTION
    for t in range(len(RF_d)):
        for v in range(NVEG):
            if RF_d[t]>Pgl[v]:
                I_d[v].append(c*Pgl[v] + (c*avErf_sat[v]/avRF_sat)*(RF_d[t]-Pgl[v]))  #- c*S[v]
                RFe_d[v].append(RF_d[t]-I_d[v][t])
            else:
                RFe_d[v].append(0)
                I_d[v].append(c*RF_d[t])

    # RETURN ARRAY VALUES TO MAIN PROGRAM
    return J, J_d, outputVAR, PET_PM_VEG, Erf, PE_PM_SOIL, E0, PET_PM_VEG_d, PE_PM_SOIL_d, E0_d, RF_d, RFint, RF_duration, n1_d, RFe_d, I_d


####################################################
# EOF #