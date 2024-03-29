# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:        startPET_PM_FAO56.py
# Purpose:
#
# Author:      Alain Francés
#
# Created:     25-11-2010
# Copyright:   (c) alf 2010
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__author__ = "Alain P. Francés <frances08512@itc.nl>"
__version__ = "0.3"
__date__ = "2012"

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import PET_P_INTER, plotPET, plotP

'''
    Reads input files for PET_PM_FAO56 and call this function

    Input:
    -------------------------------
    Two input ASCII files (space separated) are necessary:
        1 - METEO TIME SERIES (hourly data)
        2 - METEO/VEGETATION/SOIL/OPENWATER PARAMETERS

    1 - METEO TIME SERIES file format (call it inputMETEOTS.txt)
    This data file is organised in columns separated by space with the following
    format/units:
    column 1:     date: date [yyyy-mm-dd]
    column 2:     time: time [hh-mm]
    column 3:     P: rainfall [mm]
    column 4:     Ta: air temperature measured at heigth z_h [ºC]
    column 5:     RHa: relative air humidity [%] measured at heigth z_h [m]
    column 6:     Pa: air pressure [kPa]
    column 7:     u_z_m: windspeed measured at heigth z_m [m.s-1]
    column 8:     Rs: incoming solar (=shortwave) radiation [MJ.m-2.hour-1]
    Example 19 from FAO56 pag 75 (P added by myself)
    - BOF - don't copy this line in a txt file
    date       time  P Ta RHa Pa      u_z_m   Rs  # the first line is not read but is compulsory
    2010-10-01 02:00 0.0 28 90  101.325 1.9     0.0
    2010-10-01 14:00 0.0 38 52  101.325 3.3     2.450
    - EOF - don't copy this line in a txt file

    2 - METEO/VEGETATION/SOIL/OPENWATER PARAMETERS file format (call it inputPARAM.txt)
    - BOF - don't copy this line in a txt file
    #
    # - the first character in the first line is the delimitation character (DC)
    # - a line starting with DC will not be readed by the python script (it is a user comments)
    # - USE A DC at the end of folder or file names with space (i.e. "c:\My MARMITES\input.ini#")
    # - you can use as many comment and blank line as you want

    # METEO PARAMETERS
    # NMETEO: number of meteo zones
    1
    #for each meteo zone, enter the following parameters in one single line sparated by space
    # enter the following parameters in one single line sparated by space
    # it is assumed that the data acquistion is done at winter time
    # 1 - phi: latitude of the meteo station [degres, >0 hemisphere N]
    # 2 - Lm:  longitude of the meteo station [degres west from Greenwich]
    # 3 - Z: altitude of the station above sea level [m]
    # 4 - Lz:  longitude of the center of the local time zone [degres west of Greenwich]
    # 5 - FC:  for sites not located geographically in the same fuse as the one of the local time zone, indicate the time shift (generally +/-1) [h]
    # 6 - DTS: data time shift, for data NOT acquired at standart clock time [h]
    # 7 - z_m: heigth of u_z_m measurement [m]
    # 8 - z_h: heigth of humidity measurement [m]
    # Example 19 from FAO56 pag 75
    16.22 16.25 8.0 15.0 0.0 0.0 2.0 2.0

    #VEGETATION PARAMETERS
    #NVEG: number of vegetation types (NO IRRIGATION)
    3
    # to repeat NVEG times in a same line
    # 0 - VegName: name of the vegetation type [string, 1 word, no space allowed]
    # 1 - h_vd: heigth of plant dry season[m]
    # 2 - h_vw: heigth of plant wet season[m]
    # 3 - S_w_v: canopy capacity [mm]
    # 4 - C_leaf_star_v: maximum leaf conductance [m.s-1]
    # 5 - LAI_vd: leaf area index dry season [m2.m-2]
    # 6 - LAI_vw: leaf area index wet season [m2.m-2]
    # 7 - f_s_vd: shelter factor dry season []
    # 8 - f_s_vw: shelter factor wet season []
    # 9 - alfa_vd: vegetation albedo in dry season []
    # 10 - alfa_vw: vegetation albedo in wet season []
    # 11 - J_vd: starting julian day of the dry season [int 1-365]
    # 12 - J_vw: starting julian day of the wet season [int 1-365]
    # 13 - TRANS_vdw: transition period between dry and wet season [days]
    # 14 - TRANS_vwd: transition period between wet and dry season [days]
    # 15 - Zr: maximum root depth [m]
    # 16 - kTg_min: transpiration sourcing factor min [], 1>=k_T>0
    # 17 - kTg_max: transpiration sourcing factor max [], 1<k_T<=0
    # 18 - kT_f: transpiration sourcing factor f [], phi>f>wp
    # 19 - kT_s: transpiration sourcing factor slope [], 1>s>0
    # BY DEFAULT THE PROGRAM WILL COMPUTE ETref USING GRASS FAO56 PARAMETERS
    # grassFAO56 0.12 0.12 0.1 0.01 2.88 2.88 0.5 0.5 0.23 0.23 150 270 20 20 0.25 0.01 0.011 2.0
    # VEG1 (grass muelledes)
    grassMU 0.01 0.5 0.15 0.010 0.00 2.88 0.5 0.55 0.50 0.20 150 300 20 30 0.40 0.01 0.011 2.0
    # VEG2 (Q.ilex Sardon)
    Qilex   6.0 6.0 1.16 0.004 6.0 6.0 0.75 0.70 0.25 0.10 150 290 20 50 15.0 0.05 0.85 0.5
    # VEG3 (Q.pyr sardon)
    Qpyr    8.0 8.0 1.75 0.004 6.0 1.0 0.80 0.3 0.25 0.15 150 290 20 30 10.0 0.05 0.55 1.0
    
    #CROP PARAMETERS
    #NCRP: number of crops (WITH IRRIGATION)
    0 #2
    #NFIELD: number of fields and related irrigation time serie
    0 #5
    # to repeat NCROP times in a same line
    # 0 - CropType: name of the crop type [string, 1 word, no space allowed]
    # 1 - h_c: max heigth of plant[m]
    # 2 - S_w_c: canopy capacity [mm]
    # 3 - C_leaf_star_c: maximum leaf conductance [m.s-1]
    # 4 - LAI_c: max leaf area index [m2.m-2]
    # 5 - f_s_c: shelter factor []
    # 6 - alfa_c: vegetation albedo []
    # 7 - Zr_c: maximum root depth [m]
    # 8 - kTg_min_c: transpiration sourcing factor min [], 1>=k_Tu>0
    # 9 - kTg_max: transpiration sourcing factor max [], 1<k_T<=0
    # 10 - kT_f_c: transpiration sourcing factor [], phi>f>wp
    # 11 - kT_s_c: transpiration sourcing factor slope [], 1>s>0
    # CROP1 (alfafa) # CHANGE PARAMETERS WITH FAO56
    #alfafa 1.0 0.15 0.010 3.5 0.55 0.15 0.40 0.01 0.011 2.0
    # CROP2 (wheat)  # CHANGE PARAMETERS WITH FAO56
    #wheat 1.8 0.15 0.010 1.5 0.55 0.35 0.40 0.01 0.011 2.0

    # SOIL PARAMETERS
    # NSOIL: number of soil types
    2
    # to repeat NSOIL times in a same line
    # 0 - SoilType: name of the soil type [string, 1 word, no space allowed]
    # 1 - por: surface (1cm) soil porosity [m3.m-3]
    # 2 - fc: surface (1cm) soil field capacity [m3.m-3]
    # 3 - alfa_sd: soil albedo in dry season
    # 4 - alfa_sw: soil albedo in wet season
    # 5 - J_sd: starting julian day of the dry season [int 1-365]
    # 6 - J_sw: starting julian day of the wet season [int 1-365]
    # 7 - TRANS_sdw: transition period between dry and wet season [days]
    # 8 - TRANS_swd: transition period between wet and dry season [days]
    # SOIL1 (alluvium)
    alluvium 0.40 0.25 0.20 0.15 150 270 20 20
    # SOIL2 (regolith)
    regolith 0.35 0.15 0.25 0.20 150 270 20 20

    # BY DEFAULT Eo (evaporation from open water using the Penman equation) WILL BE COMPUTED
    # OPEN WATER PARAMETERS
    # 1 - alfa_w: water albedo
    # 0.06

    # PLOTTING: positive integer for plotting, 0 or negative integer for NO plotting
    0
    - EOF - don't copy this line in a txt file

    Output:
    -------------------------------
    ASCII files (comma separated) with datetime, julian day,
    Eo/ETref/PTveg/PEsoil and number of data/per day if daily data
'''

###########################################

def MMsurf(cUTIL, pathMMsurf, inputFile_TS_fn, inputFile_PAR_fn, outputFILE_fn, pathMMws, outMMsurf_fn, MMsurf_plot = 0, inputFile_IRR_TS_fn = None):

    global J_d, datenumOUT, input_TS_crop_irr_fn, kT_f_c, kT_s_c, kTg_max_c, kTg_min_c, Zr_c, inputZON_TS_PT_irr_fn, inputZON_TS_TF_irr_fn, inputZON_TS_P_irr_fn, NCROP, kT_f, kT_s, kTg_max, kTg_min, Zr, LAI_veg_d, inputZONTF_irr, inputZONP_irr, inputZONPT_irr, v, alfa_w, fc, por, C_leaf_star_v, S_w_v, z_h, z_m, FC, Lz, Z, Lm, phi, Rs, u_z_m, Pa, RHa, Ta, P, f_s_c, alfa_c, datenum_d, IRR_TS, SoilType, alfa_sw, alfa_sd, J_sw, J_sd, NSOIL, f_s_vw, f_s_vd, alfa_vw, alfa_vd, NVEG, f, data_IRR_TS, NFIELD, datenum, Rs1, u_z_m1, Pa1, RHa1, Ta1, P1, actual_day, datetime_i, datetime, time, date, NMETEO, DTS, dataMETEOTS, line, l

    def ExportResults(name, ws, row1, Dates, J, TS, ts_output = 0, n_d = [], TypeFile = "PET"):
        """
        Write the processed data in a open txt file
        INPUT:      output fluxes time series, date, Julian day, etc... and open file
        """
        global out_line
        outpathname = os.path.join(ws, name)
        outFileExport = open(outpathname, 'w')
        outFileExport.write(row1)
        for t in range(len(Dates)):
            year='%4d'%mpl.dates.num2date(Dates[t]).year
            month='%02d'%mpl.dates.num2date(Dates[t]).month
            day='%02d'%mpl.dates.num2date(Dates[t]).day
            hour='%02d'%mpl.dates.num2date(Dates[t]).hour
            minute='%02d'%mpl.dates.num2date(Dates[t]).minute
            date_t=(year + "-" + month + "-" + day + " " + hour + ":" + minute)
            if TypeFile == "PET":
                if ts_output == 0:
                    out_line =  date_t, ', ', '%3d'% J[t], ",", '%14.9G' %TS[t],'\n'  #mpl.dates.num2date(Dates[t]).isoformat()[:10]
                else:
                    out_line =  date_t, ', ', '%3d'% J[t], ",", '%14.9G' %TS[t], ",",'%2d'% n_d[t],'\n'  #mpl.dates.num2date(Dates[t]).isoformat()[:10]
            elif TypeFile == "VAR":
                out_line =  date_t, ', ', '%3d'% J[t], ",", '%14.9G' %TS[0][t], ',', '%14.9G' %TS[1][t], ',', '%14.9G' %TS[2][t], ',', '%14.9G' %TS[3][t],',', '%14.9G' %TS[4][t],',', '%14.9G' %TS[5][t],',', '%14.9G' %TS[6][t],'\n'
            elif TypeFile == "P":
                out_line =  date_t, ', ', '%3d'% J[t], ",", '%14.9G' %TS[0][t], ',', '%14.9G' %TS[1][t], ',', '%14.9G' %TS[2][t], ',', '%14.9G' %TS[3][t], ',', '%14.9G' %TS[4][t],'\n'
            elif TypeFile == "avP_Eo":
                out_line =  date_t, ', ', '%14.9G' %TS[0][t], ',', '%14.9G' %TS[1][t], ',', '\n'
            elif TypeFile == "Date":
                out_line =  date_t, ', ', '%3d'% J[t], '\n'
            for l in out_line:
                outFileExport.write(l)
        print("\n[" + outpathname + "] exported!")
        outFileExport.close()

    def ExportResults1(TS, outFileExport):
        """
        Write the processed data in a open txt file readable by  MARMITES
        INPUT:      output flux time series and open file
        """
        for i in range(len(TS)):
            out_line =  '%14.9G' %TS[i],'\n'
            for l in out_line:
                outFileExport.write(l)

    ###########################################


    #  ##### READING INPUT ###############################################
    print('\nReading MARMITESsurface INPUT FILES...')

    # METEO/VEGETATION/SOIL/WATER PARAMETERS

    inputFile = cUTIL.readFile(pathMMsurf, inputFile_PAR_fn)
    try:
        # METEO PARAMETERS
        l=0
        line = inputFile[l].split()
        NMETEO = int(line[0])
        if NMETEO < 1:
            cUTIL.ExitError(msg = '\nFATAL ERROR!\nNMETEO must be higher than 1!\n')
        phi = []
        Lm = []
        Z = []
        Lz = []
        FC = []
        DTS = []
        z_m = []
        z_h = []
        for i in range(NMETEO):
            l = l+1
            line = inputFile[l].split()
            phi.append(float(line[0])) # latitude of the meteo station [degres]
            Lm.append(float(line[1])) #  longitude of the meteo station [degres]
            Z.append(float(line[2]))
            Lz.append(float(line[3])) #  longitude of the center of the local time zone [degres] (west of Greenwich)
            FC.append(float(line[4]))
            DTS.append(float(line[5])) # data time shift, for data NOT acquired at standart clock time [h]
            z_m.append(float(line[6])) # heigth of u_z_m measurement [m]
            z_h.append(float(line[7])) # heigth of humidity measurement [m]
        #VEGETATION PARAMETERS
        l = l + 1
        line = inputFile[l].split()
        NVEG = int(line[0]) # number of vegetation types (NO IRRIGATION)
        VegName = [] # name of the vegetation type [string]
        h_vd = [] # heigth of plant dry season [m]
        h_vw = [] # heigth of plant wet season [m]
        S_w_v = [] # canopy capacity [mm]
        C_leaf_star_v = [] # bulk stomatal resistance [s.m-1]
        LAI_vd = [] # leaf area index dry season [m2.m-2]
        LAI_vw = [] # leaf area index wet season [m2.m-2]
        f_s_vd = [] # shelter factor dry season []
        f_s_vw = [] # shelter factor wet season []
        alfa_vd = [] # vegetation albedo in dry season
        alfa_vw = [] # vegetation albedo in wet season
        J_vd = [] # starting julian day of the dry season [int 1-365]
        J_vw = [] # starting julian day of the wet season [int 1-365]
        TRANS_vdw = [] # transition period between dry and wet season [days]
        TRANS_vwd = [] # transition period between wet and dry season [days]
        Zr = [] # maximum root depth [m]
        kTg_min = [] #transpiration sourcing factor min [], 1>=k_Tu>0
        kTg_max = [] #transpiration sourcing factor max [], 1<k_Tu<=0
        kT_f = [] # transpiration sourcing factor
        kT_s = []  # transpiration sourcing factor
        # input FAO56 grass parameters
        VegName.append("grassFAO56")
        h_vd.append(0.12)
        h_vw.append(0.12)
        S_w_v.append(0.1)
        C_leaf_star_v.append(0.01)
        LAI_vd.append(2.88)
        LAI_vw.append(2.88)
        f_s_vd.append(0.5)
        f_s_vw.append(0.5)
        alfa_vd.append(0.23)
        alfa_vw.append(0.23)
        J_vd.append(1)
        J_vw.append(365)
        TRANS_vdw.append(1)
        TRANS_vwd.append(1)
        Zr.append(0.25)
        kTg_min.append(0.05)
        kTg_max.append(0.10)
        kT_f.append(0.5)
        kT_s.append(1.0)
        if NVEG>0:
            for i in range(NVEG):
                l = l + 1
                e = 0
                line = inputFile[l].split()
                VegName.append(str(line[e]))
                e += 1
                h_vd.append(float(line[e]))
                e += 1
                h_vw.append(float(line[e]))
                e += 1
                S_w_v.append(float(line[e]))
                e += 1
                C_leaf_star_v.append(float(line[e]))
                e += 1
                LAI_vd.append(float(line[e]))
                e += 1
                LAI_vw.append(float(line[e]))
                e += 1
                f_s_vd.append(float(line[e]))
                e += 1
                f_s_vw.append(float(line[e]))
                e += 1
                alfa_vd.append(float(line[e]))
                e += 1
                alfa_vw.append(float(line[e]))
                e += 1
                J_vd.append(int(line[e]))
                e += 1
                J_vw.append(int(line[e]))
                e += 1
                TRANS_vdw.append(int(line[e]))
                e += 1
                TRANS_vwd.append(int(line[e]))
                e += 1
                Zr.append(float(line[e]))
                e += 1
                kTg_min.append(float(line[e]))
                e += 1
                kTg_max.append(float(line[e]))
                e += 1
                kT_f.append(float(line[e]))
                e += 1
                kT_s.append(float(line[e]))
                del e
        NVEG = NVEG + 1 # number of vegetation types + FAO56 GRASS (default)
        #CROP PARAMETERS
        l = l + 1
        line = inputFile[l].split()
        NCROP = int(line[0]) # number of crop types (WITH IRRIGATION)
        l = l + 1
        line = inputFile[l].split()
        NFIELD = int(line[0]) # number of fields and related irrigation time serie
        CropType = [] # name of the crop type [string]
        hm_c = [] # heigth of plant [m]
        S_w_c = [] # canopy capacity [mm]
        C_leaf_star_c = [] # bulk stomatal resistance [s.m-1]
        LAIm_c = [] # max leaf area index [m2.m-2]
        f_s_c = [] # shelter factor []
        alfa_c = [] # vegetation albedo in wet season
        Zr_c = [] # maximum root depth [m]
        kTg_min_c = [] #transpiration sourcing factor min [], 1>=k_Tu>0
        kTg_max_c = [] #transpiration sourcing factor min [], 1<k_Tu<=0
        kT_f_c = [] # transpiration sourcing factor
        kT_s_c = []  # transpiration sourcing factor
        if NCROP > 0:
            for i in range(NCROP):
                l = l + 1
                e = 0
                line = inputFile[l].split()
                CropType.append(str(line[e]))
                e += 1
                hm_c.append(float(line[e]))
                e += 1
                S_w_c.append(float(line[e]))
                e += 1
                C_leaf_star_c.append(float(line[e]))
                e += 1
                LAIm_c.append(float(line[e]))
                e += 1
                f_s_c.append(float(line[e]))
                e += 1
                alfa_c.append(float(line[e]))
                e += 1
                Zr_c.append(float(line[e]))
                e += 1
                kTg_min_c.append(float(line[e]))
                e += 1
                kTg_max_c.append(float(line[e]))
                e += 1
                kT_f_c.append(float(line[e]))
                e += 1
                kT_s_c.append(float(line[e]))
                del e
        # SOIL PARAMETERS
        l = l + 1
        line = inputFile[l].split()
        NSOIL = int(line[0]) # number of soil types
        SoilType = [] # name of the soil type [string]
        por = [] # surface (1cm) soil porosity [m3.m-3]
        fc = [] # surface (1cm) soil field capacity [m3.m-3]
        alfa_sd = [] # soil albedo in dry season
        alfa_sw = [] # soil albedo in wet season
        J_sd = [] # starting julian day of the dry season [int 1-365]
        J_sw = [] # starting julian day of the wet season [int 1-365]
        TRANS_sdw = [] # transition period between dry and wet season [days]
        if NSOIL>0:
            for i in range(NSOIL):
                l = l + 1
                line = inputFile[l].split()
                SoilType.append(str(line[0]))
                por.append(float(line[1]))
                fc.append(float(line[2]))
                alfa_sd.append(float(line[3]))
                alfa_sw.append(float(line[4]))
                J_sd.append(int(line[5]))
                J_sw.append(int(line[6]))
                TRANS_sdw.append(int(line[7]))
        # WATER PARAMETERS
        alfa_w = 0.06 # water albedo
    except:
        cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nError reading the input file:\n[" + inputFile_PAR_fn +"]")
    del inputFile, l, line
    print("\nMETEO/VEGETATION/CROP/SOIL PARAMETERS file imported!\n[" + inputFile_PAR_fn +"]")

    # read INPUT file
    # METEO TIME SERIES
    inputFile_TS_fn = os.path.join(pathMMsurf, inputFile_TS_fn)
    if os.path.exists(inputFile_TS_fn):
        dataMETEOTS = np.loadtxt(inputFile_TS_fn, skiprows = 1, dtype = str)
    else:
        cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe input file [" + inputFile_TS_fn + "] doesn't exist, verify name and path!")
    try:
        date = dataMETEOTS[:,0]
        time = dataMETEOTS[:,1]
        datenum = []
        datetime = []
        datenum_d = []
        datetime_i = (date[0] + " " + time[0])
        datenum_i = mpl.dates.datestr2num(datetime_i)-(DTS[0])/24.0
        actual_day = mpl.dates.num2date(datenum_i).isoformat()[:10]
        datenum_d.append(mpl.dates.datestr2num(actual_day))
        P1 = []
        Ta1 = []
        RHa1 = []
        Pa1 = []
        u_z_m1 = []
        Rs1 = []
        P = []
        Ta = []
        RHa = []
        Pa = []
        u_z_m = []
        Rs = []
        for n in range(NMETEO):
            P1.append(dataMETEOTS[:,2+n*6])
            Ta1.append(dataMETEOTS[:,3+n*6])
            RHa1.append(dataMETEOTS[:,4+n*6])
            Pa1.append(dataMETEOTS[:,5+n*6])
            u_z_m1.append(dataMETEOTS[:,6+n*6])
            Rs1.append(dataMETEOTS[:,7+n*6])
            P.append([])
            Ta.append([])
            RHa.append([])
            Pa.append([])
            u_z_m.append([])
            Rs.append([])
            for t in range(len(Ta1[n])):
                if n == 0:
                    datetime.append(date[t] + " " + time[t])
                    datenum.append(mpl.dates.datestr2num(datetime[t])-(DTS[0])/24.0)
                    if actual_day != date[t]:
                        datenum_d.append(mpl.dates.datestr2num(date[t]))
                        actual_day = date[t]
                P[n].append(float(P1[n][t]))
                Ta[n].append(float(Ta1[n][t]))
                RHa[n].append(float(RHa1[n][t]))
                Pa[n].append(float(Pa1[n][t]))
                u_z_m[n].append(float(u_z_m1[n][t]))
                Rs[n].append(float(Rs1[n][t]))
        datenum_d = np.asarray(datenum_d)
        P = np.asarray(P)
        Ta = np.asarray(Ta)
        RHa = np.asarray(RHa)
        Pa = np.asarray(Pa)
        u_z_m = np.asarray(u_z_m)
        Rs = np.asarray(Rs)
    except:
        cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nUnexpected error in the input file\n[%s]\n" % inputFile_TS_fn)
    del dataMETEOTS, date, time, datetime, datetime_i, actual_day, P1, Ta1, RHa1, Pa1, u_z_m1,Rs1
    print("\nMETEO TIME SERIES file imported!\n[" + inputFile_TS_fn +"]")
    # compute J - julian day
    YYYY = []
    MM = []
    DD = []
    HH = []
    MN = []
    J = []
    time = []
    for i, d in enumerate(datenum):
        YYYY.append(float(mpl.dates.num2date(d).year))
        MM.append(float(mpl.dates.num2date(d).month))
        DD.append(float(mpl.dates.num2date(d).day))
        HH.append('%02d'%mpl.dates.num2date(d).hour)
        MN.append('%02d'%mpl.dates.num2date(d).minute)
        J.append(DD[i] - 32 + int(275*MM[i]/9) + 2 * int(3/(MM[i] + 1)) + int(MM[i]/100-np.mod(YYYY[i],4)/4+0.975) )
        time.append('%s:%s' % (HH[i], MN[i]))
    J = np.asarray(J)
    time = np.asarray(time)
    del YYYY, MM, DD, HH, MN
    # IRRIGATION
    FIELD_startdate = []
    FIELD_enddate = []
    FIELD_growdur = []
    FIELD_wiltdur = []
    FIELD_crop = []
    if inputFile_IRR_TS_fn != None:
        # IRRIGATION TIME SERIES
        inputFile_IRR_TS_fn = os.path.join(pathMMsurf, inputFile_IRR_TS_fn)
        if os.path.exists(inputFile_IRR_TS_fn):
            data_IRR_TS = np.loadtxt(inputFile_IRR_TS_fn, skiprows = 1, dtype = str)
        else:
            cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe input file [%s] doesn't exist, verify name and path!"%inputFile_IRR_TS_fn)
        try:
            IRR_TS = []
            for n in range(NFIELD):
                IRR_TS.append(data_IRR_TS[:,2+n].astype(float))
        except:
            cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nUnexpected error in the input file\n[%s]\n" % inputFile_IRR_TS_fn)
        print("\nIRRIGATION TIME SERIES file imported!\n[" + inputFile_IRR_TS_fn +"]")
        # FIELD/CROP schedule
        FIELD_crop_schedule_fn = []
        FIELD_crop_schedule = []
        for f in range(NFIELD):
            FIELD_crop_schedule_fn.append(os.path.join(pathMMsurf, '__inputFIELD%d_crop_schedule.txt' % (f+1)))
            if os.path.exists(FIELD_crop_schedule_fn[f]):
                FIELD_crop_schedule.append(np.loadtxt(FIELD_crop_schedule_fn[f], skiprows = 1, dtype = str))
            else:
                cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe input file [%s] doesn't exist, verify name and path!"%FIELD_crop_schedule_fn[f])
        print("\nImporting FIELD/CROP schedule files...")
        try:
            for f in range(NFIELD):
                FIELD_startdate.append(FIELD_crop_schedule[f][:,0])
                FIELD_enddate.append(FIELD_crop_schedule[f][:,1])
                FIELD_growdur.append(FIELD_crop_schedule[f][:,2])
                FIELD_wiltdur.append(FIELD_crop_schedule[f][:,3])
                FIELD_crop.append(FIELD_crop_schedule[f][:,4])
                print("[%s]" % FIELD_crop_schedule_fn[f])
                for i, d in enumerate(FIELD_startdate[f]):
                    try:
                        FIELD_startdate[f][i] = mpl.dates.datestr2num(d)
                        FIELD_enddate[f][i] = mpl.dates.datestr2num(FIELD_enddate[f][i])
                        if FIELD_startdate[f][i] > FIELD_enddate[f][i]:
                            cUTIL.ErrorExit(msg = '\nFATAL ERROR!\nDates in the field/crop schedule files of field #%d are not correct!'%(i+1))
                        if i > 0:
                            if FIELD_enddate[f][i-1] > FIELD_startdate[f][i]:
                                cUTIL.ErrorExit(msg = '\nFATAL ERROR!\nDates in the field/crop schedule files of field #%d are not correct!'%(i+1))
                    except:
                        cUTIL.ErrorExit(msg = '\nFATAL ERROR!\nDates in the field/crop schedule files of field #%d are not correct!'%(i+1))
                    FIELD_growdur[f][i] = int(FIELD_growdur[f][i])
                    FIELD_wiltdur[f][i] = int(FIELD_wiltdur[f][i])
                FIELD_startdate[f] = FIELD_startdate[f].astype(float)
                FIELD_enddate[f] = FIELD_enddate[f].astype(float)
                FIELD_growdur[f] = FIELD_growdur[f].astype(float)
                FIELD_wiltdur[f] = FIELD_wiltdur[f].astype(float)
                FIELD_crop[f] = FIELD_crop[f].astype(int)
        except:
            cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nUnexpected error in the input file\n[%s]\n"%FIELD_crop_schedule_fn[f])
        print('Done!')
    else:
        IRR_TS = None

    ######################

    def paramTS(J, J_d, J_w, TRANS_dw, TRANS_wd, value_d, value_w, min_, max_, name = 'paramTSwd'):

        if value_d == value_w:
            TS = np.ones(len(J), dtype = float)*value_d
        else:
            TS = np.zeros(np.asarray(J).shape, dtype = float)
            trans_wd = (value_d - value_w)/TRANS_wd
            trans_dw = (value_w - value_d)/TRANS_dw
            for j in range(len(J)):
                if J[j] <= J_d or J[j] >= J_w:
                    # wet period
                    TS[j] = value_w
                else:
                    if J[j] < (J_d + TRANS_wd):
                        # transition wet to dry
                        TS[j] =  value_w + (J[j] - J_d)*trans_wd
                    elif J[j] > (J_w - TRANS_dw):
                        # transition dry to wet
                        TS[j] =  value_d + (J[j] - (J_w - TRANS_dw))*trans_dw
                    else:
                        # dry period
                        TS[j] =  value_d

        monthsFmt=mpl.dates.DateFormatter('%Y-%m-%d')
        plt_export_fn = os.path.join(pathMMsurf, '%s.png' % name)
        fig = plt.figure()
        ax1=fig.add_subplot(111)
        plt.plot_date(datenum,TS, 'b-')
        plt.yticks(np.linspace(min_, max_, 11))
        plt.ylim(min_*0.95,max_*1.05)
        plt.setp(ax1.get_yticklabels(), fontsize=8)
        plt.setp(ax1.get_xticklabels(), fontsize=8)
        plt.ylabel(name, fontsize=10)
        ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
        ax1.xaxis.set_major_formatter(monthsFmt)
        plt.grid(True)
        plt.savefig(plt_export_fn,dpi=150)

        return TS

    ######################

    def paramTS_f(datenum, d_start, d_end, TRANS_start, TRANS_end, crop, value, min_, max_, name = 'paramTSwd', value_min = 0.0):

        TS = np.ones((len(datenum)), dtype = float)*value_min
        i = 0
        while datenum[0] > d_start[i] and datenum[0] >= d_end[i]:
            i += 1
        for t, d in enumerate(datenum):
            trans_nc2c = (value[crop[i]-1] - value_min)/TRANS_start[i]
            trans_c2nc = (value_min - value[crop[i]-1])/TRANS_end[i]
            if d <= d_start[i] or d >= d_end[i]:
                # no crop
                TS[t] = value_min
                if d == d_end[i]:
                    i += 1
                    if i == len(crop):
                        break
            else:
                # crop
                if d < (d_start[i] + TRANS_start[i]):
                    # transition no crop to crop
                    TS[t] = value_min + (d - d_start[i])*trans_nc2c
                elif d > (d_end[i] - TRANS_end[i]):
                    # transition crop to no crop
                    TS[t] = value[crop[i]-1] + (d - (d_end[i] - TRANS_end[i]))*trans_c2nc
                else:
                    # full crop
                    TS[t] = value[crop[i]-1]

        monthsFmt=mpl.dates.DateFormatter('%Y-%m-%d')
        plt_export_fn = os.path.join(pathMMsurf, '%s.png' % name)
        fig = plt.figure()
        ax1=fig.add_subplot(111)
        plt.plot_date(datenum,TS, 'b-')
        plt.yticks(np.linspace(min_, max_, 11))
        plt.ylim(min_*0.95,max_*1.05)
        plt.setp(ax1.get_yticklabels(), fontsize=8)
        plt.setp(ax1.get_xticklabels(), fontsize=8)
        plt.ylabel(name, fontsize=10)
        ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
        ax1.xaxis.set_major_formatter(monthsFmt)
        plt.grid(True)
        plt.savefig(plt_export_fn,dpi=150)

        return TS

    ######################

    def paramTS_f_c(datenum, d_start, d_end, crop, name = 'paramTScrop'):

        TS = np.zeros((len(datenum)), dtype = float)
        i = 0
        while datenum[0] > d_start[i] and datenum[0] >= d_end[i]:
            i += 1
        for t, d in enumerate(datenum):
            if d <= d_start[i] or d >= d_end[i]:
                # no crop
                TS[t] = 0.0
            else:
                # crop
                TS[t] = float(crop[i])
            if d == d_end[i]:
                i += 1
                if i == len(crop):
                    break

        monthsFmt=mpl.dates.DateFormatter('%Y-%m-%d')
        plt_export_fn = os.path.join(pathMMsurf, '%s.png' % name)
        fig = plt.figure()
        ax1=fig.add_subplot(111)
        plt.plot_date(datenum,TS, 'b-')
        plt.setp(ax1.get_yticklabels(), fontsize=8)
        plt.setp(ax1.get_xticklabels(), fontsize=8)
        plt.ylim(0,max(crop)+1)
        plt.ylabel(name, fontsize=10)
        ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
        ax1.xaxis.set_major_formatter(monthsFmt)
        plt.grid(True)
        plt.savefig(plt_export_fn,dpi=150)

        return TS

    ######################
    print('\nComputing vegetation parameters time series...')
    alfa_v = np.zeros((NVEG, len(datenum)), dtype = float)
    f_s_v = np.zeros((NVEG, len(datenum)), dtype = float)
    LAI_v = np.zeros((NVEG, len(datenum)), dtype = float)
    h_v = np.zeros((NVEG, len(datenum)), dtype = float)
    for v in range(NVEG):
        alfa_v[v] = paramTS(J = J, J_d = J_vd[v], J_w = J_vw[v], TRANS_dw = TRANS_vdw[v], TRANS_wd = TRANS_vwd[v], value_d = alfa_vd[v], value_w = alfa_vw[v], min_ = 0.0, max_ = 1.0, name = 'albedo_VEG%d_%s' % (v, VegName[v]))
        f_s_v[v] = paramTS(J = J, J_d = J_vd[v], J_w = J_vw[v], TRANS_dw = TRANS_vdw[v], TRANS_wd = TRANS_vwd[v], value_d = f_s_vd[v], value_w = f_s_vw[v], min_ = 0.0, max_ = 1.0, name = 'f_s_VEG%d_%s' % (v, VegName[v]))
        LAI_v[v] = paramTS(J = J, J_d = J_vd[v], J_w = J_vw[v], TRANS_dw = TRANS_vdw[v], TRANS_wd = TRANS_vwd[v], value_d = LAI_vd[v], value_w = LAI_vw[v], min_ = 0.0, max_ = max(max(LAI_vd), max(LAI_vw)), name = 'LAI_VEG%d_%s' % (v, VegName[v]))
        h_v[v] = paramTS(J = J, J_d = J_vd[v], J_w = J_vw[v], TRANS_dw = TRANS_vdw[v], TRANS_wd = TRANS_vwd[v], value_d = h_vd[v], value_w = h_vw[v], min_ = 0.0, max_ = max(h_vd[v], h_vw[v]), name = 'heigth_VEG%d_%s' % (v, VegName[v]))

    alfa_s = np.zeros((NSOIL, len(datenum)), dtype = float)
    for s in range(NSOIL):
        alfa_s[s] = paramTS(J = J, J_d = J_sd[s], J_w = J_sw[s], TRANS_dw = TRANS_sdw[s], TRANS_wd = TRANS_sdw[s], value_d = alfa_sd[s], value_w = alfa_sw[s], min_ = 0.0, max_ = 1.0, name = 'albedo_SOIL%d_%s' % (s+1, SoilType[s]))

    if IRR_TS != None:
        print('\nComputing field crops parameters time series...')
        alfa_f = np.zeros((NFIELD, len(datenum)), dtype = float)
        f_s_f = np.zeros((NFIELD, len(datenum)), dtype = float)
        LAI_f = np.zeros((NFIELD, len(datenum)), dtype = float)
        C_leaf_star_f = np.zeros((NFIELD, len(datenum)), dtype = float)
        h_f = np.zeros((NFIELD, len(datenum)), dtype = float)
        S_w_f = np.zeros((NFIELD, len(datenum_d)), dtype = float)
        crop_f = np.zeros((NFIELD, len(datenum_d)), dtype = float)
        for f in range(NFIELD):
            alfa_f[f] = paramTS_f(datenum = datenum, d_start = FIELD_startdate[f], d_end = FIELD_enddate[f], TRANS_start = FIELD_growdur[f], TRANS_end = FIELD_wiltdur[f], value = alfa_c, crop = FIELD_crop[f], min_ = 0.0, max_ = 1.0, name = 'albedo_FIELD%s' % (f+1))
            f_s_f[f] = paramTS_f(datenum = datenum, d_start = FIELD_startdate[f], d_end = FIELD_enddate[f], TRANS_start = FIELD_growdur[f], TRANS_end = FIELD_wiltdur[f], value = f_s_c, crop = FIELD_crop[f], min_ = 0.0, max_ = 1.0, name = 'f_s_FIELD%s' % (f+1))
            LAI_f[f] = paramTS_f(datenum = datenum, d_start = FIELD_startdate[f], d_end = FIELD_enddate[f], TRANS_start = FIELD_growdur[f], TRANS_end = FIELD_wiltdur[f], value = LAIm_c, crop = FIELD_crop[f], min_ = 0.0, max_ = max(LAIm_c), name = 'LAI_FIELD%s' % (f+1))
            C_leaf_star_f[f] = paramTS_f(datenum = datenum, d_start = FIELD_startdate[f], d_end = FIELD_enddate[f], TRANS_start = FIELD_growdur[f], TRANS_end = FIELD_wiltdur[f], value = C_leaf_star_c, crop = FIELD_crop[f], min_ = 0.0, max_ = max(C_leaf_star_c), name = 'C_leaf_FIELD%s' % (f+1))
            h_f[f] = paramTS_f(datenum = datenum, d_start = FIELD_startdate[f], d_end = FIELD_enddate[f], TRANS_start = FIELD_growdur[f], TRANS_end = FIELD_wiltdur[f], value = hm_c, crop = FIELD_crop[f], min_ = 0.0, max_ = max(hm_c), name = 'heigth_FIELD%s' % (f+1), value_min = 0.01)
            S_w_f[f] = paramTS_f(datenum = datenum_d, d_start = FIELD_startdate[f], d_end = FIELD_enddate[f], TRANS_start = FIELD_growdur[f], TRANS_end = FIELD_wiltdur[f], value = S_w_c, crop = FIELD_crop[f], min_ = 0.0, max_ = max(S_w_c), name = 'S_w_FIELD%s' % (f+1))
            crop_f[f]  = paramTS_f_c(datenum = datenum_d, d_start = FIELD_startdate[f], d_end = FIELD_enddate[f], crop = FIELD_crop[f], name = 'crop_FIELD%s' % (f+1))
    else:
        alfa_f = []
        f_s_f = []
        LAI_f = []
        C_leaf_star_f = []
        h_f = []
        S_w_f = []
        crop_f = []

    #  ##### COMPUTING PT, PE and INTERCEPTION ##############################################
    print("\nComputing PT, PE, TF, etc...")

    inputZON_TS_P_veg_fn = "inputZONP_veg_d.txt"
    inputZONP_veg_fn = os.path.join(pathMMws, inputZON_TS_P_veg_fn)
    inputZONP_veg = open(inputZONP_veg_fn, 'w')
    inputZONP_veg.write('#\n')

    inputZON_TS_TF_veg_fn = "inputZONTF_veg_d.txt"
    inputZONTF_veg_fn = os.path.join(pathMMws, inputZON_TS_TF_veg_fn)
    inputZONTF_veg = open(inputZONTF_veg_fn, 'w')
    inputZONTF_veg.write('#\n')

    inputZON_TS_PT_veg_fn = "inputZONPT_veg_d.txt"
    inputZONPT_veg_fn = os.path.join(pathMMws, inputZON_TS_PT_veg_fn)
    inputZONPT_veg = open(inputZONPT_veg_fn, 'w')
    inputZONPT_veg.write('#\n')

    input_TS_LAI_veg_fn = "inputZONLAI_veg_d.txt"
    inputLAI_veg_fn = os.path.join(pathMMws, input_TS_LAI_veg_fn)
    inputLAI_veg = open(inputLAI_veg_fn, 'w')
    inputLAI_veg.write('#\n')

    inputZON_TS_PE_fn = "inputZONPE_d.txt"
    inputZONPE_fn = os.path.join(pathMMws, inputZON_TS_PE_fn)
    inputZONPE = open(inputZONPE_fn, 'w')
    inputZONPE.write('#\n')

    inputZON_TS_Eo_fn = "inputZONEo_d.txt"
    inputZONEo_fn = os.path.join(pathMMws, inputZON_TS_Eo_fn)
    inputZONEo = open(inputZONEo_fn, 'w')
    inputZONEo.write('#\n')

    if IRR_TS != None:
        inputZON_TS_P_irr_fn = "inputZONP_irr_d.txt"
        inputZONP_irr_fn = os.path.join(pathMMws, inputZON_TS_P_irr_fn)
        inputZONP_irr = open(inputZONP_irr_fn, 'w')
        inputZONP_irr.write('#\n')

        inputZON_TS_PT_irr_fn = "inputZONPT_irr_d.txt"
        inputZONPT_irr_fn = os.path.join(pathMMws, inputZON_TS_PT_irr_fn)
        inputZONPT_irr = open(inputZONPT_irr_fn, 'w')
        inputZONPT_irr.write('#\n')

        inputZON_TS_TF_irr_fn = "inputZONTF_irr_d.txt"
        inputZONTF_irr_fn = os.path.join(pathMMws, inputZON_TS_TF_irr_fn)
        inputZONTF_irr = open(inputZONTF_irr_fn, 'w')
        inputZONTF_irr.write('#\n')

        input_TS_crop_irr_fn = "inputZONcrop_irr_d.txt"
        inputcrop_irr_fn = os.path.join(pathMMws, input_TS_crop_irr_fn)
        inputcrop_irr = open(inputcrop_irr_fn, 'w')
        inputcrop_irr.write('#\n')
        for f in range(NFIELD):
            ExportResults1(crop_f[f], inputcrop_irr)
        inputcrop_irr.close()

    for n in range(NMETEO):
        print("\n--------------\nProcessing data of ZONE %d/%d\n--------------" % (n+1, NMETEO))
        J, J_d, outputVAR, PT_PM_VEG, Erf_VEG, PE_PM_SOIL, Eo, PT_PM_VEG_d, PE_PM_SOIL_d, Eo_d, P_veg_d, Pint_veg,                      P_veg_duration, n1_d, TF_veg_d, I_veg_d, LAI_veg_d, PT_PM_FIELD, PT_PM_FIELD_d, Erf_FIELD, P_irr_d, Pint_irr, P_irr_duration, TF_irr_d, I_irr_d = PET_P_INTER.process(
                cUTIL,
                datenum, datenum_d, J, time, pathMMsurf,\
                P[n], IRR_TS, Ta[n], RHa[n], Pa[n], u_z_m[n], Rs[n], \
                phi[n], Lm[n], Z[n], Lz[n], FC[n], z_m[n], z_h[n], \
                NVEG, VegName, S_w_v, C_leaf_star_v, alfa_v, f_s_v, LAI_v, h_v,\
                NSOIL, SoilType, por, fc, alfa_s,\
                alfa_w,\
                NFIELD ,alfa_f, f_s_f, LAI_f, C_leaf_star_f, h_f, S_w_f
                )

        #  #####  PLOTTING ##############################################
        print("\nExporting plots...")
        npdatenum = np.array(datenum)
        if MMsurf_plot == 1:
            plt.switch_backend('WXagg')
        # outputVAR = [DELTA,gama, Rs, Rs_corr, Rs0, Rnl, r]
        plot_exportVAR_fn = os.path.join(pathMMsurf,  outputFILE_fn + "_ZON" + str(n+1)+'_VAR.png')
        plotPET.plotVAR(strTitle = 'Penman-Monteith variables', x = npdatenum \
                ,y5 = PT_PM_VEG[0], y4 = Eo \
                ,y2 = outputVAR[2], y3=outputVAR[3] \
                ,y1 = outputVAR[4], y6 = outputVAR[5]\
                ,lbl_y5 = 'ETref', lbl_y4 = 'Eo' \
                ,lbl_y2 = 'Rs_meas', lbl_y3 = 'Rs_corr' \
                ,lbl_y1 = 'Rs0', lbl_y6 = 'Rnl'
                ,plot_exportVAR_fn = plot_exportVAR_fn
                , MMsurf_plot = MMsurf_plot
                )
        # PLOTTING PT and Eo
        plot_exportPT_fn = os.path.join(pathMMsurf,  outputFILE_fn + "_ZON" + str(n+1)+'_PT.png')
        plotPET.plot(x = datenum_d \
            ,y1 = PT_PM_VEG_d, y2 = Eo_d \
            ,lbl_y1 = VegName, lbl_y2 = 'Eo'
            ,plot_exportPET_fn = plot_exportPT_fn
            , MMsurf_plot = MMsurf_plot\
            ,strTitle = 'PT'
            )

        # PLOTTING PE and Eo
        plot_exportPE_fn = os.path.join(pathMMsurf,  outputFILE_fn + "_ZON" + str(n+1)+'_PE.png')
        plotPET.plot(x = datenum_d \
            ,y1 = PE_PM_SOIL_d, y2 = Eo_d \
            ,lbl_y1 = SoilType, lbl_y2 = 'Eo'
            ,plot_exportPET_fn = plot_exportPE_fn
            , MMsurf_plot = MMsurf_plot\
            ,strTitle = 'PE'
            )

        # PLOTTING P and INTERCEPTION
        # VEG
        plot_exportP_fn = os.path.join(pathMMsurf,  '%s_ZON%s_P_veg.png' % (outputFILE_fn,str(n+1)))
        plotP.plot(x = datenum_d \
                ,y1 = P_veg_d, y2 = Pint_veg, y3 = I_veg_d, y4 = TF_veg_d
                ,lbl_y1 = 'P (mm/d)',  lbl_y2 = 'Pint (mm/h/d)' \
                ,lbl_y3 = 'I (mm/d)',lbl_y4 = 'TF (mm/d)', lbl_veg = VegName\
                ,plot_exportP_fn = plot_exportP_fn
                , MMsurf_plot = MMsurf_plot\
                ,strTitle = 'Rainfall and interception'
                )
        # FIELD/CROP
        if IRR_TS != None:
            # PT
            lbl_f = []
            for f in range(NFIELD): lbl_f.append('FIELD%s'%(f+1))
            plot_exportPET_fn = os.path.join(pathMMsurf,  outputFILE_fn + "_ZON" + str(n+1)+'_PETfields.png')
            datenum_d = np.array(datenum_d)
            plotPET.plot(x = datenum_d \
                ,y1 = PT_PM_FIELD_d\
                ,lbl_y1 = lbl_f
                ,plot_exportPET_fn = plot_exportPET_fn
                , MMsurf_plot = MMsurf_plot\
                , strTitle = 'PT'
                )
            #P/I
            for f in range(NFIELD):
                plot_exportP_fn = os.path.join(pathMMsurf, '%s_ZON%s_P_irr_FIELD%d.png' % (outputFILE_fn,str(n+1),f+1))
                plotP.plot(x = datenum_d \
                        ,y1 = P_irr_d[f], y2 = Pint_irr[f] \
                        ,y3 = [I_irr_d[f]], y4 = [TF_irr_d[f]] \
                        ,lbl_y1 = 'P (mm/d)',  lbl_y2 = 'Pint (mm/h/d)' \
                        ,lbl_y3 = 'I (mm/d)',lbl_y4 = 'TF (mm/d)', lbl_veg = ['crops'] \
                        ,plot_exportP_fn = plot_exportP_fn
                        , MMsurf_plot = MMsurf_plot\
                        , strTitle = 'Rainfall and interception'
                        )

        if MMsurf_plot == 1:
            plt.switch_backend('agg')

        #  #####  EXPORTING ASCII FILES ##############################################
        try:
            if os.path.exists(pathMMsurf):
                print("\n-----------\nExporting output to ASCII files...")
        # EXPORTING RESULTS PT/ET0/PE/Eo/P/TF/INTER/
                ws = pathMMsurf
                for ts_output in range(2):
                    if ts_output == 0:
                        print("\nHourly values...")                  # hourly output
                        datenumOUT = datenum
                        dORh_suffix = "_h"  # for file naming in export
                    else:                                # daily output
                        print("\nDaily values...")
                        datenumOUT = datenum_d
                        dORh_suffix = "_d"  # for file naming in export
                        ExportResults1(P_veg_d, inputZONP_veg)
                        ExportResults1(Eo_d, inputZONEo)
                    outputfile_fn = outputFILE_fn + "_ZON" + str(n+1) + dORh_suffix
                    if ts_output != 1:
                        # OUPUT VARIABLES FOR PM EQUATION
                        name = outputfile_fn + "_VAR.int"
                        row1 = 'Date,J,DELTA,gama,Rs_meas,Rs_corr,Rs0,Rnl,r\n'
                        ExportResults(name, ws, row1, datenumOUT, J, outputVAR, TypeFile = "VAR")
                    # VEG
                    for v in range(NVEG):
                        nam_tmp = 'VEG'+ str(v)
                        if ts_output == 1:
                        # DAILY
                            # PT
                            name = outputfile_fn + "_PT_" + nam_tmp + "_" + VegName[v] + ".out"
                            row1 = 'Date,J,PT,n_d\n'
                            ExportResults(name, ws, row1, datenumOUT, J_d, PT_PM_VEG_d[v],ts_output, n1_d)
                            # P, TF, etc...
                            name = outputfile_fn + "_P_VEG" + str(v) + "_" + VegName[v] + ".out"
                            row1 = 'Date,J,P_mm,duration_day,Pint_mm_h,TF_mm,Interception_mm\n'
                            TStmp = [P_veg_d, P_veg_duration, Pint_veg, TF_veg_d[v], I_veg_d[v]]
                            ExportResults(name, ws, row1, datenumOUT, J_d, TStmp, TypeFile = "P")
                            # PT, TF and LAI for MARMITES
                            if v > 0: #to skip ETref (NVEG=0)
                                ExportResults1(PT_PM_VEG_d[v], inputZONPT_veg)
                                ExportResults1(TF_veg_d[v], inputZONTF_veg)
                        else:
                        # HOURLY
                            # PT
                            name = outputfile_fn + "_PT_" + nam_tmp + "_" + VegName[v] + ".out"
                            row1 = 'Date,J,PT\n'
                            ExportResults(name, ws, row1, datenumOUT, J, PT_PM_VEG[v], ts_output)
                            # Erf
                            name = outputfile_fn + "_Erf_" + nam_tmp + "_" + VegName[v] + ".int"
                            row1 = 'Date,J,Erf\n'
                            ExportResults(name, ws, row1, datenumOUT, J, Erf_VEG[v], ts_output)
                    # NFIELD/NCROP
                    if IRR_TS != None:
                    # P, TF, etc
                        for f in range(NFIELD):
                            if ts_output == 1:
                            # DAILY
                                # PT
                                name = "%s_PT_FIELD%d.out" % (outputfile_fn, f+1)
                                row1 = 'Date,J,PT,n_d\n'
                                ExportResults(name, ws, row1, datenumOUT, J_d, PT_PM_FIELD_d[f],ts_output, n1_d)
                                # PT for MARMITES
                                if v>0: #to skip ETref (NVEG=0)
                                    ExportResults1(PT_PM_FIELD_d[f], inputZONPT_irr)
                                # P, TF, etc...
                                name = "%s_P_FIELD%d.out" % (outputfile_fn, f+1)
                                row1 = 'Date,J,P_mm,duration_day,Pint_mm_h,TF_mm,Interception_mm\n'
                                TStmp = [P_irr_d[f], P_irr_duration[f], Pint_irr[f], TF_irr_d[f], I_irr_d[f]]
                                ExportResults(name, ws, row1, datenumOUT, J_d, TStmp,  TypeFile = "P")
                                # P for MARMITES
                                ExportResults1(P_irr_d[f], inputZONP_irr)
                                # TF for MARMITES
                                ExportResults1(TF_irr_d[f], inputZONTF_irr)
                            else:
                            # HOURLY
                                # PT
                                name = "%s_PT_FIELD%d.out" % (outputfile_fn, f+1)
                                row1 = 'Date,J,PT\n'
                                ExportResults(name, ws, row1, datenumOUT, J, PT_PM_FIELD[f], ts_output)
                                # Erf
                                name = "%s_Erf_FIELD%d.int" % (outputfile_fn, f+1)
                                row1 = 'Date,J,Erf\n'
                                ExportResults(name, ws, row1, datenumOUT, J, Erf_FIELD[f], ts_output)
                    # SOIL
                    if NSOIL>0:
                        for s in range(NSOIL):
                            name = outputfile_fn + "_PE_SOIL" + str(s+1) + "_" + SoilType[s] + ".out"
                            if ts_output == 1:
                                # PE
                                row1 = 'Date,J,PE,n_d\n'
                                ExportResults(name, ws, row1, datenumOUT, J_d, PE_PM_SOIL_d[s], ts_output,n1_d)
                                # PE for MARMITES
                                ExportResults1(PE_PM_SOIL_d[s], inputZONPE)
                            else:
                                row1 = 'Date,J,PE\n'
                                ExportResults(name, ws, row1, datenumOUT, J, PE_PM_SOIL[s], ts_output)
                    # WATER
                    name = outputfile_fn + "_Eo.out"
                    if ts_output == 1:
                        row1 = 'Date,J,Eo,n_d\n'
                        ExportResults(name, ws, row1, datenumOUT, J_d, Eo_d, ts_output,n1_d)
                    else:
                        row1 = 'Date,J,Eo\n'
                        ExportResults(name, ws, row1, datenumOUT, J, Eo, ts_output)
            else:
                print("\nNo ASCII file export required.")
        except:
            cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nError in exporting output files of MMsurf.")

    for v in range(1,NVEG):
        ExportResults1(LAI_veg_d[v], inputLAI_veg)
    inputLAI_veg.close()
    inputZONP_veg.close()
    inputZONTF_veg.close()
    inputZONPT_veg.close()
    inputZONPE.close()
    inputZONEo.close()
    if IRR_TS != None:
        inputZONP_irr.close()
        inputZONTF_irr.close()
        inputZONPT_irr.close()

    # EXPORTING ZONEVEGSOILfile
    date_fn = 'inputDATE.txt'
    outpathname = os.path.join(pathMMws, outMMsurf_fn)
    outFileExport = open(outpathname, 'w')
    outFileExport.write('#\n# NMETEO: number of meteo zones\n')
    outFileExport.write(str(NMETEO))
    outFileExport.write('\n# NVEG: number of vegetation types\n')
    outFileExport.write(str(NVEG-1))
    outFileExport.write('\n# NSOIL: number of soil types\n')
    outFileExport.write(str(NSOIL))
    outFileExport.write('\n# inputDate_fn: file with the dates\n')
    outFileExport.write(date_fn)
    outFileExport.write('\n# inputZON_TS_P_veg_fn: P zones\n')
    outFileExport.write(inputZON_TS_P_veg_fn)
    outFileExport.write('\n# inputZON_TS_TF_veg_fn: TF zones\n')
    outFileExport.write(inputZON_TS_TF_veg_fn)
    outFileExport.write('\n# inputZON_TS_PT_veg_fn: PT zones\n')
    outFileExport.write(inputZON_TS_PT_veg_fn)
    outFileExport.write('\n# input_TS_LAI_veg_fn: LAI veg time series\n')
    outFileExport.write(input_TS_LAI_veg_fn)
    outFileExport.write('\n# inputZON_TS_PE_fn: PE zones\n')
    outFileExport.write(inputZON_TS_PE_fn)
    outFileExport.write('\n# inputZON_TS_Eo_fn: Eo zones\n')
    outFileExport.write(inputZON_TS_Eo_fn)
    outFileExport.write('\n# VegName\n')
    for v in range(1,NVEG):
        outFileExport.write(VegName[v]+' ')
    outFileExport.write('\n# Zr\n')
    for v in range(1,NVEG):
        outFileExport.write(str(Zr[v])+' ')
    outFileExport.write('\n# kTg_min\n')
    for v in range(1,NVEG):
        outFileExport.write(str(kTg_min[v])+' ')
    outFileExport.write('\n# kTg_max\n')
    for v in range(1,NVEG):
        outFileExport.write(str(kTg_max[v])+' ')
    outFileExport.write('\n# kT_f\n')
    for v in range(1,NVEG):
        outFileExport.write(str(kT_f[v])+' ')
    outFileExport.write('\n# kT_s\n')
    for v in range(1,NVEG):
        outFileExport.write(str(kT_s[v])+' ')
    if IRR_TS != None:
        outFileExport.write('\n# NCROP: number of crop types\n')
        outFileExport.write(str(NCROP))
        outFileExport.write('\n# NFIELD: number of irrigation fields\n')
        outFileExport.write(str(NFIELD))
        outFileExport.write('\n# inputZON_TS_P_irr_fn: P zones\n')
        outFileExport.write(inputZON_TS_P_irr_fn)
        outFileExport.write('\n# inputZON_TS_TF_irr_fn: TF zones\n')
        outFileExport.write(inputZON_TS_TF_irr_fn)
        outFileExport.write('\n# inputZON_TS_PT_irr_fn: PT zones\n')
        outFileExport.write(inputZON_TS_PT_irr_fn)
        outFileExport.write('\n# Zr_c\n')
        for f in range(NCROP):
            outFileExport.write(str(Zr_c[f])+' ')
        outFileExport.write('\n# kTg_min_c\n')
        for f in range(NCROP):
            outFileExport.write(str(kTg_min_c[f])+' ')
        outFileExport.write('\n# kTg_max_c\n')
        for f in range(NCROP):
            outFileExport.write(str(kTg_max_c[f])+' ')
        outFileExport.write('\n# kT_f_c\n')
        for f in range(NCROP):
            outFileExport.write(str(kT_f_c[f])+' ')
        outFileExport.write('\n# kT_s_c\n')
        for f in range(NCROP):
            outFileExport.write(str(kT_s_c[f]) + ' ')
        outFileExport.write('\n# input_TS_crop_irr_fn: field crop time series\n')
        outFileExport.write(input_TS_crop_irr_fn)

    # EXPORTING DATE
    name = date_fn
    ws = pathMMws
    row1 = '#\n'
    TStmp = 0
    ExportResults(name, ws, row1, datenumOUT, J_d, TS = TStmp, TypeFile = "Date")

    return outMMsurf_fn
    ###############################################
     #   EOF   #
