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

import sys
import os
import numpy as np
import pylab
#USER INPUT: indicate in the line below the path where you copied
#            PET_PM_FAO56.py and startPET_PM_FAO56.py
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESsurf')
import PET_RF_INTER, plotPET, plotRF

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
    column 3:     RF: rainfall [mm]
    column 4:     Ta: air temperature measured at 2m [ºC]
    column 5:     RHa: relative air humidity [%] measured at heigth z_h [m]
    column 6:     Pa: air pressure [kPa]
    column 7:     u_z_m: windspeed measured at heigth z_m [m.s-1]
    column 8:     Rs: incoming solar (=shortwave) radiation [MJ.m-2.hour-1]
    Example 19 from FAO56 pag 75 (RF added by myself)
    - BOF - don't copy this line in a txt file
    date       time  RF Ta RHa Pa      u_z_m   Rs  # the first line is not read but is compulsory
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
    #NVEG: number of vegetation types
    3
    # to repeat NVEG times in a same line
    # 0 - VegType: name of the vegetation type [string, 1 word, no space allowed]
    # 1 - h: heigth of plant [m]
    # 2 - S_w: canopy capacity [mm]
    # 3 - C_leaf_star: maximum leaf conductance [m.s-1]
    # 4 - LAI_d: leaf area index dry season [m2.m-2]
    # 5 - LAI_w: leaf area index wet season [m2.m-2]
    # 6 - f_s_d: shelter factor dry season []
    # 7 - f_s_w: shelter factor wet season []
    # 8 - alfa_vd: vegetation albedo in dry season
    # 9 - alfa_vw: vegetation albedo in wet season
    # 10 - J_vd: starting julian day of the dry season [int 1-365]
    # 11 - J_vw: starting julian day of the wet season [int 1-365]
    # 12 - TRANS_vdw: transition period between dry and wet season [days]
    # 13 - Zr: maximum root depth [m]
    # 14 - k_Tu_d: dry season transpiration sourcing factor [], 1>=k_Tu>0
    # 15 - k_Tu_w: wet season transpiration sourcing factor [], 1>=k_Tu>0
    # BY DEFAULT THE PROGRAM WILL COMPUTE ETref USING GRASS FAO56 PARAMETERS
    # grassFAO56 0.12 0.1 0.01 2.88 2.88 0.5 0.5 0.23 0.23 150 270 20 0.25 1.0 1.0
    # VEG1 (grass muelledes)
    grassMU 0.5 0.15 0.01 1.2 2.88 0.5 0.55 0.30 0.20 150 270 20 0.40 1.0 1.0
    # VEG2 (Q.ilex Sardon)
    Qilex 6.0 1.16 0.004 6.0 6.0 0.75 0.70 0.25 0.10 150 270 20 6.0 0.2 0.9
    # VEG3 (Q.pyr sardon)
    Qpyr 8.0 1.75 0.004 6.0 1.0 0.80 0.3 0.25 0.15 150 270 20 4.0 0.5 0.9


    # SOIL PARAMETERS
    # NSOIL: number of soil types
    2
    # to repeat NSOIL times in a same line
    # 0 - SoilType: name of the soil type [string, 1 word, no space allowed]
    # 1 - Sy: surface (1cm) soil specific yield [m3.m-3]
    # 2 - alfa_sd: soil albedo in dry season
    # 3 - alfa_sw: soil albedo in wet season
    # 4 - J_sd: starting julian day of the dry season [int 1-365]
    # 5 - J_sw: starting julian day of the wet season [int 1-365]
    # 6 - TRANS_sdw: transition period between dry and wet season [days]
    # SOIL1 (alluvium)
    alluvium 0.05 0.20 0.15 150 270 20
    # SOIL2 (regolith)
    regolith 0.15 0.25 0.20 150 270 20

    # BY DEFAULT E0 (evaporation from open water using the Penman equation) WILL BE COMPUTED
    # OPEN WATER PARAMETERS
    # 1 - alfa_w: water albedo
    # 0.06

    # PLOTTING: positive integer for plotting, 0 or negative integer for NO plotting
    0
    - EOF - don't copy this line in a txt file

    Output:
    -------------------------------
    ASCII files (comma separated) with datetime, julian day,
    E0/ETref/PETveg/PETsoil and number of data/per day if daily data
'''

###########################################

def MMsurf(inputFOLDER_fn="", inputFile_TS_fn="", inputFile_PAR_fn="", outputFILE_fn="", pathMARM=""):


    def ExportResults(Dates, J, TS, outFileExport, ts_output = 0, n_d = [], TypeFile = "PET"):
        """
        Write the processed data in a open txt file
        INPUT:      output fluxes time series, date, Julian day, etc... and open file
        """
        for t in range(len(Dates)):
            year='%4d'%pylab.num2date(Dates[t]).year
            month='%02d'%pylab.num2date(Dates[t]).month
            day='%02d'%pylab.num2date(Dates[t]).day
            hour='%02d'%pylab.num2date(Dates[t]).hour
            minute='%02d'%pylab.num2date(Dates[t]).minute
            date_t=(year + "-" + month + "-" + day + " " + hour + ":" + minute)
            if TypeFile == "PET":
                if ts_output == 0:
                    out_line =  date_t, ', ', '%3d'% J[t], ",", '%14.9G' %TS[t],'\n'  #pylab.num2date(Dates[t]).isoformat()[:10]
                else:
                    out_line =  date_t, ', ', '%3d'% J[t], ",", '%14.9G' %TS[t], ",",'%2d'% n_d[t],'\n'  #pylab.num2date(Dates[t]).isoformat()[:10]
            elif TypeFile == "VAR":
                out_line =  date_t, ', ', '%3d'% J[t], ",", '%14.9G' %TS[0][t], ',', '%14.9G' %TS[1][t], ',', '%14.9G' %TS[2][t], ',', '%14.9G' %TS[3][t],',', '%14.9G' %TS[4][t],',', '%14.9G' %TS[5][t],',', '%14.9G' %TS[6][t],'\n'
            elif TypeFile == "RF":
                out_line =  date_t, ', ', '%3d'% J[t], ",", '%14.9G' %TS[0][t], ',', '%14.9G' %TS[1][t], ',', '%14.9G' %TS[2][t], ',', '%14.9G' %TS[3][t], ',', '%14.9G' %TS[4][t],'\n'
            elif TypeFile == "avRF_E0":
                out_line =  date_t, ', ', '%14.9G' %TS[0][t], ',', '%14.9G' %TS[1][t], ',', '\n'
            elif TypeFile == "Date":
                out_line =  date_t, '\n'
            for l in out_line:
                outFileExport.write(l)


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
    print '\n##############'
    print 'Reading MARMITESsurface INPUT FILES...'

    # METEO/VEGETATION/SOIL/WATER PARAMETERS
    inputFile = []
    inputFile_PAR_fn = os.path.join(inputFOLDER_fn, inputFile_PAR_fn)
    if os.path.exists(inputFile_PAR_fn):
        fin = open(inputFile_PAR_fn, 'r')
    else:
        print "\nThe input file [" + inputFile_PAR_fn + "] doesn't exist, verify name and path!"
        sys.exit()
    line = fin.readline().split()
    delimChar = line[0]
    try:
        for line in fin:
            line_tmp = line.split(delimChar)
            if not line_tmp == []:
                if (not line_tmp[0] == '') and (not line_tmp[0] == '\n'):
                    inputFile.append(line_tmp[0])
            else:
                raise NameError('InputFileFormat')
    except NameError:
        print "\nError in the input file, check format!\n[" + inputFile_PAR_fn +"]"
        sys.exit()
    except:
        print "\nUnexpected error in the input file\n[" + inputFile_PAR_fn +"]\n" + str(sys.exc_info()[0])
        sys.exit()
    try:
        # METEO PARAMETERS
        l=0
        line = inputFile[l].split()
        NMETEO = int(line[0])
        phi = []
        Lm = []
        Z = []
        Lz = []
        FC = []
        DTS = []
        z_m = []
        z_h = []
        if NMETEO>0:
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
        NVEG = int(line[0]) # number of vegetation types
        VegType = [] # name of the vegetation type [string]
        h = [] # heigth of plant [m]
        S = [] # canopy capacity [mm]
        C_leaf_star = [] # bulk stomatal resistance [s.m-1]
        LAI_d = [] # leaf area index dry season [m2.m-2]
        LAI_w = [] # leaf area index wet season [m2.m-2]
        f_s_d = [] # shelter factor dry season []
        f_s_w = [] # shelter factor wet season []
        alfa_vd = [] # vegetation albedo in dry season
        alfa_vw = [] # vegetation albedo in wet season
        J_vd = [] # starting julian day of the dry season [int 1-365]
        J_vw = [] # starting julian day of the wet season [int 1-365]
        TRANS_vdw = [] # transition period between dry and wet season [days]
        Zr = [] # maximum root depth [m]
        k_Tu_d = [] # dry season transpiration sourcing factor [], 1>=k_Tu>0
        k_Tu_w = [] # wet season transpiration sourcing factor [], 1>=k_Tu>0
        # input FAO56 grass parameters
        VegType.append("grassFAO56")
        h.append(0.12)
        S.append(0.1)
        C_leaf_star.append(0.01)
        LAI_d.append(2.88)
        LAI_w.append(2.88)
        f_s_d.append(0.5)
        f_s_w.append(0.5)
        alfa_vd.append(0.23)
        alfa_vw.append(0.23)
        J_vd.append(150)
        J_vw.append(270)
        TRANS_vdw.append(20)
        Zr.append(0.25)
        k_Tu_d.append(1.0)
        k_Tu_w.append(1.0)
        if NVEG>0:
            for i in range(NVEG):
                l = l + 1
                line = inputFile[l].split()
                VegType.append(str(line[0]))
                h.append(float(line[1]))
                S.append(float(line[2]))
                C_leaf_star.append(float(line[3]))
                LAI_d.append(float(line[4]))
                LAI_w.append(float(line[5]))
                f_s_d.append(float(line[6]))
                f_s_w.append(float(line[7]))
                alfa_vd.append(float(line[8]))
                alfa_vw.append(float(line[9]))
                J_vd.append(int(line[10]))
                J_vw.append(int(line[11]))
                TRANS_vdw.append(int(line[12]))
                Zr.append(float(line[13]))
                k_Tu_d.append(float(line[14]))
                k_Tu_w.append(float(line[15]))
        NVEG = NVEG + 1 # number of vegetation types + FAO56 GRASS (default)
        # SOIL PARAMETERS
        l = l + 1
        line = inputFile[l].split()
        NSOIL = int(line[0]) # number of soil types
        SoilType = [] # name of the soil type [string]
        Sy = [] # surface (1cm) soil specific yield [m3.m-3]
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
                Sy.append(float(line[1]))
                alfa_sd.append(float(line[2]))
                alfa_sw.append(float(line[3]))
                J_sd.append(int(line[4]))
                J_sw.append(int(line[5]))
                TRANS_sdw.append(int(line[6]))
        # WATER PARAMETERS
        alfa_w = 0.06 # water albedo
    except:
        print "\nError reading the input file:\n[" + inputFile_PAR_fn +"]"
        sys.exit()
    fin.close()
    del line, line_tmp, inputFile, l
    print "\nMETEO/VEGETATION/SOIL/OPENWATER PARAMETERS file imported!\n[" + inputFile_PAR_fn +"]"

    # read INPUT file
    # METEO TIME SERIES
    separator = ' ' # '\t' for tab delimited data file
    try:
        inputFile_TS_fn = os.path.join(inputFOLDER_fn, inputFile_TS_fn)
        if os.path.exists(inputFile_TS_fn):
            dataMETEOTS = np.loadtxt(inputFile_TS_fn, skiprows = 1, dtype = str)
            date = dataMETEOTS[:,0]
            time = dataMETEOTS[:,1]
            datenum = []
            datetime = []
            datenum_d = []
            datetime_i = (date[0] + " " + time[0])
            datenum_i = pylab.datestr2num(datetime_i)-(DTS[0])/24.0
            actual_day = pylab.num2date(datenum_i).isoformat()[:10]
            datenum_d.append(pylab.datestr2num(actual_day))
            RF1 = []
            Ta1 = []
            RHa1 = []
            Pa1 = []
            u_z_m1 = []
            Rs1 = []
            RF = []
            Ta = []
            RHa = []
            Pa = []
            u_z_m = []
            Rs = []
            if NMETEO>0:
                for n in range(NMETEO):
                    RF1.append(dataMETEOTS[:,2+n*6])
                    Ta1.append(dataMETEOTS[:,3+n*6])
                    RHa1.append(dataMETEOTS[:,4+n*6])
                    Pa1.append(dataMETEOTS[:,5+n*6])
                    u_z_m1.append(dataMETEOTS[:,6+n*6])
                    Rs1.append(dataMETEOTS[:,7+n*6])
                    RF.append([])
                    Ta.append([])
                    RHa.append([])
                    Pa.append([])
                    u_z_m.append([])
                    Rs.append([])
                    for t in range(len(Ta1[n])):
                        if n == 0:
                            datetime.append(date[t] + " " + time[t])
                            datenum.append(pylab.datestr2num(datetime[t])-(DTS[0])/24.0)
                            if actual_day <> date[t]:
                                datenum_d.append(pylab.datestr2num(date[t]))
                                actual_day = date[t]
                        RF[n].append(float(RF1[n][t]))
                        Ta[n].append(float(Ta1[n][t]))
                        RHa[n].append(float(RHa1[n][t]))
                        Pa[n].append(float(Pa1[n][t]))
                        u_z_m[n].append(float(u_z_m1[n][t]))
                        Rs[n].append(float(Rs1[n][t]))
        else:
            print "\nThe input file [" + inputFile_TS_fn + "] doesn't exist, verify name and path!"
            sys.exit()
    except:
            print "\nUnexpected error in the input file\n[" + inputFile_TS_fn +"]\n" + str(sys.exc_info()[1])
            sys.exit()
    del dataMETEOTS, date, time, datetime, datetime_i, actual_day, RF1, Ta1, RHa1, Pa1, u_z_m1,Rs1
    print "\nMETEO TIME SERIES file imported!\n[" + inputFile_TS_fn +"]"


    #  ##### COMPUTING PET and INTERCEPTION ##############################################
    print "\n###############"
    print "Computing PET..."

    if NMETEO>0:

        inputZON_TS_RF_fn = "inputZONRF.txt"
        inputZONRF_fn = os.path.join(pathMARM, inputZON_TS_RF_fn)
        inputZONRF = open(inputZONRF_fn, 'w')

        inputZON_TS_PET_fn = "inputZONPET.txt"
        inputZONPET_fn = os.path.join(pathMARM, inputZON_TS_PET_fn)
        inputZONPET = open(inputZONPET_fn, 'w')

        inputZON_TS_RFe_fn = "inputZONRFe.txt"
        inputZONRFe_fn = os.path.join(pathMARM, inputZON_TS_RFe_fn)
        inputZONRFe = open(inputZONRFe_fn, 'w')

        inputZON_TS_PE_fn = "inputZONPE.txt"
        inputZONPE_fn = os.path.join(pathMARM, inputZON_TS_PE_fn)
        inputZONPE = open(inputZONPE_fn, 'w')

        inputZON_TS_E0_fn = "inputZONE0.txt"
        inputZONE0_fn = os.path.join(pathMARM, inputZON_TS_E0_fn)
        inputZONE0 = open(inputZONE0_fn, 'w')

        for n in range(NMETEO):
            print "\n###############"
            print "Processing data of ZONE" + str(n+1) + "/%s" % str(NMETEO)
            J, J_d, outputVAR, PET_PM_VEG, Erf, PE_PM_SOIL, E0, PET_PM_VEG_d, PE_PM_SOIL_d, E0_d, RF_d, RFint, RF_duration, n1_d, RFe_d, I_d = \
                PET_RF_INTER.process                                                               (datenum, datenum_d, \
                    RF[n], Ta[n], RHa[n], Pa[n], u_z_m[n], Rs[n], \
                    phi[n], Lm[n], Z[n], Lz[n], FC[n], z_m[n], z_h[n], \
                    NVEG, VegType, h, S, C_leaf_star, LAI_d, LAI_w,\
                    f_s_d, f_s_w, alfa_vd, alfa_vw, J_vd, J_vw, TRANS_vdw,\
                    NSOIL, SoilType, Sy, alfa_sd, alfa_sw, J_sd, J_sw, TRANS_sdw, alfa_w)

            #  #####  PLOTTING ##############################################
            print "\n###############"
            print "Exporting plots..."
            npdatenum = np.array(datenum)
            # outputVAR = [DELTA,gama, Rs, Rs_corr, Rs0, Rnl, r]
            plot_exportVAR_fn = os.path.join(inputFOLDER_fn,  outputFILE_fn + "_ZON" + str(n+1)+'_VAR.png')
            plotPET.plotVAR(strTitle = 'Penman-Monteith variables', x = npdatenum \
                    ,y5 = PET_PM_VEG[0], y4 = E0 \
                    ,y2 = outputVAR[2], y3=outputVAR[3] \
                    ,y1 = outputVAR[4], y6 = outputVAR[5]\
                    ,lbl_y5 = 'ETref', lbl_y4 = 'E0' \
                    ,lbl_y2 = 'Rs_meas', lbl_y3 = 'Rs_corr' \
                    ,lbl_y1 = 'Rs0', lbl_y6 = 'Rnl'
                    ,plot_exportVAR_fn = plot_exportVAR_fn
                    )
            # PLOTTING PET
            plot_exportPET_fn = os.path.join(inputFOLDER_fn,  outputFILE_fn + "_ZON" + str(n+1)+'_PET.png')
            npdatenum_d = np.array(datenum_d)
            plotPET.plot(strTitle = 'PET veg&soil', x = npdatenum_d \
                ,y1 = PET_PM_VEG_d, y2 = E0_d, y3 = PE_PM_SOIL_d \
                ,lbl_y1 = VegType, lbl_y3 = SoilType
                ,plot_exportPET_fn = plot_exportPET_fn
                )
            # PLOTTING RF and INTERCEPTION
            plot_exportRF_fn = os.path.join(inputFOLDER_fn,  outputFILE_fn + "_ZON" + str(n+1)+'_RF.png')
            plotRF.plot(strTitle = 'Rainfall and interception', x = npdatenum_d \
                    ,y1 = RF_d, y2 = RFint, y3 = I_d, y4 = RFe_d
                    ,lbl_y1 = 'RF (mm/d)',  lbl_y2 = 'RFint (mm/h/d)' \
                    ,lbl_y3 = 'I (mm/d)',lbl_y4 = 'RFe (mm/d)', lbl_veg = VegType\
                    ,plot_exportRF_fn = plot_exportRF_fn
                    )

            #  #####  EXPORTING ASCII FILES ##############################################
            try:
                if os.path.exists(inputFOLDER_fn):
                    print "\n###############"
                    print "Exporting output to ASCII files..."
            # EXPORTING RESULTS PET/ET0/PE/E0/RF/RFe/INTER/
                    for i in range(2):
                        ts_output = i
                        if ts_output == 0:
                            print "\nHourly output..."                  # hourly output
                            datenumOUT = datenum
                            dORh_suffix = "_h"  # for file naming in export
                        else:                                # daily output
                            print "\nDaily ouput..."
                            datenumOUT = datenum_d
                            dORh_suffix = "_d"  # for file naming in export
                            ExportResults1(RF_d, inputZONRF)
                            ExportResults1(E0_d, inputZONE0)
                        outputfile_fn = outputFILE_fn + "_ZON" + str(n+1) + dORh_suffix
                        if ts_output <> 1:
                            # OUPUT VARIABLES FOR PM EQUATION
                            outpathname = os.path.join(inputFOLDER_fn, outputfile_fn + "_VAR.int")
                            outFileExport = open(outpathname, 'w')
                            outFileExport.write('Date,J,DELTA,gama,Rs_meas,Rs_corr,Rs0,Rnl,r\n')
                            ExportResults(datenumOUT, J, outputVAR, outFileExport, TypeFile = "VAR")
                            outFileExport.close()
                        #VEG
                        for v in range(NVEG):
                            if ts_output == 1:
                                # PET
                                outpathname = os.path.join(inputFOLDER_fn, outputfile_fn + "_PET_VEG" + str(v) + "_" + VegType[v] + ".out")
                                outFileExport = open(outpathname, 'w')
                                outFileExport.write('Date,J,PET,n_d\n')
                                ExportResults(datenumOUT, J_d, PET_PM_VEG_d[v],outFileExport, ts_output, n1_d)
                                # RF, RFe, etc
                                outpathname = os.path.join(inputFOLDER_fn, outputfile_fn + "_RF_VEG" + str(v) + "_" + VegType[v] + ".out")
                                outFileExport = open(outpathname, 'w')
                                outFileExport.write('Date,J,RF_mm,duration_day,RFint_mm_h,RFe_mm,Interception_mm\n')
                                TStmp = [RF_d, RF_duration, RFint, RFe_d[v], I_d[v]]
                                ExportResults(datenumOUT, J_d, TStmp, outFileExport, TypeFile = "RF")
                                outFileExport.close()
                                # PET and Rfe for MARMITES
                                if v>0: #to skip ETref (NVEG=0)
                                    ExportResults1(PET_PM_VEG_d[v], inputZONPET)
                                    ExportResults1(RFe_d[v], inputZONRFe)
                                print "\n[" + outpathname + "] exported!"
                            else:
                                outpathname = os.path.join(inputFOLDER_fn, outputfile_fn + "_PET_VEG" + str(v) + "_" + VegType[v] + ".out")
                                outFileExport = open(outpathname, 'w')
                                outFileExport.write('Date,J,PET\n')
                                ExportResults(datenumOUT, J, PET_PM_VEG[v],outFileExport, ts_output)
                                outpathname1 = os.path.join(inputFOLDER_fn, outputfile_fn + "_Erf_VEG" + str(v) + "_" + VegType[v] + ".int")
                                outFileExport1 = open(outpathname1, 'w')
                                outFileExport1.write('Date,J,Erf\n')
                                ExportResults(datenumOUT, J, Erf[v],outFileExport1, ts_output)
                                outFileExport1.close()
                            outFileExport.close()
                            print "\n[" + outpathname + "] exported!"
                        # SOIL
                        if NSOIL>0:
                            for s in range(NSOIL):
                                outpathname = os.path.join(inputFOLDER_fn, outputfile_fn + "_PE_SOIL" + str(s+1) + "_" + SoilType[s] + ".out")
                                outFileExport = open(outpathname, 'w')
                                if ts_output == 1:
                                    # PE
                                    outFileExport.write('Date,J,PE,n_d\n')
                                    ExportResults(datenumOUT, J_d, PE_PM_SOIL_d[s],outFileExport, ts_output,n1_d)
                                    # PE for MARMITES
                                    ExportResults1(PE_PM_SOIL_d[s], inputZONPE)
                                else:
                                    outFileExport.write('Date,J,PE\n')
                                    ExportResults(datenumOUT, J, PE_PM_SOIL[s],outFileExport, ts_output)
                                outFileExport.close()
                                print "\n[" + outpathname + "] exported!"
                        # WATER
                        outpathname = os.path.join(inputFOLDER_fn, outputfile_fn + "_E0.out")
                        outFileExport = open(outpathname, 'w')
                        if ts_output == 1:
                            outFileExport.write('Date,J,E0,n_d\n')
                            ExportResults(datenumOUT, J_d, E0_d,outFileExport, ts_output,n1_d)
                        else:
                            outFileExport.write('Date,J,E0\n')
                            ExportResults(datenumOUT, J, E0,outFileExport, ts_output)
                        outFileExport.close()
                        print "\n[" + outpathname + "] exported!"
                else:
                    print "\nNo ASCII file export required."
            except:
                print "\nError in output file, some output files were not exported."


    inputZONRF.close()
    inputZONRFe.close()
    inputZONPET.close()
    inputZONPE.close()
    inputZONE0.close()

    # EXPORTING ZONEVEGSOILfile
    MARMsurf_fn = "inputMARMsurf.txt"
    date_fn = 'inputDATE.txt'
    outpathname = os.path.join(pathMARM, MARMsurf_fn)
    outFileExport = open(outpathname, 'w')
    outFileExport.write('#\n# NMETEO: number of meteo zones\n')
    outFileExport.write(str(NMETEO))
    outFileExport.write('\n# NVEG: number of vegetation types\n')
    outFileExport.write(str(NVEG-1))
    outFileExport.write('\n# NSOIL: number of soil types\n')
    outFileExport.write(str(NSOIL))
    outFileExport.write('\n# inputDate_fn: file with the dates\n')
    outFileExport.write(date_fn)
    outFileExport.write('\n# inputZON_TS_RF_fn: RF zones\n')
    outFileExport.write(inputZON_TS_RF_fn)
    outFileExport.write('\n# inputZON_TS_PET_fn: RF zones\n')
    outFileExport.write(inputZON_TS_PET_fn)
    outFileExport.write('\n# inputZON_TS_RFe_fn: RF zones\n')
    outFileExport.write(inputZON_TS_RFe_fn)
    outFileExport.write('\n# inputZON_TS_PE_fn: RF zones\n')
    outFileExport.write(inputZON_TS_PE_fn)
    outFileExport.write('\n# inputZON_TS_E0_fn: RF zones\n')
    outFileExport.write(inputZON_TS_E0_fn)
    outFileExport.write('\n# TRANS_vdw\n')
    for v in range(1,len(TRANS_vdw)):
        outFileExport.write(str(TRANS_vdw[v])+' ')
    outFileExport.write('\n# Zr\n')
    for v in range(1,len(Zr)):
        outFileExport.write(str(Zr[v])+' ')
    outFileExport.write('\n# k_Tu_slp\n')
    for v in range(1,len(k_Tu_d)):
        outFileExport.write(str(k_Tu_d[v])+' ')
    outFileExport.write('\n# k_Tu_inter\n')
    for v in range(1,len(k_Tu_w)):
        outFileExport.write(str(k_Tu_w[v])+' ')
    outFileExport.write('\n# TRANS_sdw\n')
    for v in range(len(TRANS_sdw)):
        outFileExport.write(str(TRANS_sdw[v])+' ')

    # EXPORTING DATE
    outpathname = os.path.join(pathMARM, date_fn)
    outFileExport = open(outpathname, 'w')
    TStmp = 0
    ExportResults(datenumOUT, J_d, TStmp, outFileExport, TypeFile = "Date")
    outFileExport.close()

    return MARMsurf_fn
    ###############################################
     #   EOF   #
