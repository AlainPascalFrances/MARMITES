# # TEMP SCRIPT to process RF data from RF logger INST to HOUR # #
# compute hourly RF based on instantaneous measurement of RF

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import MARMITESprocess_v3 as MMproc

#####################

def ExportResults(Dates, TS, outFileExport):
    """
    Write the processed data in a open txt file
    INPUT:      output fluxes time series, date, Julian day, etc... and open file
    """
    for t in range(len(Dates)):
        year='%4d'%mpl.dates.num2date(Dates[t]).year
        month='%02d'%mpl.dates.num2date(Dates[t]).month
        day='%02d'%mpl.dates.num2date(Dates[t]).day
        hour='%02d'%mpl.dates.num2date(Dates[t]).hour
        minute='%02d'%mpl.dates.num2date(Dates[t]).minute
        date_t=(year + "-" + month + "-" + day + " " + hour + ":" + minute)
        out_line =  date_t, ' ', '%14.9G' %TS[0][t],'\n'
        for l in out_line:
            outFileExport.write(l)

#####################

inputFOLDER_fn = r'E:\00code_ws\SARDON\MMsurf_ws'
inputFILE_fn = 'RF_trab.txt'
outputFILE_fn = 'RF_trab.dRF'

inputFILE_fn = os.path.join(inputFOLDER_fn, inputFILE_fn)
outputFILE_fn = os.path.join(inputFOLDER_fn, outputFILE_fn)

if os.path.exists(inputFILE_fn):
    dataMETEOTS = np.loadtxt(inputFILE_fn, skiprows = 1, dtype = str)
    date = dataMETEOTS[:,0]
    time = dataMETEOTS[:,1]
    datenum = []
    datetime = []
    datenum_d = []
    datetime_i = (date[0] + " " + time[0])
    datenum_i = mpl.dates.datestr2num(datetime_i)
    actual_day = mpl.dates.num2date(datenum_i).isoformat()[:10]
    datenum_d.append(mpl.dates.datestr2num(actual_day))
    RF  = []
    RF1 = dataMETEOTS[:,2]
    for t in range(len(RF1)):
        datetime.append(date[t] + " " + time[t])
        datenum.append(mpl.dates.datestr2num(datetime[t]))
        if actual_day <> date[t]:
            datenum_d.append(mpl.dates.datestr2num(date[t]))
            actual_day = date[t]
        RF.append(float(RF1[t]))

RF_h =  []
datenumRF_h = []
RF_h_tmp = 0.0
actual_dayhour = mpl.dates.num2date(datenum[0]).isoformat()[:13]
for t in range(len(datenum)):
    if actual_dayhour == mpl.dates.num2date(datenum[t]).isoformat()[:13]:
        RF_h_tmp = RF_h_tmp + RF[t]
    else:
        RF_h.append(RF_h_tmp)
        datenumRF_h.append(mpl.dates.datestr2num(actual_dayhour)+1.0/24)
        RF_h_tmp = RF[t]
        actual_dayhour = mpl.dates.num2date(datenum[t]).isoformat()[:13]
RF_h.append(RF_h_tmp)
datenumRF_h.append(mpl.dates.datestr2num(actual_dayhour)+1.0/24)
# export E0 during RF > RFtreshold with datehour (intermediate file)
outpathname = os.path.join(inputFOLDER_fn, outputFILE_fn)
outFileExport = open(outpathname, 'w')
outFileExport.write('Date Time RF_mm\n')
TStmp = [RF_h]
ExportResults(datenumRF_h, TStmp, outFileExport)
outFileExport.close()
print "\n[" + outpathname + "] exported!"