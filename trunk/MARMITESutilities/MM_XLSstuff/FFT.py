# # TEMP SCRIPT to process RF data from RF logger INST to HOUR # #
# compute hourly RF based on instantaneous measurement of RF

import os
import scikits.timeseries.lib.interpolate as ts_inter
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import MARMITESprocess_v3 as MMproc
import ols

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
inputFILE_fn = '__meteoTB1.txt'
outputFILE_fn = '__meteoTB_filled.txt'

inputFILE_fn = os.path.join(inputFOLDER_fn, inputFILE_fn)
outputFILE_fn = os.path.join(inputFOLDER_fn, outputFILE_fn)

if os.path.exists(inputFILE_fn):
    dataMETEOTS = np.loadtxt(inputFILE_fn, skiprows = 1, dtype = str, )
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
    Ta  = []
    RH  = []
    WS  = []
    Sin = []
    RF1 = dataMETEOTS[:,2]
    Ta1 = dataMETEOTS[:,3]
    RH1 = dataMETEOTS[:,4]
    WS1  = dataMETEOTS[:,6]
    Sin1 = dataMETEOTS[:,7]
    for t in range(len(RF1)):
        datetime.append(date[t] + " " + time[t])
        datenum.append(mpl.dates.datestr2num(datetime[t]))
        if actual_day <> date[t]:
            datenum_d.append(mpl.dates.datestr2num(date[t]))
            actual_day = date[t]
        RF.append(float(RF1[t]))
        Ta.append(float(Ta1[t]))
        RH.append(float(RH1[t]))
        WS.append(float(WS1[t]))
        Sin.append(float(Sin1[t]))

Ta_m = np.ma.masked_values(Ta, -999.99, atol = 0.09)
#X = ts_inter.interp_masked1d(Ta_m, 'cubic')  #{‘constant’, ‘linear’, ‘cubic’, quintic’}
#plt.plot(datenum,X)
#plt.plot(datenum,Ta_m, 'o')
#plt.show()

x = np.zeros((len(datenum),4), dtype = np.float)
x[:,0] = Ta_m
x[:,1] = RH
x[:,2] = WS
x[:,3] = Sin
mymodel = ols.ols(x[:,0],x[:,2:],y_varnm = 'Ta',x_varnm = ['WS','Sin'])
mymodel.p               # return coefficient p-values
mymodel.summary()       # print results
print
mymodel = ols.ols(x[:,1],x[:,2:],y_varnm = 'RH',x_varnm = ['WS','Sin'])
mymodel.p               # return coefficient p-values
mymodel.summary()       # print results

##X = np.fft.fft(Ta_m)
##Y = np.zeros(len(Ta))
###Y[important frequencies] = X[important frequencies]
##plt.plot(datenum,X)
##plt.show()
##
##X1 = np.fft.fft(WS)
##Y1 = np.zeros(len(WS))
##plt.plot(datenum,X1)
##plt.show()
##
##X2 = np.fft.fft(Sin)
##Y2 = np.zeros(len(Sin))
##plt.plot(datenum,X2)
##plt.show()

##
##
##
##
##
##
##
##
##RF_h =  []
##datenumRF_h = []
##RF_h_tmp = 0.0
##actual_dayhour = mpl.dates.num2date(datenum[0]).isoformat()[:13]
##for t in range(len(datenum)):
##    if actual_dayhour == mpl.dates.num2date(datenum[t]).isoformat()[:13]:
##        RF_h_tmp = RF_h_tmp + RF[t]
##    else:
##        RF_h.append(RF_h_tmp)
##        datenumRF_h.append(mpl.dates.datestr2num(actual_dayhour)+1.0/24)
##        RF_h_tmp = RF[t]
##        actual_dayhour = mpl.dates.num2date(datenum[t]).isoformat()[:13]
##RF_h.append(RF_h_tmp)
##datenumRF_h.append(mpl.dates.datestr2num(actual_dayhour)+1.0/24)
### export E0 during RF > RFtreshold with datehour (intermediate file)
outpathname = os.path.join(inputFOLDER_fn, outputFILE_fn)
##outFileExport = open(outpathname, 'w')
##outFileExport.write('Date Time RF_mm\n')
##TStmp = [RF_h]
##ExportResults(datenumRF_h, TStmp, outFileExport)
##outFileExport.close()
print "\n[" + outpathname + "] exported!"