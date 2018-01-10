# # TEMP SCRIPT to process RF data from RF logger INST to HOUR # #
# compute hourly RF based on instantaneous measurement of RF

import os
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
    RF  = []
    Ta  = []
    RH  = []
    WS  = []
    Sin = []
    RF1 = (dataMETEOTS[:,0])
    Ta1 = (dataMETEOTS[:,1])
    RH1 = (dataMETEOTS[:,2])
    WS1  = (dataMETEOTS[:,4])
    Sin1 = (dataMETEOTS[:,5])
    for t in range(len(RF1)):
        RF.append(float(RF1[t]))
        Ta.append(float(Ta1[t]))
        RH.append(float(RH1[t]))
        WS.append(float(WS1[t]))
        Sin.append(float(Sin1[t]))

x = np.zeros((len(RF),5), dtype = np.float)
x[:,0] = Ta
x[:,1] = RH
x[:,2] = WS
x[:,3] = Sin
x[:,4] = RF
mymodel = ols.ols(x[:,0],x[:,1:4],y_varnm = 'Ta',x_varnm = ['RH','WS','Sin'])
mymodel.p               # return coefficient p-values
mymodel.summary()       # print results
print
mymodel = ols.ols(x[:,1],x[:,2:5],y_varnm = 'RH',x_varnm = ['WS','Sin','RF'])
mymodel.p               # return coefficient p-values
mymodel.summary()       # print results

print 'DONE!'