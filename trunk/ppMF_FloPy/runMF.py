#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      alf
#
# Created:     02-05-2011
# Copyright:   (c) alf 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
""" Post processing of MODFLOW data after MARMITES run."""

import pylab
import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
sys.path.append(r'E:\00code\MARMITES\trunk\ppMF_FloPy')
import ppMODFLOW_flopy as ppMF
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESutilities')
import MARMITESprocess as MMproc
sys.path.append(r'E:\00code\MARMITES\trunk\MARMITESutilities\MARMITESplot')
import MARMITESplot as MMplot
import CreateColors
import h5py

MF_ws = r'E:\00code\00ws\zz_TESTS\MARMITESv2_r13c6l2\MF_ws'
MF_ini_fn = '_inputMF_flopy.ini'
MM_ws = r'E:\00code\00ws\zz_TESTS\MARMITESv2_r13c6l2'
inputDate_fn = 'inputDATE.txt'

MM_ws = r'E:\00code\00ws\zz_TESTS\MARMITESv2_r13c6l2'
MM_fn = '_inputMM.ini'

inputFile = MMproc.readFile(MM_ws,MM_fn)

l=0
try:
    # # ECHO ON/OFF (1 - interpreter verbose BUT NO report, 0 - NO interpreter verbose BUT report)
    verbose = int(inputFile[l].strip())
    l += 1
    # output plot (1 is YES, 0 is NO)
    plot_out  = int(inputFile[l].strip())
    l += 1
    #run MARMITESsurface  (1 is YES, 0 is NO)
    MMsurf_yn = int(inputFile[l].strip())
    l += 1
    # Define MARMITESsurface folder
    MMsurf_ws = inputFile[l].strip()
    l += 1
    # METEO TIME SERIES file name
    inputFile_TS_fn = inputFile[l].strip()
    l += 1
    # METEO/VEGETATION/SOIL/WATER PARAMETERS file name
    inputFile_PAR_fn = inputFile[l].strip()
    l += 1
    # ouputprefix
    outputFILE_fn = inputFile[l].strip()
    l += 1
    # ZONEVEGSOILfile
    outMMsurf_fn = inputFile[l].strip()
    l += 1
    # Define MODFLOW ws folders
    MF_ws = inputFile[l].strip()
    l += 1
    MF_ini_fn = inputFile[l].strip()
    l += 1
    plot_freq =  int(inputFile[l].strip())
    l += 1
    #GRID (ll means lower left)
    xllcorner = float(inputFile[l].strip())
    l += 1
    yllcorner = float(inputFile[l].strip())
    l += 1
    gridMETEO_fn = inputFile[l].strip()
    l += 1
    gridSOIL_fn = inputFile[l].strip()
    l += 1
    gridSOILthick_fn = inputFile[l].strip()
    l += 1
    gridIRR_fn = inputFile[l].strip()
    l += 1
    gridPONDhmax_fn =  inputFile[l].strip()
    l += 1
    gridPONDw_fn =  inputFile[l].strip()
    l += 1
    SOILparam_fn = inputFile[l].strip()
    l += 1
    IRR_fn = inputFile[l].strip()
    l += 1
    inputObs_fn = inputFile[l].strip()
    l += 1
    inputObsHEADS_fn = inputFile[l].strip()
    l += 1
    inputObsSM_fn = inputFile[l].strip()
    l += 1
    convcrit = float(inputFile[l].strip())
    l += 1
    ccnum = int(inputFile[l].strip())
    l += 1
    chunks = int(inputFile[l].strip())
except:
    print '\nType error in the input file %s' % (MM_fn)
    sys.exit()
del inputFile

MF_ws = os.path.join(MM_ws,MF_ws)
if os.path.exists(MM_ws):
    if os.path.exists(MF_ws):
        print ('\nMARMITES workspace:\n%s\n\nMODFLOW workspace:\n%s' % (MM_ws, MF_ws))
    else:
        print "The folder %s doesn't exist!\nPlease create it and run the model again." % (MF_ws)
else:
    print "The folder %s doesn't exist!\nPlease create it and run the model again." % (MM_ws)

print'\n##############'
print 'Importing MODFLOW configuration file'
nrow, ncol, delr, delc, nlay, nper, perlen, nstp, hnoflo, hdry, laytyp, lenuni, itmuni = ppMF.ppMFini(MF_ws, MF_ini_fn, out = 'MM')

# compute cbc conversion factor from volume to mm
if lenuni == 1:
    conv_fact = 304.8
elif lenuni == 2:
    conv_fact = 1000.0
elif lenuni == 3:
    conv_fact = 10.0
else:
    print 'FATAL ERROR!\nDefine the length unit in the MODFLOW ini file!\n (see USGS Open-File Report 00-92)'
    sys.exit()
    # TODO if lenuni<>2 apply conversion factor to delr, delc, etc...
if laytyp[0]==0:
    print 'FATAL ERROR!\nThe first layer cannot be confined type!\nChange your parameter laytyp in the MODFLOW lpf package.\n(see USGS Open-File Report 00-92)'
    sys.exit()
if itmuni <> 4:
    print 'FATAL ERROR! Time unit is not in days!'
    sys.exit()

print'\n##############'
print 'Reading MARMITESsurf RUN'

try:
    inputFile = MMproc.readFile(MM_ws,outMMsurf_fn)
except:
    print 'The file %s is not correct, run MARMITESsurface first!' % (outMMsurf_fn)

l=0
TRANS_vdw = []
Zr = []
k_Tu_slp = []
k_Tu_inter = []
TRANS_sdw = []
try:
    NMETEO = int(inputFile[l].strip())
    l += 1
    NVEG = int(inputFile[l].strip())
    l += 1
    NSOIL = int(inputFile[l].strip())
    l += 1
    inputDate_fn = str(inputFile[l].strip())
    l += 1
    inputZON_TS_RF_fn = str(inputFile[l].strip())
    l += 1
    inputZON_TS_PET_fn = str(inputFile[l].strip())
    l += 1
    inputZON_TS_RFe_fn = str(inputFile[l].strip())
    l += 1
    inputZON_TS_PE_fn = str(inputFile[l].strip())
    l += 1
    inputZON_TS_E0_fn = str(inputFile[l].strip())
    l += 1
    line = inputFile[l].split()
    for v in range(NVEG):
        TRANS_vdw.append(int(line[v]))
    l += 1
    line = inputFile[l].split()
    for v in range(NVEG):
        Zr.append(float(line[v]))
    l += 1
    line = inputFile[l].split()
    for v in range(NVEG):
        k_Tu_slp.append(float(line[v]))
    l += 1
    line = inputFile[l].split()
    for v in range(NVEG):
        k_Tu_inter.append(float(line[v]))
    l += 1
    line = inputFile[l].split()
    for v in range(NSOIL):
        TRANS_sdw.append(int(line[v]))
except:
    print '\nType error in file [' + inputFile_fn + ']'
    sys.exit()
del inputFile

if isinstance(nper, str):
    try:
        MFtime_fn = nper.split()[0].strip()
        perlenmax = int(nper.split()[1].strip())
    except:
        print '\nError in nper format of the MODFLOW ini file!\n'
        sys.exit()
    nper, perlen, nstp, inputZON_TS_RF_fn, inputZON_TS_PET_fn, inputZON_TS_RFe_fn, inputZON_TS_PE_fn, inputZON_TS_E0_fn = ppMF.ppMFtime(MM_ws, MF_ws, MFtime_fn, perlenmax, inputDate_fn, inputZON_TS_RF_fn, inputZON_TS_PET_fn, inputZON_TS_RFe_fn, inputZON_TS_PE_fn, inputZON_TS_E0_fn, NMETEO, NVEG, NSOIL)
else:
    MFtime_fn = None

print'\n##############'
print 'Running MODFLOW'
top_array, ibound, h5_MF_fn = ppMF.ppMF(MF_ws, MF_ini_fn, rch_input = 0.0001, rch_dft = 0.0001)

h5_MF = h5py.File(h5_MF_fn)

cbc_nam = []
for c in h5_MF['cbc_nam']:
    cbc_nam.append(c.strip())
iDRN = cbc_nam.index('DRAINS')
iSTO = cbc_nam.index('STORAGE')
iRCH = cbc_nam.index('RECHARGE')
iWEL = cbc_nam.index('WELLS')

h_MF = np.zeros((sum(perlen), nrow, ncol, nlay))
cbc_MF = np.zeros((sum(perlen), len(cbc_nam), nrow, ncol, nlay))
t = 0
ncell = 0
for n in range(nper):
    if perlen[n]>1:
        for x in range(perlen[n]):
            h_MF[t,:,:,:] = h5_MF['heads'][n,:,:,:]
            for i in range(nrow):
                for j in range(ncol):
                    cbc_MF[t,:,i,j,:] = conv_fact*h5_MF['cbc'][n,:,i,j,:]/(delr[j]*delc[i])
                    if ibound[i,j,0]<>0.0:
                        ncell += 1
            t += 1
    else:
        h_MF[t,:,:,:] = h5_MF['heads'][n,:,:,:]
        for i in range(nrow):
            for j in range(ncol):
                cbc_MF[t,:,i,j,:] = conv_fact*h5_MF['cbc'][n,:,i,j,:]/(delr[j]*delc[i])
                if ibound[i,j,0]<>0.0:
                        ncell += 1
        t += 1
h_MF_m = np.ma.masked_values(h_MF, hnoflo, atol = 0.09)
del h_MF
top_array_m = np.ma.masked_values(top_array, hnoflo, atol = 0.09)
h5_MF.close()

# #############################
# ###  OUTPUT EXPORT   ########
# #############################

print '\n##############\nMARMITES exporting...'

MM_PROCESS = MMproc.PROCESS(MM_ws               = MM_ws,
                            MF_ws                    = MF_ws,
                            nrow                     = nrow,
                            ncol                     = ncol,
                            xllcorner                = xllcorner,
                            yllcorner                = yllcorner,
                            cellsizeMF               = delr[0],
                            nstp                     = nstp,
                            hnoflo                   = hnoflo
                            )

# SOIL PARAMETERS
_nsl, _nam_soil, _st, _slprop, _Sm, _Sfc, _Sr, _Si, _Ks = MM_PROCESS.inputSoilParam(SOILparam_fn = SOILparam_fn, NSOIL = NSOIL)
_nslmax = max(_nsl)

# READ date of input files
inputDate_fn=os.path.join(MM_ws, inputDate_fn)
if os.path.exists(inputDate_fn):
    inputDate=np.loadtxt(inputDate_fn, dtype = str)
    inputDate = inputDate[:,0]
    inputDate = pylab.datestr2num(inputDate)
    for i in range(1,len(inputDate)):
        #__________________Check date consistency________________#
        difDay=inputDate[i]-inputDate[i-1]
        if (difDay !=1.0):
            print 'DifDay = ' + str(difDay)
            raise ValueError, 'The dates of the input data (RF and PET) are not sequencial, check your daily time step!\nError in date %s ' % str(inputDate[i])
else:
    raise ValueError, "\nThe file %s doesn't exist!!!" % inputDate_fn

# READ observations time series (heads and soil moisture)
obsCHECK = 1  # 0: no obs, 1: obs
if obsCHECK==1:
    print "\nReading observations time series (hydraulic heads and soil moisture)..."
    obs, outpathname, obs_h, obs_S = MM_PROCESS.inputObs(
                            inputObs_fn = inputObs_fn,
                            inputObsHEADS_fn = inputObsHEADS_fn,
                            inputObsSM_fn = inputObsSM_fn,
                            inputDate   = inputDate,
                            _nslmax    = _nslmax,
                            MFtime_fn  = MFtime_fn
                            )
    # Write first output in a txt file
    outPESTheads_fn      = 'PESTheads.dat'
    outPESTsm_fn         = 'PESTsm.dat'
    outPESTheads=open(os.path.join(MM_ws,outPESTheads_fn), 'w')
    outPESTsm=open(os.path.join(MM_ws,outPESTsm_fn), 'w')
else:
    print "\nNo reading of observations time series (hydraulic heads and soil moisture) required."

print '\nExporting ASCII files and plots...'
h_MF_m = np.ma.masked_values(h_MF_m, hdry, atol = 1E+25)
hmax = []
hmin = []
DRNmax = []
DRNmin = []
cbcmax = []
cbcmin = []
for L in range(nlay):
    hmax.append(np.nanmax(h_MF_m[:,:,:,:].flatten()))
    hmin.append(np.nanmin(h_MF_m[:,:,:,:].flatten()))
    DRNmax.append(np.nanmax(-cbc_MF[:,iDRN,:,:,:]).flatten())
    DRNmin.append(np.nanmin(-cbc_MF[:,iDRN,:,:,:]).flatten())
    cbcmax.append(np.nanmax(-cbc_MF[:,:,:,:,:]).flatten())
    cbcmin.append(np.nanmin(-cbc_MF[:,:,:,:,:]).flatten())
for o in range(len(obs.keys())):
        npa_m_tmp = np.ma.masked_values(obs_h[o], hnoflo, atol = 0.09)
        hmax.append(np.nanmax(npa_m_tmp.flatten()))
        hmin.append(np.nanmin(npa_m_tmp.flatten()))
hmax = float(np.ceil(np.nanmax(hmax)))
hmin = float(np.floor(np.nanmin(hmin)))
DRNmax = float(np.ceil(np.nanmax(DRNmax)))
DRNmin = float(np.floor(np.nanmin(DRNmin)))
cbcmax = float(np.ceil(np.nanmax(cbcmax)))
cbcmin = float(np.floor(np.nanmin(cbcmin)))
# plot water budget
plt_export_fn = os.path.join(MM_ws, '_plt_UNSATandGWbudgets.png')
index_cbc = [iRCH, iSTO, iDRN, iWEL]
flxlbl = []
flxlst = []
for l in range(nlay):
    for x in range(len(index_cbc)):
        flxlst.append(cbc_MF[:,index_cbc[x],:,:,l].sum()/sum(perlen)/ncell)
        flxlbl.append('GW_' + cbc_nam[index_cbc[x]] + '_L' + str(l+1))
colors_flx = CreateColors.main(hi=0, hf=180, numbcolors = len(flxlbl))
MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = 'MARMITES and MODFLOW water flux budget for the whole catchment', fluxmax = cbcmax, fluxmin = cbcmin)
del flxlst

# exporting ASCII files and plots at observations cells
if obsCHECK == 1:
    colors_nsl = CreateColors.main(hi=00, hf=180, numbcolors = (_nslmax+1))
    for o in range(len(obs.keys())):
        i = obs.get(obs.keys()[o])['i']
        j = obs.get(obs.keys()[o])['j']
        # export ASCII file at piezometers location
        #TODO extract heads at piezo location and not center of cell
        MM_PROCESS.ExportResultsPEST(i, j, inputDate, _nslmax, h_MF_m[:,i,j,0], obs_h[o][0,:], obs_S[o], outPESTheads, outPESTsm, obs.keys()[o])
        # plot water budget at each obs. cell
        plt_title = obs.keys()[o]
        plt_export_fn = os.path.join(MM_ws, '_plt_'+ obs.keys()[o] + 'UNSATandGWbudgets.png')
        flxlst =[]
        for l in range(nlay):
            for x in range(len(index_cbc)):
                flxlst.append(cbc_MF[:,index_cbc[x],i,j,l].sum()/sum(perlen))
        MMplot.plotGWbudget(flxlst = flxlst, flxlbl = flxlbl, colors_flx = colors_flx, plt_export_fn = plt_export_fn, plt_title = plt_title, fluxmax = cbcmax, fluxmin = cbcmin)
        del flxlst
    # output for PEST
    outPESTheads.close()
    outPESTsm.close()
    del obs, obs_h, obs_S

del MM_PROCESS

if plot_out == 1:
    # plot heads (grid + contours), DRN, etc... at specified TS
    TSlst = []
    TS = 0
    while TS < len(h_MF_m):
        TSlst.append(TS)
        TS += plot_freq
    TSlst.append(len(h_MF_m)-1)
    for TS in TSlst:
        # plot heads [m]
        V=[]
        for L in range(nlay):
            V.append(h_MF_m[TS,:,:,L])
        MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'hydraulic heads elevation (m)', msg = 'DRY', plt_title = 'HEADS', MM_ws = MM_ws, interval_type = 'arange', interval_diff = 0.5, Vmax = hmax, Vmin = hmin)
        del V
        # plot diff between drain elevation and heads elevation [m]
        DrnHeadsLtop = top_array_m - h_MF_m[TS,:,:,0]
        DrnHeadsLtop_m = np.ma.masked_greater(DrnHeadsLtop,0.0)
        V = [DrnHeadsLtop_m]
        diffMin = 0
        diffMax = np.nanmin(DrnHeadsLtop_m)
        MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = 1, V = V,  cmap = plt.cm.RdYlGn, CBlabel = 'diff. between DRN elev and hyd. heads elev. (m)', msg = ' - no drainage', plt_title = 'HEADSDRNdiff', MM_ws = MM_ws, interval_type = 'linspace', interval_num = 5, Vmax = diffMin, Vmin = diffMax, fmt='%.3G')
        # plot GW drainage [mm]
        V = []
        for L in range(nlay):
            V.append(-cbc_MF[TS,iDRN,:,:,L])
        MMplot.plotLAYER(TS, ncol = ncol, nrow = nrow, nlay = nlay, nplot = nlay, V = V,  cmap = plt.cm.Blues, CBlabel = 'groundwater drainage (mm/time step)', msg = '- no drainage', plt_title = 'DRN', MM_ws = MM_ws, interval_type = 'linspace', interval_num = 5, Vmin = DRNmin, Vmax = DRNmax, fmt='%.3G')
        del DrnHeadsLtop, DrnHeadsLtop_m, V

del cbc_MF, h_MF_m, MM, MM_S
del top_array_m, gridSOIL, inputDate
del hmax, hmin, DRNmax, DRNmin, cbcmax, cbcmin

del nrow, ncol, delr, delc, nlay, nstp, nper, hnoflo, hdry, ibound, laytyp, h5_MF_fn, top_array, lenuni



print "\nMODFLOW done!"

