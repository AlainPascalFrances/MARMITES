# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:        startMARMITES
# Purpose:
#
# Author:      frances08512
#
# Created:     01-03-2011
# Copyright:   (c) frances08512 2010
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__author__ = "Alain P. Francés <frances08512@itc.nl>"
__version__ = "0.1"
__date__ = "March 2011"

import pylab
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
from mf import *
import mfreadbinaries as mfrdbin
#####################################

def convASCIIraster2array(filenameIN, arrayOUT, cellsizeMF, nrow, ncol):
    '''
    Read ESRI/MODFLOW ASCII file and convert to numpy array
    '''

    # Load the grid files
    if os.path.exists(filenameIN):
        fin = open(filenameIN, 'r')
    else:
        raise ValueError, "The file %s doesn't exist!!!" % filenameIN
#        fout = open(outfile, 'w')

    # Read the header
    line = fin.readline().split()
    ncol = int(line[1])

    line = fin.readline().split()
    nrow = int(line[1])

    line = fin.readline().split()
    xllcorner = float(line[1])

    line = fin.readline().split()
    yllcorner = float(line[1])

    line = fin.readline().split()
    cellsizeEsriAscii = float(line[1])

    line = fin.readline().split()
    NODATA_value = float(line[1])

    # Process the file
#    print "\nConverting %s to np.array..." % (filenameIN)
    for row in range(nrow):
        for col in range(ncol):
            if col == 0: line = fin.readline().split()
            arrayOUT[row,col]=line[col]

    # verify grid consistency between MODFLOW and ESRI ASCII
    if arrayOUT.shape[0] != nrow or arrayOUT.shape[1] != ncol or cellsizeMF != cellsizeEsriAscii:
        raise BaseException, '\nError in consistency between the MODFLOW grid and the input gridof the file %s.\nCheck the cell size and the number of rows, columns and cellsize' % filenameIN

    return arrayOUT

    fin.close()
#####################################

def ppMF(model_ws = '', MM_ws = ''):

    messagemanual="Please read the manual!\n(that by the way still doesn't exist...)"

    # read input file (called _input.ini in the MARMITES workspace
    # the first character on the first line has to be the character used to comment
    # the file can contain any comments as the user wish, but the sequence of the input has to be respected
    inputFileMF_fn = os.path.join(model_ws, "_inputMF_flopy.ini")
    inputFile = []
    if os.path.exists(inputFileMF_fn):
        fin = open(inputFileMF_fn, 'r')
    else:
        raise ValueError, "Input file doesn't exist, verify name and path!"
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
        raise ValueError, 'Error in the input file, check format!\n%s' % (messagemanual)
        sys.exit()
    except e:
        print "Unexpected error in the input file:\n", sys.exc_info()[0]
        print messagemanual
        sys.exit()
    l=0
    try:
        SP_d = 0 # daily stress period
        modelname =  str(inputFile[l].strip())
        l = l + 1
        namefile_ext = str(inputFile[l].strip())
        l = l + 1
        exe_name = str(inputFile[l].strip())
        l = l + 1
        dum_sssp1 = int(inputFile[l].strip())
        l = l + 1
        ext_dis = str(inputFile[l].strip())
        l = l + 1
        nlay = int(inputFile[l].strip())
        l = l + 1
        ncol = int(inputFile[l].strip())
        l = l + 1
        nrow = int(inputFile[l].strip())
        l = l + 1
        nper = int(inputFile[l].strip())
        if nper <0:
            nper = abs(nper)
            SP_d = 1
        l = l + 1
        itmuni = int(inputFile[l].strip())
        l = l + 1
        lenuni = int(inputFile[l].strip())
        l = l + 1
        laycbd = []
        laycbd_tmp =  inputFile[l].split()
        for i in range(nlay):
            laycbd.append(int(laycbd_tmp[i]))
        l = l + 1
        delr = []
        delr_tmp =  inputFile[l].split()
        for i in range(ncol):
            delr.append(int(delr_tmp[i]))
        l = l + 1
        delc = []
        delc_tmp =  inputFile[l].split()
        for i in range(nrow):
            delc.append(int(delc_tmp[i]))
        l = l + 1
        top_fn = str(inputFile[l].strip())
        l = l + 1
        botm_fn = []
        botm_tmp =  inputFile[l].split()
        for i in range(nlay):
            botm_fn.append(str(botm_tmp[i]))
        l = l + 1
        if SP_d == 0:
            perlen = []
            perlen_tmp =  inputFile[l].split()
            for i in range(nper):
                perlen.append(int(perlen_tmp[i]))
            l = l + 1
            nstp = []
            nstp_tmp =  inputFile[l].split()
            for i in range(nper):
                nstp.append(int(nstp_tmp[i]))
            l = l + 1
            tsmult = []
            tsmult_tmp =  inputFile[l].split()
            for i in range(nper):
                tsmult.append(int(tsmult_tmp[i]))
            l = l + 1
            Ss_tr = []
            Ss_tr_tmp =  inputFile[l].split()
            for i in range(nper):
                Ss_tr.append(str(Ss_tr_tmp[i]))
            l = l + 1
        else:
            perlen = []
            nstp = []
            tsmult = []
            Ss_tr = []
            for t in range(nper):
                perlen.append(1)
                nstp.append(1)
                tsmult.append(1)
                Ss_tr.append('TR')
            l = l + 4
        ext_bas = str(inputFile[l].strip())
        l = l + 1
        ibound_fn = []
        ibound_tmp =  inputFile[l].split()
        for i in range(nlay):
            ibound_fn.append(str(ibound_tmp[i]))
        l = l + 1
        strt_fn = []
        strt_tmp =  inputFile[l].split()
        for i in range(nlay):
            strt_fn.append(str(strt_tmp[i]))
        l = l + 1
        hnoflo = float(inputFile[l].strip())
        l = l + 1
        ext_lpf = str(inputFile[l].strip())
        l = l + 1
        ilpfcb = int(inputFile[l].strip())
        l = l + 1
        hdry = eval(inputFile[l].strip())
        l = l + 1
        nplpf = int(inputFile[l].strip())
        l = l + 1
        laytyp = []
        laytyp_tmp =  inputFile[l].split()
        for i in range(nlay):
            laytyp.append(int(laytyp_tmp[i]))
        l = l + 1
        layavg = []
        layavg_tmp =  inputFile[l].split()
        for i in range(nlay):
            layavg.append(int(layavg_tmp[i]))
        l = l + 1
        chani = []
        chani_tmp =  inputFile[l].split()
        for i in range(nlay):
            chani.append(int(chani_tmp[i]))
        l = l + 1
        layvka = []
        layvka_tmp =  inputFile[l].split()
        for i in range(nlay):
            layvka.append(int(layvka_tmp[i]))
        l = l + 1
        laywet = []
        laywet_tmp =  inputFile[l].split()
        for i in range(nlay):
            laywet.append(int(laywet_tmp[i]))
        l = l + 1
        hk_fn = []
        hk_tmp =  inputFile[l].split()
        for i in range(nlay):
            hk_fn.append(str(hk_tmp[i]))
        l = l + 1
        vka_fn = []
        vka_tmp =  inputFile[l].split()
        for i in range(nlay):
            vka_fn.append(str(vka_tmp[i]))
        l = l + 1
        ss_fn = []
        ss_tmp =  inputFile[l].split()
        for i in range(nlay):
            ss_fn.append(str(ss_tmp[i]))
        l = l + 1
        sy_fn = []
        sy_tmp =  inputFile[l].split()
        for i in range(nlay):
            sy_fn.append(str(sy_tmp[i]))
        l = l + 1
        ext_oc = str(inputFile[l].strip())
        l = l + 1
        ihedfm = int(inputFile[l].strip())
        l = l + 1
        iddnfm = int(inputFile[l].strip())
        l = l + 1
        ext_cbc = str(inputFile[l].strip())
        l = l + 1
        ext_heads = str(inputFile[l].strip())
        l = l + 1
        ext_ddn = str(inputFile[l].strip())
        l = l + 1
        ext_rch = str(inputFile[l].strip())
        l = l + 1
        nrchop = int(inputFile[l].strip())
        l = l + 1
        rch_tmp =  inputFile[l].split()
        try:
            rch_fn=float(rch_tmp[0])
        except ValueError:
            rch_fn = str(rch_tmp[0])
    except:
        print "Unexpected error in the input file:\n", sys.exc_info()[0]
        sys.exit()
    fin.close()

    # 1 - reaf asc file and convert in np.array

    print "\nImporting ESRI ASCII files..."

    top_path = os.path.join(model_ws, top_fn)
    top_array = np.zeros((nrow,ncol))
    top_array = convASCIIraster2array(top_path, top_array, cellsizeMF = delr[0], nrow = nrow, ncol=ncol)
    top_array = list(top_array)

    botm_array = np.zeros((nrow,ncol, nlay))
    for l in range(nlay):
        botm_path = os.path.join(model_ws, botm_fn[l])
        botm_array[:,:,l] = convASCIIraster2array(botm_path, botm_array[:,:,l], cellsizeMF = delr[0], nrow = nrow, ncol=ncol)
    botm_array = list(botm_array)

    ibound_array = np.zeros((nrow,ncol, nlay), dtype = int)
    for l in range(nlay):
        ibound_path = os.path.join(model_ws, ibound_fn[l])
        ibound_array[:,:,l] = convASCIIraster2array(ibound_path, ibound_array[:,:,l], cellsizeMF = delr[0], nrow = nrow, ncol=ncol)
    ibound_array = list(ibound_array)

    strt_array = np.zeros((nrow,ncol, nlay))
    for l in range(nlay):
        strt_path = os.path.join(model_ws, strt_fn[l])
        strt_array[:,:,l] = convASCIIraster2array(strt_path, strt_array[:,:,l], cellsizeMF = delr[0], nrow = nrow, ncol=ncol)
    strt_array = list(strt_array)

    hk_array = np.zeros((nrow,ncol, nlay))
    for l in range(nlay):
        hk_path = os.path.join(model_ws, hk_fn[l])
        hk_array[:,:,l] = convASCIIraster2array(hk_path, hk_array[:,:,l], cellsizeMF = delr[0], nrow = nrow, ncol=ncol)
    hk_array = list(hk_array)

    vka_array = np.zeros((nrow,ncol, nlay))
    for l in range(nlay):
        vka_path = os.path.join(model_ws, vka_fn[l])
        vka_array[:,:,l] = convASCIIraster2array(vka_path, vka_array[:,:,l], cellsizeMF = delr[0], nrow = nrow, ncol=ncol)
    vka_array = list(vka_array)

    ss_array = np.zeros((nrow,ncol, nlay))
    for l in range(nlay):
        ss_path = os.path.join(model_ws, ss_fn[l])
        ss_array[:,:,l] = convASCIIraster2array(ss_path, ss_array[:,:,l], cellsizeMF = delr[0], nrow = nrow, ncol=ncol)
    ss_array = list(ss_array)

    sy_array = np.zeros((nrow,ncol, nlay))
    for l in range(nlay):
        sy_path = os.path.join(model_ws, sy_fn[l])
        sy_array[:,:,l] = convASCIIraster2array(sy_path, sy_array[:,:,l], cellsizeMF = delr[0], nrow = nrow, ncol=ncol)
    sy_array = list(sy_array)

# RECHARGE
    if isinstance(rch_fn,float):
        rch_array = rch_fn
    else:
        rch_array = []
        for n in range(nper):
            rch_path = os.path.join(model_ws, rch_fn + str(n+1) + '.asc')
            rch_array.append(convASCIIraster2array(rch_path, arrayOUT = np.zeros((nrow,ncol)), cellsizeMF = delr[0], nrow = nrow, ncol=ncol))

# DRAIN
    layer_row_column_elevation_cond = [[]]
    ri=0
    for r in top_array:
        ri = ri+1
        ci=0
        for v in r:
            ci=ci+1
            layer_row_column_elevation_cond[0].append([1,ri,ci,v,1e6])
    # in layer 2, cell outlet, elevation is bottom of layer 2
    layer_row_column_elevation_cond[0].append([2,13,3,22,1])

# average for 1st SS stress period
    if dum_sssp1 == 1:
        if isinstance(rch_fn,str):
            rch_array = np.asarray(rch_array)
            rch_SS = np.zeros((nrow,ncol))
            for n in range(nper):
                rch_SS = rch_SS + rch_array[n,:,:]*perlen[n]/sum(perlen)
            rch_array = list(rch_array)
            rch_array.insert(0, rch_SS)
        nper = nper + 1
        perlen.insert(0,1)
        nstp.insert(0,1)
        tsmult.insert(0,1)
        Ss_tr.insert(0, 'SS')

    # 2 - create the modflow packages files

    # MFfile initialization
    mf = modflow(modelname = modelname, exe_name = exe_name, namefile_ext = namefile_ext, version = 'mf2000', model_ws= model_ws)
    # dis initialization
    for i in range(nper):
        if Ss_tr[i] == 'TR':
            Ss_tr[i] = False
        elif Ss_tr[i] == 'SS':
            Ss_tr[i] = True
        else:
            print '\nVariable Ss_tr from the DIS package is not correct, check the MODFLOW manual'
            sys.exit()
    # dis package
    dis = mfdis(model = mf, nrow = nrow, ncol = ncol, nlay = nlay, nper = nper, delr = delr, delc = delc, laycbd = laycbd, top = top_array, botm = botm_array, perlen = perlen, nstp = nstp, tsmult = tsmult, itmuni = itmuni, lenuni = lenuni, steady = Ss_tr, extension = ext_dis)
    # bas package
    bas = mfbas(model = mf, ibound = ibound_array, strt = strt_array, hnoflo = hnoflo, extension = ext_bas)
    # lpf initialization
    lpf = mflpf(model = mf, hdry = hdry, laytyp = laytyp, layavg = layavg, chani = chani, layvka = layvka, laywet = laywet, hk = hk_array, vka = vka_array, ss = ss_array, sy = sy_array, extension=ext_lpf)
    # rch initialization
    rch = mfrch(mf, irchcb=lpf.ilpfcb, nrchop=nrchop, rech=rch_array, extension = ext_rch)
    # drn package initialization
    drn = mfdrn(model = mf, idrncb=lpf.ilpfcb, layer_row_column_elevation_cond = layer_row_column_elevation_cond)
    # ghb package initialization
#    ghb = mfghb(model=mf, igbhcb = lpf.ilpfcb, layer_row_column_head_cond = [[2,13,3,90,50]])
    # output control initialization
    oc = mfoc(mf, ihedfm=ihedfm, iddnfm=iddnfm, item2=[[0,1,1,1]], item3=[[0,0,1,0]], extension=[ext_oc,ext_cbc,ext_heads,ext_ddn])
    # preconditionned conjugate-gradient initialization
    pcg = mfpcg(mf)
    # write packages files
    dis.write_file()
    bas.write_file()
    lpf.write_file()
    rch.write_file()
    drn.write_file()
#    ghb.write_file()
    oc.write_file()
    pcg.write_file()
    h_fn = os.path.join(model_ws, modelname + "." + ext_heads)
    if os.path.exists(h_fn):
        os.remove(h_fn)
    cbc_fn = os.path.join(model_ws, modelname + "." + ext_cbc)
    if os.path.exists(cbc_fn):
        os.remove(cbc_fn)
    # run MODFLOW and read the heads back into Python
    print "\nMODFLOW run...\n"
    mf.write_name_file()
    mf.run_model(pause = False)
    # extract heads
    h = mfrdbin.mfhdsread(mf, 'LF95').read_all(h_fn)
    h_1 = h[1]
    if len(h_1)<sum(perlen):
        raise ValueError, '\nMODFLOW error!\nCheck the MODFLOW list file in folder:\n%s!' % model_ws
    print ''
    # extract cell-by-cell budget
    cbc = mfrdbin.mfcbcread(mf, 'LF95').read_all(cbc_fn)
    cbc_1 = cbc[1]
    cbc_nam = cbc[2]
    if dum_sssp1 == 1:
        h_1 = h_1[1:]
        cbc_1 = cbc_1[1:]
        nper = nper - 1
        perlen = perlen[1:]
    h1 = np.asarray(h_1)
    cbc1 = np.asarray(cbc_1)
    top_array = np.asarray(top_array)

    return SP_d, nrow, ncol, delr, delc, nlay, perlen, nper, np.asarray(top_array), hnoflo, hdry, np.asarray(ibound_array), laytyp[0], h1, cbc1, cbc_nam, top_array, inputFileMF_fn