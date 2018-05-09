    # -*- coding: utf-8 -*-

__author__ = "Alain Francés <frances08512@itc.nl>"
__version__ = "0.3"
__date__ = "2012"

import os
import numpy as np
import matplotlib as mpl

class clsPROCESS:
    def __init__(self, cUTIL, MM_ws, MM_ws_out, MF_ws, nrow, ncol, nlay, xllcorner, yllcorner, cellsizeMF, hnoflo):
        self.cUTIL = cUTIL
        self.MM_ws = MM_ws
        self.MM_ws_out = MM_ws_out
        self.MF_ws = MF_ws
        self.nrow= nrow
        self.ncol= ncol
        self.nlay= nlay
        self.xllcorner= xllcorner
        self.yllcorner=yllcorner
        self.cellsizeMF=cellsizeMF
        self.hnoflo=hnoflo
        self.smMM = []
        self.smMMname = []

    ######################

    def inputEsriAscii(self, grid_fn, datatype, stdout = None, report = None):
        try:
            if datatype == np.int or datatype == int:
                grid_fn = int(grid_fn)
            else:
                grid_fn = float(grid_fn)
            grid_out = np.ones((self.nrow,self.ncol), dtype = datatype)*grid_fn
        except:
            grid_fn=os.path.join(self.MM_ws,grid_fn)
            grid_out=np.zeros([self.nrow,self.ncol], dtype = datatype)
            grid_out = self.convASCIIraster2array(grid_fn,grid_out, stdout = stdout, report = report)
        return grid_out
        del grid_out

    ######################

    def checkarray(self, var, dtype = np.float, stdout = None, report = None):
        try:
            if len(var)>1:
                lst_out = []
                for v in var:
                    if dtype == np.int or dtype == int:
                        lst_out.append(int(v))
                    else:
                        lst_out.append(float(v))
            else:
                if dtype == np.int or dtype == int:
                    lst_out = int(var[0])
                else:
                    lst_out = float(var[0])
        except:
            array = np.zeros((len(var),self.nrow,self.ncol), dtype = dtype)
            l = 0
            for v in var:
                if isinstance(v, str):
                    array_path = os.path.join(self.MF_ws, v)
                    array[l,:,:] = self.convASCIIraster2array(array_path, array[l,:,:], stdout = stdout, report = report)
                else:
                    self.cUTIL.ErrorExit('\nFATAL ERROR!\nMODFLOW ini file incorrect, check files or values %s' % var, stdout = stdout, report = report)
                l += 1
            if len(var)>1:
                lst_out = list(array)
            elif l>1:
                lst_out = list(array[0,:,:])
            else:
                lst_out = array[0,:,:]

        return lst_out

    ######################

    def float2array(self, array):
        if self.nlay < 2:
            if isinstance(array, list):
                array = (np.asarray(array)).reshape((1, self.nrow, self.ncol))
            else:
                array = np.asarray([array])
        else:
            array = np.asarray(array)
        if np.asarray(array).shape[0] == self.nlay and len(array.shape) == 1:
            array_tmp = np.ones([self.nlay, self.nrow, self.ncol], dtype = np.float)
            for l, e in enumerate (array):
                array_tmp[l,:,:] *= e  #*ibound[l,:,:]
            array = array_tmp
        return array
        
    ######################

    def convASCIIraster2array(self, filenameIN, arrayOUT, stdout = None, report = None):
        '''
        Read ESRI/MODFLOW ASCII file and convert to numpy array
        '''

        # Load the grid files
        global NODATA_value, ncol_tmp, nrow_tmp, cellsizeEsriAscii, fin
        if os.path.exists(filenameIN):
            fin = open(filenameIN, 'r')
        else:
            self.cUTIL.ErrorExit("\nFATAL ERROR!\nThe file %s doesn't exist!!!" % filenameIN, stdout = stdout, report = report)

        # test if it is ESRI ASCII file or PEST grid
        line = fin.readline().split()
        testfile = line[0]
        if isinstance(testfile, str):
            # Read the header
            ncol_tmp = int(line[1])
            line = fin.readline().split()
            nrow_tmp = int(line[1])
            line = fin.readline().split()
            line = fin.readline().split()
            line = fin.readline().split()
            cellsizeEsriAscii = float(line[1])
            line = fin.readline().split()
            NODATA_value = line[1]
        elif isinstance(testfile, float):
            ncol_tmp = int(line[0])
            nrow_tmp = int(line[1])

        # verify grid consistency between MODFLOW and ESRI ASCII
        if arrayOUT.shape[0] != self.nrow or arrayOUT.shape[1] != self.ncol or self.cellsizeMF != cellsizeEsriAscii:
            self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nMODFLOW grid and ESRI ASCII grid from file %s don't correspond!\nCheck the cell size and the number of rows, columns and cellsize." % filenameIN, stdout = stdout, report = report)

        # Process the file
#        print "\nConverting %s to np.array" % (filenameIN)
        linecol = 0
        for row in range(nrow_tmp):
        #   if (row % 100) == 0: print ".",
            for col in range(ncol_tmp):
                # Get new data if necessary
                if col == 0 or linecol == len(line):
                    line = fin.readline().split()
                    linecol = 0
                if line[linecol] == NODATA_value:
                    arrayOUT[row,col] = self.hnoflo
                else:
                    arrayOUT[row,col] = line[linecol]
                linecol += 1

        fin.close()
        del line, fin

        return arrayOUT
        del arrayOUT

    ######################

    def procMF(self, cMF, h5_MF, ds_name, ds_name_new, conv_fact = 1.0, index = 0):

        # cbc format is : (kstp), kper, textprocess, nrow, ncol, nlay
        t = 0
        h5_MF.create_dataset(name = ds_name_new, data = np.zeros((sum(cMF.perlen), cMF.nlay, cMF.nrow, cMF.ncol)))
        for n in range(cMF.nper):
            if cMF.perlen[n] != 1:
                for x in range(cMF.perlen[n]):
                    if ds_name == 'heads':
                        h5_MF[ds_name_new][t,:,:,:] = h5_MF['heads'][n,:,:,:]
                    else:
                        array_tmp = h5_MF[ds_name][n,:,:,:,index]
                        if cMF.reggrid == 1 and ds_name != 'heads':
                            h5_MF[ds_name_new][t,:,:,:] = conv_fact*array_tmp/(cMF.delr[0]*cMF.delc[0])
                        else:
                            for i in range(cMF.nrow):
                                for j in range(cMF.ncol):
                                    h5_MF[ds_name_new][t,i,j,:] = conv_fact*array_tmp[i,j,:]/(cMF.delr[j]*cMF.delc[i])
                        del array_tmp
                    t += 1
            else:
                if ds_name == 'heads':
                    h5_MF[ds_name_new][t,:,:,:] = h5_MF['heads'][n,:,:,:]
                else:
                    array_tmp = h5_MF[ds_name][n,:,:,:,index]
                    if cMF.reggrid == 1:
                        h5_MF[ds_name_new][t,:,:,:] = conv_fact*array_tmp/(cMF.delr[0]*cMF.delc[0])
                    else:
                        for i in range(cMF.nrow):
                            for j in range(cMF.ncol):
                                h5_MF[ds_name_new][t,:,i,j] = conv_fact*array_tmp[:,i,j]/(cMF.delr[j]*cMF.delc[i])
                    del array_tmp
                t += 1

    ######################

    def procMM(self, cMF, h5_MM, ds_name, ds_name_new, conv_fact = 1.0):

        t = 0
        h5_MM.create_dataset(name = ds_name_new, data = np.zeros((sum(cMF.perlen), cMF.nrow, cMF.ncol)))
        for n in range(cMF.nper):
            array_tmp = h5_MM[ds_name][n,:,:]
            if cMF.perlen[n] != 1:
                for x in range(cMF.perlen[n]):
                    h5_MM[ds_name_new][t,:,:] = conv_fact*array_tmp
                    t += 1
            else:
                h5_MM[ds_name_new][t,:,:] = conv_fact*array_tmp
                t += 1

    ######################

    def inputSP(self, NMETEO, NVEG, NSOIL, nper,
                inputZON_SP_RF_veg_fn, inputZON_SP_RFe_veg_fn, inputZON_SP_LAI_veg_fn,
                inputZON_SP_PT_fn, inputZON_SP_PE_fn,
                inputZON_SP_E0_fn,
                NFIELD = None,
                inputZON_SP_RF_irr_fn = None, inputZON_SP_RFe_irr_fn = None,
                inputZON_SP_PT_irr_fn = None, input_SP_crop_irr_fn = None,
                stdout = None, report = None):

        # READ input ESRI ASCII rasters vegetation
        global crop_irr_SP, PT_irr_zonesSP, RFe_irr_zonesSP, RF_irr_zonesSP, crop_irr_tmp, PT_irr_tmp, RFe_irr_tmp, RF_irr_tmp, PE_tmp, LAI_veg_tmp, PT_veg_tmp, RFe_veg_tmp, E0, RF_veg
        gridVEGarea_fn=[]
        for v in range(NVEG):
            gridVEGarea_fn.append(os.path.join(self.MM_ws,'inputVEG' + str(v+1)+'area.asc'))
        gridVEGarea=np.zeros([NVEG,self.nrow,self.ncol], dtype=float)
        for v in range(NVEG):
            grid_tmp=np.zeros([self.nrow,self.ncol], dtype=float)
            gridVEGarea[v,:,:]=self.convASCIIraster2array(gridVEGarea_fn[v],grid_tmp, stdout = stdout, report = report)
        gridVEGareatot = np.add.accumulate(gridVEGarea, axis = 0)
        area100_test = gridVEGareatot > 100.0
        if area100_test.sum() > 0:
            self.cUTIL.ErrorExit(msg = '\nFATAL ERROR!\nThe total area of the vegetation in one cell cannot exceed 100.0%!', stdout = stdout, report = report)

        # READ RF and E0 for each zone
        # RF
        RF_veg_fn=os.path.join(self.MM_ws, inputZON_SP_RF_veg_fn)
        if os.path.exists(RF_veg_fn):
            RF_veg = np.loadtxt(RF_veg_fn)
        else:
            self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe file %s doesn't exist!!!" % RF_veg_fn, stdout = stdout, report = report)
        RF_veg_zonesSP = np.zeros([NMETEO,nper], dtype=float)
        # E0
        E0_fn=os.path.join(self.MM_ws, inputZON_SP_E0_fn)
        if os.path.exists(E0_fn):
            E0 = np.loadtxt(E0_fn)
        else:
            self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe file %s doesn't exist!!!" % E0_fn, stdout = stdout, report = report)
        E0_zonesSP = np.zeros([NMETEO,nper], dtype=float)
        for n in range(NMETEO):
            for t in range(nper):
                RF_veg_zonesSP[n,t]=RF_veg[n*nper+t]
                E0_zonesSP[n,t]=E0[n*nper+t]

        # READ PT/RFe for each zone and each vegetation
        # RFe
        RFe_veg_fn = os.path.join(self.MM_ws, inputZON_SP_RFe_veg_fn)
        if os.path.exists(RFe_veg_fn):
            RFe_veg_tmp = np.loadtxt(RFe_veg_fn)
        else:
            self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe file %s doesn't exist!!!" % RFe_veg_fn, stdout = stdout, report = report)
        RFe_veg_zonesSP = np.zeros([NMETEO,NVEG,nper], dtype=float)
        # LAI
        LAI_veg_fn = os.path.join(self.MM_ws, inputZON_SP_LAI_veg_fn)
        if os.path.exists(LAI_veg_fn):
            LAI_veg_tmp = np.loadtxt(LAI_veg_fn)
        else:
            self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe file %s doesn't exist!!!" % LAI_veg_fn, stdout = stdout, report = report)
        LAI_veg_zonesSP = np.zeros([NVEG,nper], dtype=float)
        # PT
        PT_veg_fn = os.path.join(self.MM_ws, inputZON_SP_PT_fn)
        if os.path.exists(PT_veg_fn):
            PT_veg_tmp = np.loadtxt(PT_veg_fn)
        else:
            self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe file %s doesn't exist!!!" % PT_veg_fn, stdout = stdout, report = report)
        PT_veg_zonesSP = np.zeros([NMETEO,NVEG,nper], dtype=float)
        for n in range(NMETEO):
            for v in range(NVEG):
                for t in range(nper):
                    #structure is [number of zones, number of vegetation type, time]
                    RFe_veg_zonesSP[n,v,t] = RFe_veg_tmp[t+(n*NVEG+v)*nper]
                    PT_veg_zonesSP[n,v,t]  = PT_veg_tmp[t+(n*NVEG+v)*nper]
        for v in range(NVEG):
            for t in range(nper):
            #structure is [number of zones, number of vegetation type, time]
                LAI_veg_zonesSP[v,t] = LAI_veg_tmp[t+v*nper]

        # READ PE for each zone and each soil
        PE_fn = os.path.join(self.MM_ws, inputZON_SP_PE_fn)
        if os.path.exists(PE_fn):
            PE_tmp = np.loadtxt(PE_fn)
        else:
            self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe file %s doesn't exist!!!" % PE_fn, stdout = stdout, report = report)
        PE_zonesSP=np.zeros([NMETEO,NSOIL,nper], dtype=float)
        for n in range(NMETEO):
            for s in range(NSOIL):
                for t in range(nper):
                    PE_zonesSP[n,s,t] = PE_tmp[t+(n*NSOIL+s)*nper]
                    #structure is [number of zones, number of soil type, time]

        # read IRRIGATION RF and RFe
        if NFIELD != None:
            # RF_irr
            RF_irr_fn = os.path.join(self.MM_ws, inputZON_SP_RF_irr_fn)
            if os.path.exists(RF_irr_fn):
                RF_irr_tmp = np.loadtxt(RF_irr_fn)
            else:
                self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe file %s doesn't exist!!!" % RF_irr_fn, stdout = stdout, report = report)
            RF_irr_zonesSP = np.zeros([NMETEO,NFIELD,nper], dtype=float)
            # RFe irr
            RFe_irr_fn = os.path.join(self.MM_ws, inputZON_SP_RFe_irr_fn)
            if os.path.exists(RFe_irr_fn):
                RFe_irr_tmp = np.loadtxt(RFe_irr_fn)
            else:
                self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe file %s doesn't exist!!!" % RFe_irr_fn, stdout = stdout, report = report)
            RFe_irr_zonesSP = np.zeros([NMETEO,NFIELD,nper], dtype=float)
            # PT irr
            PT_irr_fn = os.path.join(self.MM_ws, inputZON_SP_PT_irr_fn)
            if os.path.exists(PT_irr_fn):
                PT_irr_tmp = np.loadtxt(PT_irr_fn)
            else:
                self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe file %s doesn't exist!!!" % PT_irr_fn, stdout = stdout, report = report)
            PT_irr_zonesSP = np.zeros([NMETEO,NFIELD,nper], dtype=float)
            # crop irr
            crop_irr_fn = os.path.join(self.MM_ws, input_SP_crop_irr_fn)
            if os.path.exists(crop_irr_fn):
                crop_irr_tmp = np.loadtxt(crop_irr_fn)
            else:
                self.cUTIL.ErrorExit(msg = "\nFATAL ERROR!\nThe file %s doesn't exist!!!" % crop_irr_fn, stdout = stdout, report = report)
            crop_irr_SP = np.zeros([NFIELD,nper], dtype=int)
            for n in range(NMETEO):
                for f in range(NFIELD):
                    for t in range(nper):
                        #structure is [number of zones, number of field, time]
                        RF_irr_zonesSP[n,f,t]  = RF_irr_tmp[t+(n*NFIELD+f)*nper]
                        RFe_irr_zonesSP[n,f,t] = RFe_irr_tmp[t+(n*NFIELD+f)*nper]
                        PT_irr_zonesSP[n,f,t]  = PT_irr_tmp[t+(n*NFIELD+f)*nper]
            for f in range(NFIELD):
                for t in range(nper):
                    #structure is [number of field, time]
                    crop_irr_SP[f,t] = crop_irr_tmp[t+f*nper]

        if NFIELD == None:
            return gridVEGarea, RF_veg_zonesSP, E0_zonesSP, PT_veg_zonesSP, RFe_veg_zonesSP, LAI_veg_zonesSP, PE_zonesSP
        else:
            return gridVEGarea, RF_veg_zonesSP, E0_zonesSP, PT_veg_zonesSP, RFe_veg_zonesSP, LAI_veg_zonesSP, PE_zonesSP, RF_irr_zonesSP, RFe_irr_zonesSP, PT_irr_zonesSP, crop_irr_SP

    ######################

    def inputSoilParam(self, SOILparam_fn, NSOIL, stdout = None, report = None):

        # Soils parameter initialisation
        nam_soil=[]
        nsl=[]
        st =[]
        slprop=[]
        Sm =[]
        Sfc=[]
        Sr =[]
        Si =[]
        Ks =[]

        # soil parameters file
        inputFile = self.cUTIL.readFile(self.MM_ws,SOILparam_fn)
        SOILzones=int(int(inputFile[0]))
        if SOILzones>NSOIL:
            print '\nWARNING!\n' + str(SOILzones) + ' soil parameters groups in file [' + SOILparam_fn + ']\n Only ' + str(NSOIL) + ' PE time series found.'
        # trick to initialise the reading position in the next loop
        nsl.append(0)
        for i in range(SOILzones):
            nsl.append(int(inputFile[i+1]))
            if nsl[i+1]<1:
                self.cUTIL.ErrorExit(msg = '\nFATAL ERROR!\nThe model requires at least 1 soil layer!', stdout = stdout, report = report)

        # soil parameter definition for each soil type
        nslst = SOILzones+1
        nam_soil = []
        st = []
        for z in range(SOILzones):
            nam_soil.append(inputFile[nslst].strip())
            nslst += 1
            st.append(inputFile[nslst].strip())
            nslst += 1
            slprop.append([])
            Sm.append([])
            Sfc.append([])
            Sr.append([])
            Si.append([])
            Ks.append([])
            for ns in range(nsl[z+1]):
                slprop[z].append(float(inputFile[nslst]))
                nslst += 1
                Sm[z].append(float(inputFile[nslst]))
                nslst += 1
                Sfc[z].append(float(inputFile[nslst]))
                nslst += 1
                Sr[z].append(float(inputFile[nslst]))
                nslst += 1
                Si[z].append(float(inputFile[nslst]))
                nslst += 1
                Ks[z].append(float(inputFile[nslst]))
                nslst += 1
                if not(Sm[z][ns]>Sfc[z][ns]>Sr[z][ns]) or not(Sm[z][ns]>=Si[z][ns]>=Sr[z][ns]):
                    self.cUTIL.ErrorExit('\nFATAL ERROR!\nSoils parameters are not valid!\nThe conditions are Sm>Sfc>Sr and Sm>Si>Sr!', stdout = stdout, report = report)
            if sum(slprop[z][0:nsl[z+1]])>1:
                self.cUTIL.ErrorExit('\nFATAL ERROR!\nThe sum of the soil layers proportion of %s is >1!\nCorrect your soil data input!\n' % nam_soil[z], stdout = stdout, report = report)
            if sum(slprop[z][0:nsl[z+1]])<1:
                self.cUTIL.ErrorExit('\nFATAL ERROR!\nThe sum of the soil layers proportion of %s is <1!\nCorrect your soil data input!\n' % nam_soil[z], stdout = stdout, report = report)

        return nsl[1:len(nsl)], nam_soil, st, slprop, Sm, Sfc, Sr, Si, Ks
        del nsl, nam_soil, st, slprop, Sm, Sfc, Sr, Si, Ks

    ######################

    def inputObs(self, inputObs_fn, inputObsHEADS_fn, inputObsSM_fn, inputObsRo_fn, inputDate, _nslmax, nlay, stdout = None, report = None):
        '''
        observations cells for soil moisture and heads (will also compute SATFLOW)
        '''

        # read coordinates and SATFLOW parameters
        global name
        inputFile = self.cUTIL.readFile(self.MM_ws,inputObs_fn)

        # define a dictionnary of observations,  format is: Name (key) x y i j hi h0 RC STO
        obs = {}
        obs_list = []
        for o in range(len(inputFile)):
            line = inputFile[o].split()
            name = line[0]
            obs_list.append(name)
            x = float(line[1])
            y = float(line[2])
            lay = int(line[3])
            try:
                hi  = float(line[4])
                h0  = float(line[5])
                RC  = float(line[6])
                STO = float(line[7])
            except:
                hi = h0 = RC = STO =  self.hnoflo
            try:
                lbl  = str(line[0])
            except:
                lbl = ''
            # verify if coordinates are inside MODFLOW grid
            if (x <= self.xllcorner or
               x > (self.xllcorner+self.ncol*self.cellsizeMF) or
               y <= self.yllcorner or
               y > (self.yllcorner+self.nrow*self.cellsizeMF)):
                   print 'WARNING!\nObservation point %s has coordinates outside the MODFLOW grid and will not be considered.' % name
            else:
                if lay > nlay or lay < 1:
                    print 'WARNING!\nLayer %s of observation point %s is not valid (corrected to layer 1)!\nCheck your file %s (layer number should be between 1 and the number of layer of the MODFLOW model, in this case %s).' % (lay, name, inputObs_fn, nlay)
                    lay = 0
                else:
                    lay = lay - 1
                # compute the coordinates in the MODFLOW grid
                #TODO use the PEST utilities for space extrapolation
                i = int(self.nrow - np.ceil((y-self.yllcorner)/self.cellsizeMF))
                j = int(np.ceil((x-self.xllcorner)/self.cellsizeMF)-1)
                #  read obs time series
                obsh_fn = os.path.join(self.MM_ws, '%s_%s.txt' % (inputObsHEADS_fn, name))
                if os.path.exists(obsh_fn):
                    obs_h, obs_h_yn = self.verifObs(inputDate, obsh_fn, obsnam = name, stdout = stdout, report = report)
                else:
                    obs_h = []
                    obs_h_yn = 0
                obssm_fn=os.path.join(self.MM_ws, '%s_%s.txt' % (inputObsSM_fn, name))
                if os.path.exists(obssm_fn):
                    obs_sm, obs_sm_yn = self.verifObs(inputDate, obssm_fn, _nslmax, obsnam = name, stdout = stdout, report = report)
                else:
                    obs_sm = []
                    obs_sm_yn = 0
                obsRo_fn=os.path.join(self.MM_ws, '%s_%s.txt' % (inputObsRo_fn, name))
                if os.path.exists(obsRo_fn):
                    obs_Ro, obs_Ro_yn = self.verifObs(inputDate, obsRo_fn, obsnam = name, stdout = stdout, report = report)
                else:
                    obs_Ro = []
                    obs_Ro_yn = 0
                obs[name] = {'x':x,'y':y,'i': i, 'j': j, 'lay': lay, 'hi':hi, 'h0':h0, 'RC':RC, 'STO':STO, 'lbl':lbl, 'obs_h':obs_h, 'obs_h_yn':obs_h_yn, 'obs_SM':obs_sm, 'obs_sm_yn':obs_sm_yn, 'obs_Ro':obs_Ro, 'obs_Ro_yn':obs_Ro_yn}

        #  read catchment obs time series
        obsh_fn = os.path.join(self.MM_ws, '%s_catchment.txt' % inputObsHEADS_fn)
        if os.path.exists(obsh_fn):
            obs_h, obs_h_yn = self.verifObs(inputDate, obsh_fn, obsnam = name, stdout = stdout, report = report)
        else:
            obs_h = []
            obs_h_yn = 0
        obssm_fn = os.path.join(self.MM_ws, '%s_catchment.txt' % inputObsSM_fn)
        if os.path.exists(obssm_fn):
            obs_sm, obs_sm_yn = self.verifObs(inputDate, obssm_fn, _nslmax, obsnam = name, stdout = stdout, report = report)
        else:
            obs_sm = []
            obs_sm_yn = 0
        obsRo_fn = os.path.join(self.MM_ws, '%s_catchment.txt' % inputObsRo_fn)
        if os.path.exists(obsRo_fn):
            obs_Ro, obs_Ro_yn = self.verifObs(inputDate, obsRo_fn, obsnam = name, stdout = stdout, report = report)
        else:
            obs_Ro = []
            obs_Ro_yn = 0
        obs_catch = {}
        obs_catch['catch'] = {'obs_h':obs_h, 'obs_h_yn':obs_h_yn, 'obs_SM':obs_sm, 'obs_sm_yn':obs_sm_yn, 'obs_Ro':obs_Ro, 'obs_Ro_yn':obs_Ro_yn}
        obs_catch_list = [obs_h_yn, obs_sm_yn, obs_Ro_yn]

        return obs, obs_list, obs_catch, obs_catch_list
        del inputObs_fn, inputObsHEADS_fn, inputObsSM_fn, inputDate, _nslmax
        del obs

    ######################

    def verifObs(self, inputDate, filename, _nslmax = 0, obsnam = 'unknown location', stdout = None, report = None):
        '''
        Import and process data and parameters
        '''

        global obs_yn, obsDate, obsValue
        if os.path.exists(filename):
            try:
                obsData=np.loadtxt(filename, dtype = str)
                obsDate = mpl.dates.datestr2num(obsData[:,0])
                obsValue = []
                for l in range(1,len(obsData[0])):
                    obsValue.append(obsData[:,l].astype(float))
                del obsData
                if len(obsValue)<_nslmax:
                    for l in range(_nslmax-len(obsValue)):
                        obsValue.append(self.hnoflo)
            except:
                self.cUTIL.ErrorExit(msg = '\nFATAL ERROR!\nFormat of observation file uncorrect!\n%s' % filename, stdout = stdout, report = report)
            obsOutput = np.ones([len(obsValue),len(inputDate)], dtype=float)*self.hnoflo
            obs_yn = 0
            if (obsDate[len(obsDate)-1] < inputDate[0]) or (obsDate[0] > inputDate[len(inputDate)-1]):
                obsOutput = []
                print '\nObservations of file %s outside the modeling period!' % filename
            else:
                for l in range(len(obsValue)):
                    if not isinstance(obsValue[l], float):
                        j=0
                        if inputDate[0] >= obsDate[0]:
                           while obsDate[j]<inputDate[0]:
                                j += 1
                        for i in range(len(inputDate)):
                            if j<len(obsDate):
                                if inputDate[i]==obsDate[j]:
                                    if obsValue[l][j] > 0.0:
                                        obsOutput[l,i]=obsValue[l][j]
                                    else:
                                        obsOutput[l,i]=self.hnoflo
                                    j += 1
                                    obs_yn = 1
                                else:
                                    obsOutput[l,i]=self.hnoflo
                            else:
                                obsOutput[l,i]=self.hnoflo
                    else:
                        obsOutput[l,:] = np.ones([len(inputDate)])*self.hnoflo
        else:
            obsOutput = []
        return obsOutput, obs_yn
        del inputDate, obsOutput

    ######################

    def outputEAgrd(self, outFile_fn, outFolder = []):

        if outFolder == []:
            outFolder = self.MM_ws

        outFile=open(os.path.join(outFolder, outFile_fn), 'w')

        outFile=self.writeHeaderESRIraster(outFile)

        return outFile
        del outFile

    ######################

    def writeHeaderESRIraster(self, file_asc):
        '''
        Write ESRI ASCII raster header
        '''
        file_asc.write('ncols  ' + str(self.ncol)+'\n' +
                       'nrows  '  + str(self.nrow)+'\n' +
                       'xllcorner  '+ str(self.xllcorner)+'\n' +
                        'yllcorner  '+ str(self.yllcorner)+'\n' +
                        'cellsize   '+ str(self.cellsizeMF)+'\n' +
                        'NODATA_value  '+ str(self.hnoflo)+'\n')
        return file_asc
        del file_asc

    #####################################

    def ExportResultsMM4PEST(self, i, j, inputDate, _nslmax, results_S, index_S, obs_S, obsname):
        """
        Export the processed data in a txt file
        INPUTS:      output fluxes time series and date
        OUTPUT:      ObsName.txt
        """
        self.smMM.append([])
        len_smMM = len(self.smMM)
        for t in range(len(inputDate)):
            out_date = mpl.dates.num2date(inputDate[t]).strftime("%d/%m/%Y")
            for l in range(_nslmax):
                try:
                    obs_S_tmp = obs_S[l,t]
                except:
                    obs_S_tmp = -1.0
                if results_S[t,l,index_S.get('iSsoil_pc_s')] > 0.0 and obs_S_tmp > 0.0:
                    self.smMM[len_smMM-1].append((obsname+'SM_l'+str(l+1)).ljust(12,' ')+ out_date.ljust(12,' ')+ ' 00:00:00 ' + str(results_S[t,l,index_S.get('iSsoil_pc_s')]).ljust(10,' ') + '\n')
        del i, j, _nslmax, results_S, index_S, obs_S, obsname

    #####################################

    def ExportResultsPEST(self, i, j, inputDate, _nslmax, obs_h, obs_S, outPESTheads, outPESTsm, obsname):
        """
        Export the obs data in a txt file in PEST format
        INPUTS:      output fluxes time series and date
        OUTPUT:      PESTxxx.smp
        """
        self.smMM.append([])
        len_smMM = len(self.smMM)
        for t in range(len(inputDate)):
            date = mpl.dates.num2date(inputDate[t]).strftime("%d/%m/%Y")
            try:
                if obs_h[t] != self.hnoflo:
                    outPESTheads.write(obsname.ljust(10,' ')+ date.ljust(10,' ')+ ' 00:00:00 ' + str(obs_h[t]).ljust(10,' ') + '\n')
            except:
                pass
            try:
                for l in range (_nslmax):
                    if obs_S[l,t] != self.hnoflo:
                        self.smMM[len_smMM-1].append((obsname+'SM_l'+str(l+1)).ljust(10,' ')+ date.ljust(10,' ')+ ' 00:00:00 ' + str(obs_S[l,t]).ljust(10,' ') + '\n')
            except:
                pass
    #####################################

    def compRMSE (self, sim, obs) :
        def sqre_diff (v, w) :
            return (v - w) ** 2
        s = len(sim)
        v = sum(map(sqre_diff, sim, obs))
        return np.sqrt(v / s)

    #####################################

    def compE (self, sim, obs, hnoflo) :
        def sqre_diff (v, w) :
            return (v - w) ** 2
        #s = len(sim)
        if sum(sqre_diff(obs, self.compAVGE(obs)))>0.0:
            v = 1 - sum(map(sqre_diff, sim, obs))/sum(sqre_diff(obs, self.compAVGE(obs)))
        else:
            v = hnoflo
        return v

    #####################################

    def compAVGE(self, x):
        assert len(x) > 0
        return float(sum(x)) / len(x)

    #####################################

    def compR(self, x, y, hnoflo):
        # pearson correlation coefficient r
        assert len(x) == len(y)
        n = len(x)
        assert n > 0
        avg_x = self.compAVGE(x)
        avg_y = self.compAVGE(y)
        diffprod = 0
        xdiff2 = 0
        ydiff2 = 0
        for idx in range(n):
            xdiff = x[idx] - avg_x
            ydiff = y[idx] - avg_y
            diffprod += xdiff * ydiff
            xdiff2 += xdiff * xdiff
            ydiff2 += ydiff * ydiff
        if np.sqrt(xdiff2 * ydiff2)> 0.0 :
            return diffprod / np.sqrt(xdiff2 * ydiff2)
        else:
            return hnoflo

    #####################################

    def compCalibCrit(self, sim, obs, hnoflo): 
        rmse = None
        rsr  = None
        nse  = None
        r    = None
        a = np.array([sim,obs])
        a = np.transpose(a)
        if hnoflo > 0:
            b = a[~(a > hnoflo - 1000.0).any(1)]
        else:
            b = a[~(a < hnoflo + 1000.0).any(1)]
        if b[:,0] <> []:
            rmse = (self.compRMSE(b[:,0], b[:,1]))
            if np.std(b[:,1]) > 0:
                rsr = (rmse/(np.std(b[:,1])))
            nse = (self.compE(b[:,0], b[:,1], hnoflo))
            r = (self.compR(b[:,0], b[:,1], hnoflo))
        return rmse, rsr, nse, r

#####################################

    def compCalibCritObs(self, Spc, h_MF, Sobs, hobs, hnoflo, obs_name, nsl, h_MM = None):
        
        global l
        rmseSM = None
        rmseHEADS = None
        rmseHEADSc = None
        rsrSM = None
        rsrHEADS = None
        rsrHEADSc = None
        nseSM = None
        nseHEADS = None
        nseHEADSc = None
        rSM = None
        rHEADS = None
        rHEADSc = None

        testSM = 0
        Sobs_m = []
        Spc1full = []
        for l in range(nsl):
            Spc1full.append(Spc[:,l])
            try:
                Sobs_m.append(np.ma.masked_values(Sobs[l,:], hnoflo, atol = 0.09))
                testSM += 1
            except:
                Sobs_m.append([])
        if testSM > 0 or hobs <> []:
            print 'RMSE/RSR/NSE/r'
            if testSM > 0:
                rmseSM = []
                rsrSM = []
                nseSM = []
                rSM = []
                try:
                    for l, (y, y_obs) in enumerate(zip(Spc1full, Sobs_m)):
                        if y_obs <> []:
                            rmse, rsr, nse, r = self.compCalibCrit(y,y_obs, hnoflo)
                            rmseSM.append(100.0*rmse)
                            rsrSM.append(rsr)
                            nseSM.append(nse)
                            rSM.append(r)
                            del rmse, rsr, nse, r
                            if rmseSM[l] != None:
                                print 'SM layer %d: %.1f %% / %.2f / %.2f / %.2f' % (l+1, rmseSM[l], rsrSM[l], nseSM[l], rSM[l])
                except:
                    print 'SM layer %d: error' % (l+1)
            if hobs <> []:
                try:
                    rmse, rsr, nse, r = self.compCalibCrit(h_MF,hobs, hnoflo)
                    rmseHEADS = [rmse]
                    rsrHEADS  = [rsr]
                    nseHEADS  = [nse]
                    rHEADS    = [r]
                    del rmse, rsr, nse, r
                    if rmseHEADS[0] != None:
                        print 'h: %.2f m / %.2f / %.2f / %.2f' % (rmseHEADS[0], rsrHEADS[0], nseHEADS[0], rHEADS[0])
                except:
                    print 'h: error'
                if h_MM is not None:
                    try:
                        rmse, rsr, nse, r = self.compCalibCrit(h_MM,hobs,hnoflo)
                        rmseHEADSc = [rmse]
                        rsrHEADSc = [rsr]
                        nseHEADSc = [nse]
                        rHEADSc = [r]
                        del rmse, rsr, nse, r
                        if rmseHEADSc[0] != None:
                            print 'hcorr: %.2f m / %.2f / %.2f / %.2f' % (rmseHEADSc[0], rsrHEADSc[0], nseHEADSc[0], rHEADSc[0])
                    except:
                        print 'hcorr: error'                    

        return rmseHEADS, rmseHEADSc, rmseSM, rsrHEADS, rsrHEADSc, rsrSM, nseHEADS, nseHEADSc, nseSM, rHEADS, rHEADSc, rSM

##################
if __name__ == "__main__":
    print '\nWARNING!\nStart MARMITES-MODFLOW models using the script startMARMITES_v3.py\n'

# EOF