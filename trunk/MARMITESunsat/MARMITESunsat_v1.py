# -*- coding: utf-8 -*-
"""
DLPWATFLUX stands for Distributed depth-wise Lumped-Parameter model for
spatio-temporal assessment of WATer FLUXes in the unsaturated zone.
The familiar and friendly name of DLPWATFLUX is MARMITES, a french word to
design a big cooking pot used by sorcerers for all kinds of experiments!
The main objective of DPLWATFLUX development is to partition the rainfall
in the several fluxes from the unsaturated and saturated zone (source
reservoir). It applies the concepts enunciated by Lubcynski (2010):                        ## TO BE UPDATED## 20100823
ET = Es + I + ETu + ETg
ETu = Eu + Tu
ETg = Eg + Tg
DLPWATFLUX is a modular algorithm that computes on a daily temporal
basis: interception, surface storage and runoff, evaporation from the
unsaturated and saturated zone (from MODFLOW), soil moisture storage,
aquifer recharge.
Input driving forces are rainfall and potential evapotranspiration daily
time series.
Calibration is made against soil moisture (DPLWATFLUX) and hydraulic
heads (MODFLOW).
It is link to MODFLOW2000 (mf2k) through PEST (Doherty 2010).
DLPWATFLUX imports and exports mf2k packages for spatial and temporal
discretization, intial heads, recharge, etc.
Thanks to Toews 2007 and flopy.                                         ## TO BE UPDATED## 20100823
For more documentation, see Frances et al 2010.                         ## TO BE UPDATED## 20100823
References:
Lub 2010
Frances et al 2010                                                      ## TO BE UPDATED## 20101015
"""

__author__ = "Alain Francés <frances08512@itc.nl>"
__version__ = "1.0"
__date__ = "November 2010"

import os
import numpy as np
import pylab

class UNSAT:

    """
    unsat class: compute the water balance in soil (depth-wise, each reservoir
    correspond to a soil horizon (A, B, C and rock)
    ## TO BE UPDATED""
    _______________________________________________________________________________

    INPUTS
            PARAMETERS
                MAXIL       Interception looses
                Sm          Max. Soil moiture storage
                Sfc         Soil moiture storage at field capacity
                Sr          Residual Soil moiture storage
                Si          Inicial soil moisture
                Ks          Saturated hydraulic condutivity
                SUSTm     Max. surface storage
            STATE VARIABLES
                RF           Daily rainfall
                PET         Daily evapotranspiration
    OUTPUTS
            RFe              Daily Excess rainfall
            ETu             Daily evapotranspiration
            S               Daily soil moisture
            Rp              Daily percolation
            SUST            Daily ponding
            Qs              Daily runoff

    Provide daily values for the LINRES module
    ______________________________________________________________________________

    """

    def __init__(self, hnoflo):
        self.hnoflo = hnoflo


    def _partition(self, SUSTprev,PET,E0,Sprev,
                        D, RFe,Sm,Sfc,Sr,Ks,SUSTm):

        '''
        Ponding and surface runoff function
        '''

        def pond(s_tmp,D,Sm,SUSTm):
            countPOND_tmp = 0
            countRunoff_tmp = 0
            if (s_tmp-Sm)>0.000001:
                sust_tmp =D*(s_tmp-Sm)
                if SUSTm>0.0:
                    countPOND_tmp = 1
                if sust_tmp-SUSTm >= 0.000001:
                    qs_tmp =(sust_tmp-SUSTm)
                    sust_tmp = SUSTm
                    countRunoff_tmp = 1
                else:
                    qs_tmp =(0.0)
            else:
                sust_tmp =(0.0)
                qs_tmp =(0.0)
            return sust_tmp, qs_tmp, countPOND_tmp, countRunoff_tmp
            del sust_tmp, qs_tmp, s_tmp, countPOND_tmp, countRunoff_tmp

        def perc(s_tmp,D,Sm,Sfc,Ks):
            '''
            Percolation function
            #TODO see also SWAT PERCOLATION pag 150 chap 2:3.2
            '''
            # Percent. of gravitational water
            Sg=(s_tmp-Sfc)/(Sm-Sfc)
            if s_tmp-Sfc<=0.000001:
                rp_tmp= (0.0)
            elif (s_tmp-Sfc)>0.000001 and (s_tmp-Sm)<=0.000001:
                if (Ks*Sg-D*(s_tmp-Sfc))> 0.000001:
                    rp_tmp= D*(s_tmp-Sfc)
                else:
                    rp_tmp= Ks*Sg
            elif (s_tmp-Sm)>0.000001:
                if Ks-D*(Sm-Sfc)>0.000001:
                    rp_tmp= D*(Sm-Sfc)
                else:
                    rp_tmp = Ks
            return rp_tmp
            del Sg, s_tmp, rp_tmp

        def evp(s_tmp,pet,D,Sm,Sr):
            '''
            Actual evapotranspiration function
            '''
            # Percent. of soil saturation
            Se=(s_tmp-Sr)/(Sm-Sr)
            if s_tmp-Sr<=0.000001:
                evp_tmp= 0.0
            elif s_tmp-Sr>0.000001 and s_tmp-Sm<=0.000001:
                if (pet*Se-(D*(s_tmp-Sr)))> 0.000001:
                    evp_tmp= D*(s_tmp-Sr)
                else:
                    evp_tmp= pet*Se
            elif (s_tmp-Sm)>0.000001:
                if pet-D*(Sm-Sr)>0.000001:
                    evp_tmp= D*(Sm-Sr)
                else:
                    evp_tmp= pet
            return evp_tmp
            del s_tmp, pet

        # MAIN
        # soil reservoir water content initialisation
        # test if SUST can infiltrate in soil
        # compute Es from SUST and PET
        if E0-SUSTprev>=0.000001:
            S_tmp = Sprev*D + RFe
            ETs_tmp = SUSTprev
        else:
            S_tmp = Sprev*D + RFe + SUSTprev-E0
            ETs_tmp = E0

        # SUST and Qs
        PONDtmp = pond(S_tmp/D,D,Sm,SUSTm)
        SUST_tmp=PONDtmp[0]
        Qs_tmp=PONDtmp[1]
        countPOND_tmp=PONDtmp[2]
        countRunoff_tmp=PONDtmp[3]
        S_tmp=S_tmp-(SUST_tmp+Qs_tmp)

        # Rp
        Rp_tmp=perc(S_tmp/D,D,Sm,Sfc,Ks)
        S_tmp=S_tmp-Rp_tmp

        # ETu
        ETu_tmp=evp(S_tmp/D,PET,D,Sm,Sr)
        S_tmp=S_tmp-ETu_tmp

        return (SUST_tmp, Qs_tmp, countPOND_tmp, countRunoff_tmp, Rp_tmp, ETu_tmp, S_tmp, ETs_tmp)
        del SUST_tmp, Qs_tmp, Rp_temp, ETu_tmp, S_tmp

    def run(self, i, j,
                  Sm, Sfc, Sr, Si, D, Ks, SUSTm,
                  ELEV, HEADS,
                  RF, E0, PETveg, RFeveg, PEsoil, VEGarea,
                  perlen, AqType, hdry):

        # Output initialisation
        Ttotal=len(RF)
        # PET for the vegetation patchwork
        PET_tot=np.zeros([len(PETveg[0])], dtype=float)
        # RFe for the vegetation patchwork
        RFe_tot=np.zeros([len(RFeveg[0])], dtype=float)
        # Surface storage initialisation
        SUST=np.zeros([Ttotal], dtype=np.float)
        # Surface runoff initialisation
        Qs=np.zeros([Ttotal], dtype=np.float)
        #Deep percolation initialisation
        Rp=np.zeros([Ttotal], dtype=np.float)                             #Soil moisture storage initialisation
        S=np.zeros([Ttotal], dtype=np.float)                              #Actual evapotranspiration from unsaturated zone initialisation
        ETu=np.zeros([Ttotal], dtype=np.float)
        #Actual evapotranspiration from surface initialisation
        Es=np.zeros([Ttotal], dtype=np.float)
        #Flood frequency (saturated overland flow)
        countFLOOD=np.zeros(Ttotal, dtype=int)
        #soil partial saturation by GW frequency
        countSATpart=np.zeros(Ttotal, dtype=int)
        #MASS BALANCE
        MB=np.zeros(Ttotal, dtype=float)
        # hortonian overland flow
        countPOND=np.zeros(Ttotal, dtype=int)
        countRunoff=np.zeros(Ttotal, dtype=int)
        INTER=np.zeros([Ttotal], dtype=np.float)

        if D < 0.005: #correct soil thickness less than 5cm
            D=0.005
        Dbot= (ELEV-D)*1000        # compute elevation of the bottom soil layer
        D = D*1000

        # PROCESSING THE WHOLE DATA SET
        for t in range(int(perlen)):    # t: current time step

            # Preprocessing of PET/PE/INTER and RFe for different vegetation
            SOILarea = 100
            for v in range(len(PETveg)):
                if VEGarea[v]<>self.hnoflo:
                    PET_tot[t]=PET_tot[t]+PETveg[v,t]*VEGarea[v]/100
                    SOILarea = SOILarea - VEGarea[v]
            PET_tot[t] = PET_tot[t] + PEsoil[t]*SOILarea/100

            SOILarea = 100
            for v in range(len(RFeveg)):
                if VEGarea[v]<>self.hnoflo:
                    RFe_tot[t]=RFe_tot[t]+RFeveg[v,t]*VEGarea[v]/100
                    SOILarea = SOILarea - VEGarea[v]
            RFe_tot[t] = RFe_tot[t] + RF[t]*SOILarea/100

            INTER[t] = RF[t] - RFe_tot[t]

            # test if first time step, if yes use Si
            if t>0:
                Sprev=S[t-1]
                SUSTprev=SUST[t-1]
            else:
                Sprev=Si
                SUSTprev=0

            # handle drycell
            if HEADS[t]==hdry:
                HEADStmp = Dbot - 1000
            else:
                HEADStmp=HEADS[t]*1000

            if AqType==1:
                # AQUIFER UNCONFINED, water table can rise in the soil and above surface
                # heads below soil bottom
                if HEADStmp-Dbot<=0.000001:
                    if countFLOOD[t]-countFLOOD[t-1]==-1:
                        Sprev=Sm
                    SUSTtmp, Qstmp, countPONDtmp, countRunofftmp, Rptmp, ETutmp, Stmp, ETstmp=self._partition(SUSTprev,PET_tot[t],E0[t],Sprev, D, RFe_tot[t],Sm,Sfc,Sr,Ks,SUSTm)
                else:
                    # heads above soil bottom and below surface
                    if (HEADStmp-Dbot>0.000001 and ELEV-HEADStmp>0.000001):
                        if countFLOOD[t]-countFLOOD[t-1]==-1:
                            Sprev=Sm
                        countSATpart[t] = 1
                        D=ELEV-HEADStmp
                        SUSTtmp, Qstmp, Rptmp, ETutmp, Stmp, ETstmp=self._partition(SUSTprev,PET_tot[t],E0[t],Sprev,D, RFe_tot[t],Sm,Sfc,Sr,Ks,SUSTm)
                    # heads above surface
                    else:
                        countFLOOD[t]=1
                        SUSTtmp=HEADStmp-ELEV+SUSTprev
                        if SUSTtmp > SUSTm:
                            Qstmp =(SUSTtmp-SUSTm)
                            SUSTtmp = SUSTm
                        else:
                            Qstmp =(0.0)
                        ETutmp=0.0
                        Rptmp=0.0
                        Stmp=Sm*D
            else:
            # AQUIFER CONFINED, water table cannot rise in the soil or above topography
                SUSTtmp, Qstmp, countPONDtmp, countRunofftmp, Rptmp, ETutmp, Stmp, ETstmp=self._partition(SUSTprev,PET_tot[t],E0[t],Sprev, D, RFe_tot[t],Sm,Sfc,Sr,Ks,SUSTm)

            # fill the table and compute water balance
            SUST[t]=SUSTtmp
            Qs[t]=Qstmp
            Rp[t]=Rptmp
            ETu[t]=ETutmp
            Es[t]=ETstmp
            countPOND[t]=countPONDtmp
            countRunoff[t]=countRunofftmp
            S[t]=Stmp/D
            if countFLOOD[t]==0:
                MB[t]=RF[t]+SUSTprev-INTER[t]-Rp[t]-Qs[t]-SUST[t]-ETu[t]-Es[t]-((S[t]-Sprev)*D)
            else:
                MB[t]=HEADStmp-ELEV-SUST[t]-Qs[t]+SUSTprev

        return RF, PET_tot, RFe_tot, SUST, Qs, ETu, S, Rp, countFLOOD, Es, countSATpart, MB, INTER, countPOND, countRunoff, E0
        del RF, PET_tot, RFe_tot, SUST, Qs, ETu, S, Rp, countFLOOD, Es, countSATpart, MB, R, countPOND, countRunoff, Spercent, INTER

#######################################################

class LINRES:
    """
    LINRES: LINear REServoirs
    Calculate recharge R in function o percolation (Rp processed in SOMOS)
    _______________________________________________________________________________

    INPUTS
            PARAMETERS
                n           Number of reservoirs
                f           Unsaturated recession constant
            STATE VARIABLES
                Rp          Daily percolation
    OUTPUTS
            R               Daily recharge

    Provide daily reharge for the SATFLOW module
    ______________________________________________________________________________
    ______________________________________________________________________________
    """

#    def __init__(self):
#        self.n = n
#        self.f = f

    def run(self, Rp, n, f):

        Y=np.zeros((len(Rp), n+1), float)

        # Initialization of the first line of the array
        Y[0,0] = (1+f)*Rp[0]/f

        #  loop
        for t in range(1,len(Rp)):
            for i in range(n+1):
                if i==0:
                    Y[t,i] = (1+f)*Rp[t]/f
                else:
                    for j in range(0,i+1):
                        Y[t,i] = Y[t,i]+((1+f)**(-j))*Y[t-1,i-j]
                    Y[t,i] = Y[t,i] * f/(1+f)
        R = Y[:,n]

        return R

        del R, Y

#######################################################

class SATFLOW:
    """
    SATFLOW: SATurated FLOW
    Calculate water level fluctuations in function of recharge (R processed in LINRES)
    _______________________________________________________________________________

    INPUTS
            PARAMETERS
                hi
                h0          base level
                RC          Saturated recession constant
                STO         Aquifer Storage capacity
            STATE VARIABLES
                R           Daily recharge
    OUTPUTS
            h               Daily water level
    ______________________________________________________________________________
    ______________________________________________________________________________
    """

##    def __init__(self):
##        self.hi = hi
##        self.h0 = h0
##        self.RC = RC
##        self.STO = STO

    def run(self, R, hi, h0, RC, STO):
        h1=np.zeros([len(R)], dtype=np.float)
        h_tmp = (hi*1000 + R[0]/STO -hi*1000/RC)
        h1[0]=h_tmp
        for t in range(1,len(R)):
            h_tmp = h1[t-1] + R[t]/STO -h1[t-1]/RC
            h1[t]=h_tmp

        h=np.zeros([len(R)], dtype=np.float)
        for t in range(0,len(R)):
            h_tmp= (h1[t] + h0*1000)/1000
            h[t]=h_tmp

        return h

        del h1, h

#######################################################

class process:
    def __init__(self, MARM_ws, MF_ws, nrow, ncol, xllcorner, yllcorner, cellsizeMF, perlen, hnoflo):
        self.MARM_ws = MARM_ws
        self.MF_ws = MF_ws
        self.nrow= nrow
        self.ncol= ncol
        self.xllcorner= xllcorner
        self.yllcorner=yllcorner
        self.cellsizeMF=cellsizeMF
        self.perlen=perlen
        self.hnoflo=hnoflo


    def inputEsriAscii(self, grid_fn, datatype):

        grid_fn=os.path.join(self.MARM_ws,grid_fn)

        grid_out=np.zeros([self.nrow,self.ncol], dtype = datatype)

        grid_out=self.convASCIIraster2array(grid_fn,grid_out)
        return grid_out

    def convASCIIraster2array(self, filenameIN, arrayOUT):
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
        ncol_tmp = int(line[1])

        line = fin.readline().split()
        nrow_tmp = int(line[1])

        line = fin.readline().split()
        xllcorner_tmp = float(line[1])

        line = fin.readline().split()
        yllcorner_tmp = float(line[1])

        line = fin.readline().split()
        cellsizeEsriAscii = float(line[1])

        line = fin.readline().split()
        NODATA_value = float(line[1])

        # Process the file
        print "\nConverting %s to np.array..." % (filenameIN)
        for row in range(nrow_tmp):
        #   if (row % 100) == 0: print ".",
            for col in range(ncol_tmp):
                # Get new data if necessary
                if col == 0: line = fin.readline().split()
                if line[col] == NODATA_value:
                    arrayOUT[row,col] = self.hnoflo
                else:
                    arrayOUT[row,col] = line[col]

        # verify grid consistency between MODFLOW and ESRI ASCII
        if arrayOUT.shape[0] != self.nrow or arrayOUT.shape[1] != self.ncol or self.cellsizeMF != cellsizeEsriAscii:
            raise BaseException, '\nError in consistency between the MODFLOW grid and the input gridof the file %s.\nCheck the cell size and the number of rows, columns and cellsize' % filenameIN

        return arrayOUT

        fin.close()

    def inputTS(self,
                outputFILE_fn,
                SOILparam_fn, NMETEO, NVEG, NSOIL,
                inputDate_fn, inputZON_TS_RF_fn,
                inputZON_TS_PET_fn, inputZON_TS_RFe_fn,
                inputZON_TS_PE_fn, inputZON_TS_E0_fn
                ):   #IRR_fn

        # READ date of input files (RF and PET)
        inputDate_fn=os.path.join(self.MARM_ws, inputDate_fn)
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
        if len(inputDate)<>sum(self.perlen):
            raise ValueError, 'The number of time steps in MF (%i) is not the same as the number of days (%i) of the input data (RF and PET).\n' % (int(sum(self.perlen)), int(len(inputDate)))

        # soil parameters
        SOILparam_fn=os.path.join(self.MF_ws,SOILparam_fn)   #in MF folder for PEST
        if os.path.exists(SOILparam_fn):
            SOILparam = np.loadtxt(SOILparam_fn, dtype=float)
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % SOILparam_fn
        SOILzones=int(SOILparam[0])
        if SOILzones>NSOIL:
            print 'WARNING:\n' + str(SOILzones) + ' soil parameters groups in file [' + SOILparam_fn + ']\n Only ' + str(NSOIL) + ' PE time serie(s) found.'

        # Soils parameter initialisation
        Sm=np.zeros([SOILzones], dtype=np.float)
        Sfc=np.zeros([SOILzones], dtype=np.float)
        Sr=np.zeros([SOILzones], dtype=np.float)
        Si=np.zeros([SOILzones], dtype=np.float)
        Ks=np.zeros([SOILzones], dtype=np.float)
        SUSTm=np.zeros([SOILzones], dtype=np.float)
        n=np.zeros([SOILzones], dtype=np.int)
        f=np.zeros([SOILzones], dtype=np.float)

        # soil parameter definition for each soil zone
        for z in range(SOILzones):
            Sm[z]=SOILparam[2+9*z]
            Sfc[z]=SOILparam[3+9*z]
            Sr[z]=SOILparam[4+9*z]
            Si[z]=SOILparam[5+9*z]
            Ks[z]=SOILparam[6+9*z]
            SUSTm[z]=SOILparam[7+9*z]
            n[z]=SOILparam[8+9*z]
            f[z]=SOILparam[9+9*z]

        # READ input ESRI ASCII rasters vegetation
        gridVEGarea_fn=[]
        for v in range(NVEG):
            gridVEGarea_fn.append(os.path.join(self.MARM_ws,'inputVEG' + str(v+1)+'area.asc'))
        gridVEGarea=np.zeros([NVEG,self.nrow,self.ncol], dtype=float)
        for v in range(NVEG):
            gridtmp=np.zeros([self.nrow,self.ncol], dtype=float)
            gridVEGarea[v,:,:]=self.convASCIIraster2array(gridVEGarea_fn[v],gridtmp)

        # READ RF for each zone
        RF_fn=os.path.join(self.MARM_ws, inputZON_TS_RF_fn)
        if os.path.exists(RF_fn):
            RF = np.loadtxt(RF_fn)
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % RF_fn
        RFzonesTS=np.zeros([NMETEO,sum(self.perlen)], dtype=float)
        for i in range(NMETEO):
            for t in range(int(sum(self.perlen))):
                RFzonesTS[i,t]=RF[i*sum(self.perlen)+t]

        E0_fn=os.path.join(self.MARM_ws, inputZON_TS_E0_fn)
        if os.path.exists(E0_fn):
            E0 = np.loadtxt(E0_fn)
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % E0_fn
        E0zonesTS=np.zeros([NMETEO,sum(self.perlen)], dtype=float)
        for i in range(NMETEO):
            for t in range(int(sum(self.perlen))):
                E0zonesTS[i,t]=E0[i*sum(self.perlen)+t]

        ### READ IRR for each zone
        ##IRR_fn=os.path.join(self.MARM_ws,IRR_fn)
        ##if os.path.exists(IRR_fn):
        ##    IRR = np.loadtxt(IRR_fn)
        ##else:
        ##    raise ValueError, "\nThe file %s doesn't exist!!!" % IRR_fn
        ##IRRzonesnumb=int(IRR[0])
        ##IRRzones=np.zeros([IRRzonesnumb,sum(dis.self.perlen)], dtype=float)
        ##for i in range(IRRzonesnumb):
        ##    for j in range(sum(dis.self.perlen)):
        ##        IRRzones[i,j]=IRR[i*sum(dis.self.perlen)+j+1]


        # READ PET for each zone and each vegetation
        PET_fn = os.path.join(self.MARM_ws, inputZON_TS_PET_fn)
        if os.path.exists(PET_fn):
            PETtmp = np.loadtxt(PET_fn)
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % PET_fn
        PETvegzonesTS=np.zeros([NMETEO,NVEG,int(sum(self.perlen))], dtype=float)
        for z in range(NMETEO):
            for v in range(NVEG):
                for t in range(int(sum(self.perlen))):
                    PETvegzonesTS[z,v,t]=PETtmp[z*v*int(sum(self.perlen))+t]
                    #structure is [number of zones, number of vegetation type, time]
        PETvegzonesTS=np.asarray(PETvegzonesTS)

        # READ RFe for each zone and each vegetation
        RFe_fn = os.path.join(self.MARM_ws, inputZON_TS_RFe_fn)
        if os.path.exists(RFe_fn):
            RFetmp = np.loadtxt(RFe_fn)
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % RFe_fn
        RFevegzonesTS=np.zeros([NMETEO,NVEG,int(sum(self.perlen))], dtype=float)
        for z in range(NMETEO):
            for v in range(NVEG):
                for t in range(int(sum(self.perlen))):
                    RFevegzonesTS[z,v,t]=RFetmp[z*v*int(sum(self.perlen))+t]
                    #structure is [number of zones, number of vegetation type, time]
        RFevegzonesTS=np.asarray(RFevegzonesTS)


        # READ PE for each zone and each soil
        PE_fn = os.path.join(self.MARM_ws, inputZON_TS_PE_fn)
        if os.path.exists(PE_fn):
            PEtmp = np.loadtxt(PE_fn)
        else:
            raise ValueError, "\nThe file %s doesn't exist!!!" % PE_fn
        PEsoilzonesTS=np.zeros([NMETEO,NSOIL,int(sum(self.perlen))], dtype=float)
        for z in range(NMETEO):
            for v in range(NSOIL):
                for t in range(int(sum(self.perlen))):
                    PEsoilzonesTS[z,v,t]=PEtmp[z*v*int(sum(self.perlen))+t]
                    #structure is [number of zones, number of vegetation type, time]
        PEsoilzonesTS=np.asarray(PEsoilzonesTS)

        return Sm, Sr, Sfc, Si, Ks, SUSTm, n, f, gridVEGarea, RFzonesTS, E0zonesTS, PETvegzonesTS, RFevegzonesTS, PEsoilzonesTS, inputDate


    def inputObs(self, inputObs_fn, inputDate):
        '''
        observations cells for soil moisture and heads (will also compute SATFLOW)
        '''

        # read coordinates and SATFLOW parameters
        filenameIN = os.path.join(self.MARM_ws,inputObs_fn)
        if os.path.exists(filenameIN):
            fin = open(filenameIN, 'r')
            lines = fin.readlines()
            fin.close()
        else:
            raise BaseException, "The file %s doesn't exist!!!" % filenameIN

        # define a dictionnary of observatins,  format is: Name (key) x y i j hi h0 RC STO
        obs = {}
        for i in range(1,len(lines)):
            line = lines[i].split()
            name = line[0]
            x = float(line[1])
            y = float(line[2])
            try:
                hi  = float(line[3])
                h0  = float(line[4])
                RC  = float(line[5])
                STO = float(line[6])
            except:
                hi = h0 = RC = STO =  self.hnoflo
            # verify if coordinates are inside MODFLOW grid
            if x < self.xllcorner or \
               x > (self.xllcorner+self.ncol*self.cellsizeMF) or \
               y < self.yllcorner or \
               y > (self.yllcorner+self.nrow*self.cellsizeMF):
                   raise BaseException, 'The coordinates of the observation point %s are not inside the MODFLOW grid' % name
            # compute the cordinates in the MODFLOW grid
            #TODO use the PEST utilities for space extrapolation
            i = self.nrow - np.ceil((y-self.yllcorner)/self.cellsizeMF)
            j = np.ceil((x-self.xllcorner)/self.cellsizeMF) - 1
            obs[name] = {'x':x,'y':y,'i': i, 'j': j, 'hi':hi, 'h0':h0, 'RC':RC, 'STO':STO}

        outpathname=[]
        #  read obs time series
        for o in range(len(obs.keys())):
            outpathname.append(os.path.join(self.MARM_ws,'output'+obs.keys()[o]+'.txt'))
        obs_h=[]
        obs_sm=[]
        for o in range(len(obs.keys())):
            obsh_fn=os.path.join(self.MARM_ws,'inputObsHEADS_' + obs.keys()[o] +'.txt')
            obssm_fn=os.path.join(self.MARM_ws,'inputObsSM_' + obs.keys()[o] + '.txt')
            obs_h.append(self.verifObs(inputDate, obsh_fn))
            obs_sm.append(self.verifObs(inputDate, obssm_fn))

        return obs, outpathname, obs_h, obs_sm


    def verifObs(self, inputDate, filename):
        '''
        Import and process data and parameters
        '''

        try:
            #_________________Open file___________________________#
            obsOutput=np.zeros([len(inputDate)], dtype=float)
            if os.path.exists(filename):
                obsData=np.loadtxt(filename, dtype = str)
                obsDate = pylab.datestr2num(obsData[:,0])
                obsValue = obsData[:,1].astype(float)
                if obsDate[0]<inputDate[0]:
                    print 'WARNING, Obs data starting before RF and PET,\n these obs data will not be plotted correctly\n'
                if len(inputDate)<len(obsDate):
                    print 'WARNING, there is more Obs data than RF and PET data,\n these obs data will not be plotted correctly\n'
                j=0
                for i in range(len(inputDate)):
                    if j<len(obsDate):
                        if inputDate[i]==obsDate[j]:
                            obsOutput[i]=obsValue[j]
                            j=j+1
                        else:
                            obsOutput[i]=self.hnoflo
                    else:
                        obsOutput[i]=self.hnoflo
            else:
                for i in range(len(inputDate)):
                    obsOutput[i]=self.hnoflo

        except (ValueError, TypeError, KeyboardInterrupt, IOError), e:
            obsOutput = 0
            print e

        return obsOutput

        del obsOutput, obsData, obsDate, obsValue


    def outputEAgrd(self, outFile_fn, outFolder = []):

        if outFolder == []:
            outFolder = self.MARM_ws

        outFile=open(os.path.join(outFolder, outFile_fn), 'w')

        outFile=self.writeHeaderESRIraster(outFile)

        gridout=np.zeros([self.nrow,self.ncol], dtype=float)

        return outFile, gridout

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

    def ExportResults(self, inputDate, TS, hmeas, Smeas, outFileExport, outPESTheads, outPESTsm, obsname):
        """
        Export the processed data in a txt file
        INPUTS:      output fluxes time series and date
        OUTPUT:      output.txt
        """
        # Write the rest#
        for t in range(len(inputDate)):
            year='%4d'%pylab.num2date(inputDate[t]).year
            month='%02d'%pylab.num2date(inputDate[t]).month
            day='%02d'%pylab.num2date(inputDate[t]).day
            date=(day+"/"+month+"/"+year)
            if hmeas[t]!=self.hnoflo:
                hmeas_tmp = str(hmeas[t])
                outPESTheads.write(obsname.ljust(10,' ')+ date.ljust(14,' ')+ '00:00:00        '+ str(TS[15][t])+ '    \n')
            else:
                hmeas_tmp = str(self.hnoflo)
            if Smeas[t]!=self.hnoflo:
                Smeas_tmp = str(Smeas[t])
                outPESTsm.write((obsname+'SM').ljust(10,' ')+ date.ljust(14,' ')+ '00:00:00        '+ str(TS[6][t]) + '    \n')
            else:
                Smeas_tmp = str(self.hnoflo)
            # 'Date,RF,PET,RFe,Inter,ETu,Es,S,SUST,SUSTcount,Qs, Qscount,   Rp,R,hSATFLOW,hMF, hmeas,Smeas,SSATpart,FLOODcount, MB\n'
            # 0[results[i,j,iRF,:], 1results[i,j,iPET,:], 2results[i,j,iRFe,:],3results[i,j,iINTER,:], 4results[i,j,iETu,:], 5results[i,j,iETs,:],6results[i,j,iS,:], 7results[i,j,iSUST,:],8results[i,j,iPOND,:],9results[i,j,iQs,:],10results[i,j,iRunoff,:],11results[i,j,iFLOOD,:],12results[i,j,iRp,:], 13results[i,j,iR,:], 14h_satflow, 15heads[:,i,j,0], 16results[i,j,iSATpart,:],17results[i,j,iMB,:], 18results[i,j,iE0,:]]
            out_line =  pylab.num2date(inputDate[t]).isoformat()[:10], ',', str(TS[0][t]), ',', str(TS[18][t]), ',', str(TS[1][t]), ',', str(TS[2][t]), ',' ,str(TS[3][t]), ',',str(TS[4][t]), ',',str(TS[5][t]), ',',str(TS[6][t]), ',',str(TS[7][t]), ',', str(TS[8][t]), ',', str(TS[9][t]),  ',', str(TS[10][t]),  ',', str(TS[12][t]), ',', str(TS[13][t]), ',', str(TS[14][t]), ',', str(TS[15][t]), ',', hmeas_tmp,',', Smeas_tmp,',', str(TS[16][t]), ',', str(TS[11][t]),',', str(TS[17][t]),'\n'
            for l in out_line:
                outFileExport.write(l)

# EOF