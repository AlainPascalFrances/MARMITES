# -*- coding: utf-8 -*-
##=========================================================================================##
##===========================|------------------------------------|========================##
##===========================|               ReSo                 |========================##
##===========================|    Recharge and Soil Moisture      |========================##
##===========================|         from ADAS Station          |========================##
##===========================|------------------------------------|========================##
##=========================================================================================##
##=========================================================================================##

import numpy
from pylab import *
from matplotlib.dates import MonthLocator, DateFormatter
from matplotlib.ticker import FormatStrFormatter
from sys import *
import os

class ReSo:
    """ ReSo from 'Recharge and Soil Moisture' computes on a daily temporal basis aquifer
    recharge using as input rainfall and potential evapotranspiration.
    It is a modular software, check below the several modules and their inputs.
    Equations are based on:
    Van der Lee J. and Gehrels, J. (1990)
    Modelling Aquifer Recharge – Introduction to the Lumped Parameter Model EARTH
    Free University of Amsterdam, The Netherlands
    """

##=========================================================================================##
##===========================|-----------------|===========================================##
##===========================| SOMOS FUNCTION  |===========================================##
##===========================|-----------------|===========================================##
##=========================================================================================##

    def SOMOS(self, MAXIL, Sm, Sfc, Sr, Si, D, Ks, SUSTm, P, PET):

        """
        SOMOS: SOil MOisture Storage
        Calculate water balance in soil
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
                    P           Daily rainfall
                    PET         Daily evapotranspiration
        OUTPUTS
                Pe              Daily Excess rainfall
                ETa             Daily evapotranspiration
                S               Daily soil moisture
                Rp              Daily percolation
                SUST            Daily ponding
                Qs              Daily runoff

        Provide daily values for the LINRES module
        ______________________________________________________________________________
        ______________________________________________________________________________
        """

    ##________________________________FUNCTION CORE________________________________##


        #__________________FIRST STEP - INITIAL CONDITIONS WITH Si_________________#

        Pe=[]   # Precipitation excess inicialization
        SUST=[] # Surface storage inicialization
        Qs=[]   # Surface runoff inicialization
        Rp=[]   # Deep percolation inicialization
        ETa=[]  # Actual evapotranspiration inicialization
        S=[]    # Soil moisture storage inicialization
        MB = 0.0


        for t in range(len(P)):

            if t>0:
                Sprev = S[t-1]
                SUSTprev = SUST[t-1]
            else:
                Sprev = Si
                SUSTprev = 0

            if P[t]>MAXIL:
                Pe.append(P[t]-MAXIL)
            else:
                Pe.append(0.0)

            S_tmp = Sprev*D + Pe[t] + SUSTprev

            # SUST and Qs
            pond_tmp = self.pond(S_tmp/D,D,Sm,SUSTm)
            SUST.append(pond_tmp[0])
            Qs.append(pond_tmp[1])
            S_tmp=S_tmp-(SUST[t]+Qs[t])

            # FS - Rp
            Rp.append(self.perc(S_tmp/D,D,Sm,Sfc,Ks))
            S_tmp=S_tmp-Rp[t]

            # FS - ETa
            ETa.append(self.evp(S_tmp/D,PET[t],D,Sm,Sr))
            S_tmp=S_tmp-ETa[t]

            S.append(S_tmp/D)

            MB += P[t]+SUSTprev-(P[t]-Pe[t])-Rp[t]-Qs[t]-SUST[t]-ETa[t]-((S[t]-Sprev)*D)
            
        return P, PET, Pe, SUST, Qs, ETa, S, Rp, MB

        del P, PET, Pe, SUST, Qs, ETa, S, Rp

##=========================================================================================##
##________________________________SOMOS FUNCTIONS DEFINITION_______________________________##

    #_________________________Ponding and surface runoff function________________#

    def pond(self,s,D,Sm,SUSTm):
        if s>Sm:
            sust_tmp =D*(s-Sm)
            if sust_tmp > SUSTm:
                qs_tmp =(sust_tmp-SUSTm)
                sust_tmp = SUSTm
            else:
                qs_tmp =(0.0)
        else:
            sust_tmp =(0.0)
            qs_tmp =(0.0)
        return sust_tmp, qs_tmp

    #_________________________Actual evapotranspiration function__________________#

    def evp(self, s,pet,D,Sm,Sr):
        Se=(s-Sr)/(Sm-Sr)    # Percent. of soil saturation
        if s<Sr:
            return 0.0
        elif pet>((Sm-Sr)*D):
            return Se*((Sm-Sr)*D)
        else:
            return Se*pet

    #_________________________Percolation function_________________________________#

    def perc(self,s,D,Sm,Sfc,Ks):
        Sg=(s-Sfc)/(Sm-Sfc)  # Percent. of gravitational water
        if s<=Sfc:
            return (0.0)
        elif s<=Sm:
            if (Ks*Sg) > (D*(s-Sfc)):
                return (D*(s-Sfc))
            else:
                return (Ks*Sg)
        else:
            if Ks>(D*(Sm-Sfc)):
                return (D*(Sm-Sfc))
            else:
                return (Ks)
    #____________________________________________________________________________#

#________________________________END OF FUNCTIONS DEFINITION__________________##

##=========================================================================================##
##===========================|-----------------|===========================================##
##===========================| LINRES FUNCTION |===========================================##
##===========================|-----------------|===========================================##
##=========================================================================================##
    def LINRES(self, Rp, n, f):
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

        Y=zeros((len(Rp), n+1), float)

        # Initialization of the first line of the array
        Y[0,0] = (1+f)*Rp[0]/f

        # Main loop
        for t in range(1,len(Rp)):
            for i in range(0,n+1):
                if i==0:
                    Y[t,i] = (1+f)*Rp[t]/f
                else:
                    for j in range(0,i+1):
                        Y[t,i] = Y[t,i]+((1+f)**(-j))*Y[t-1,i-j]
                    Y[t,i] = Y[t,i] * f/(1+f)
        R = Y[:,n]

        return R

        del R, Y

##=========================================================================================##
##===========================|------------------|==========================================##
##===========================| SATFLOW FUNCTION |==========================================##
##===========================|------------------|==========================================##
##=========================================================================================##

    def SATFLOW(self, R, hi, h0, RC, STO):
        """
        SATFLOW: SATurated FLOW
        Calculate water level fluctuations in function of recharge (R processed in LINRES)
        _______________________________________________________________________________

        INPUTS
                PARAMETERS
                    hi
                    h0          Base level
                    RC          Saturated recession constant
                    STO         Aquifer Storage capacity
                STATE VARIABLES
                    R           Daily recharge
        OUTPUTS
                h               Daily water level
        ______________________________________________________________________________
        ______________________________________________________________________________
        """

        h1=[]
        h_tmp = (hi*1000 + R[0]/STO -hi*1000/RC)
        h1.append(h_tmp)
        for t in range(1,len(R)):
            h_tmp = h1[t-1] + R[t]/STO -h1[t-1]/RC
            h1.append(h_tmp)

        h=[]
        for t in range(0,len(R)):
            h_tmp= (h1[t] + h0*1000)/1000
            h.append(h_tmp)

        return h

        del h1, h

##=========================================================================================##
##====================| Correct flooding effect ( h> piezometer elevation |================##
##=========================================================================================##

    def Elevation(self, STO, RC, h0, SUSTm, Elev, PET, R, SUST, Qs, h):

        for i in range(len(h)):
            if h[i]>Elev:
                h[i]=Elev
                Rexcdt_tmp=R[i]
                R[i]=STO*(1000*(Elev-h0)+1000*(h[i-1]-h0)*(1-RC)/RC)
                Rexcdt=Rexcdt_tmp-R[i]+SUST[i-1]-PET[i]
                if Rexcdt > SUSTm:
                    qs_tmp =(Rexcdt-SUSTm)
                    sust_tmp = SUSTm
                else:
                    qs_tmp =(0.0)
                    if Rexcdt>0:
                        sust_tmp=Rexcdt
                    else:
                        sust_tmp=0
                SUST[i]=sust_tmp
                Qs[i]=qs_tmp

        return R, SUST, Qs, h

        del R, SUST, Qs, h

##=========================================================================================##
##===========================| EXPORT RESULTS |============================================##
##=========================================================================================##

    def ExportResults(self, DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, Output_path):
        """
        OUT: OUTput data
        Export the processed data in a txt file
        _______________________________________________________________________________

        INPUTS
                    The same as the GRAPH module
        OUTPUT
                    output.txt      txt output file
        _______________________________________________________________________________
        _______________________________________________________________________________
        """
        #__________________Write first output in a txt file______________________#


        outFile=open(Output_path, 'w')
        outFile.write('Date,P,PET,Pe,ETa,S,SUST,Qs,Rp,R,h,hmeas,Smeas\n')
        #__________________         Write the rest        ______________________#

        for t in range(0,len(S)):
            if hmeas[t]!=-999:
                hmeas_tmp = str(hmeas[t])
            else:
                hmeas_tmp = ""
            if Smeas[t]!=-999:
                Smeas_tmp = str(Smeas[t])
            else:
                Smeas_tmp = ""
            out_line = num2date(DateInput[t]).isoformat()[:10], ',', str(P[t]), ',', str(PET[t]), ',', str(Pe[t]), ',' ,str(ETa[t]), ',',str(S[t]), ',', str(SUST[t]), ',',str(Qs[t]), ',', str(Rp[t]),  ',', str(R[t]),  ',', str(h[t]), ',',  hmeas_tmp, ',', Smeas_tmp,'\n'
            for l in out_line:
                outFile.write(l)
        outFile.close()

##=========================================================================================##
##===========================| Import and process data and parameters |====================##
##=========================================================================================##

    def processInput(self, strDateFormat, input_path, hmeas_path, Smeas_path):

        ErrorType = 999
        try:
            #_________________Open input file for date, P ad PET___________________________#
            if os.path.exists(input_path):
                input_ar=numpy.loadtxt(input_path, converters={0:date2num})
                DateInput = input_ar[:,0]
                for i in range(1,input_ar.shape[0]):
                    #__________________Check date consistency________________#
                    difDay=DateInput[i]-DateInput[i-1]
                    if (difDay !=1.0):
                        sys.exit()
                P = input_ar[:,1]
                PET = input_ar[:,2]
                #_________________Open measured piezometric levels___________________________#
                hmeas=[]
                if os.path.exists(hmeas_path):
                    hmeas_ar=numpy.loadtxt(hmeas_path, converters={0:date2num})
                    Datehmeas = hmeas_ar[:,0]
                    if Datehmeas[0]<DateInput[0]:
                        ErrorType = 101;
                    if len(DateInput)<len(Datehmeas):
                        ErrorType=102
                    hmeas_tmp = hmeas_ar[:,1]
                    j=0
                    for i in range(len(DateInput)):
                        if j<len(Datehmeas):
                            if DateInput[i]==Datehmeas[j]:
                                if isinstance(float(hmeas_tmp[j]), float):
                                    hmeas.append(hmeas_tmp[j])
                                else:
                                    hmeas.append(-999)
                                j=j+1
                            else:
                                hmeas.append(-999)
                        else:
                            hmeas.append(-999)
                else:
                    for i in range(len(DateInput)):
                        hmeas.append(-999)
                #_________________Open measured soil moisture___________________________#
                Smeas=[]
                if os.path.exists(Smeas_path):
                    Smeas_ar=numpy.loadtxt(Smeas_path, converters={0:date2num})
                    DateSmeas=Smeas_ar[:,0]
                    if DateSmeas[0]<DateInput[0]:
                        ErrorType = 103
                    if len(DateInput)<len(DateSmeas):
                        ErrorType = 104
                    Smeas_tmp = Smeas_ar[:,1]
                    j=0
                    for i in range(len(DateInput)):
                        if j<len(DateSmeas):
                            if DateInput[i]==DateSmeas[j]:
                               if isinstance(float(Smeas_tmp[j]), float):
                                    Smeas.append(Smeas_tmp[j])
                               else:
                                    Smeas.append(-999)
                               j=j+1
                            else:
                                Smeas.append(-999)
                        else:
                            Smeas.append(-999)
                else:
                    for i in range(len(DateInput)):
                        Smeas.append(-999)
            else:
                ErrorType = -1
                DateInput, P , PET, hmeas, Smeas = 0,0,0,0,0

        except (ValueError, TypeError, KeyboardInterrupt) as e:
            DateInput, P , PET, hmeas, Smeas = 0,0,0,0,0
            ErrorType = e
        except (SystemExit):
            DateInput, P , PET, hmeas, Smeas = 0,0,0,0,0
            ErrorType = 100

        return DateInput, P , PET, hmeas, Smeas, ErrorType

        del DateInput, P , PET, hmeas, Smeas

##=========================================================================================##
##=============================| MAIN LOOP |===============================================##
##=========================================================================================##

def runReSo(strDateFormat, param, input_path, hmeas_path, Smeas_path):

    dataset=ReSo()

    DateInput, P , PET, hmeas, Smeas, ErrorType = dataset.processInput(strDateFormat, input_path, hmeas_path, Smeas_path)
    #__________________Check input file____________________________________________#
    if ErrorType == -1 or ErrorType==100 or isinstance(ErrorType, str):
        DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas = 0,0,0,0,0,0,0,0,0,0,0,0,0
        return DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, ErrorType
        sys.exit()

    # get input data from GUI
    try:
        MAXIL=float(param[0])
        Sm=float(param[1])
        Sfc=float(param[2])
        Sr=float(param[3])
        Si=float(param[4])
        Ks=float(param[5])
        D=float(param[6])
        SUSTm=float(param[7])
        n=int(float(param[8]))
        f=float(param[9])
        RC=float(param[10])
        STO=float(param[11])
        hi=float(param[12])
        h0=float(param[13])
        Elev=float(param[14])
    except (ValueError, TypeError):
        ErrorType = -3
        DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas = 0,0,0,0,0,0,0,0,0,0,0,0,0
        return DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, ErrorType
        del DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, ErrorType
        sys.exit()


    #__________________Check parameters______________________________________________#
    if Sm>Sfc>Sr and Sm>Si>=Sr:
        pass
    else:
        ErrorType = 0
        DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas = 0,0,0,0,0,0,0,0,0,0,0,0,0
        return DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, ErrorType
        del DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, ErrorType
        sys.exit()

    # cal functions for reservoirs calculations
    print( "\nData verified succesfully!")
    P, PET, Pe, SUST, Qs, ETa, S, Rp, MB = dataset.SOMOS(MAXIL, Sm, Sfc, Sr, Si, D, Ks, SUSTm, P, PET)
    print("\nSOMOS done!")
    print("MASS BALANCE: %.2f" % MB)
    R = dataset.LINRES(Rp, n, f)
    print("\nLINRES done!")
    h = dataset.SATFLOW(R, hi, h0, RC, STO)
    print("\nSATFLOW done!")

    ## check if h> elevation of piezometer

    R, SUST, Qs, h = dataset.Elevation(STO, RC, h0, SUSTm, Elev, PET, R, SUST, Qs, h)

    return DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, ErrorType

    del DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, ErrorType


def callExportResults(DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, Output_path):
    dataset=ReSo()
    dataset.ExportResults(DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, Output_path)
    del DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas