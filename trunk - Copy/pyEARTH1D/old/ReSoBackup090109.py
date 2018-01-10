# -*- coding: cp1252 -*-
##=========================================================================================##
##===========================|------------------------------------|========================##
##===========================|               ReSo                 |========================##
##===========================|    Recharge and Soil Moisture      |========================##
##===========================|         from ADAS Station          |========================##
##===========================|------------------------------------|========================##
##=========================================================================================##
##=========================================================================================##

from pylab import *
from matplotlib.dates import MonthLocator, DateFormatter
from matplotlib.ticker import FormatStrFormatter
from sys import *
import os
import datetime
import wx

class ReSo: 
    """ ReSo from 'Recharge and Soil Moisture' computes on a daily temporal basis aquifer
    recharge using as input rainfall and potential evapotranspiration.
    It is a modular software, check below the several modules and their inputs.
    Equations are based on:
    Van der Lee J. and Gehrels, J. (1990)
    Modelling Aquifer Recharge – Introduction to the Lumped Parameter Model EARTH
    Free University of Amsterdam, The Netherlands
    """
##    def __init__(self, name, x, y):
##        self.name=name
##        self.x=x
##        self.y=y


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

    ##________________________________FUNCTIONS DEFINITION____________________________##

        #_________________________Actual evapotranspiration function__________________#
        
        def evp(s,pet):
            Se=(s-Sr)/(Sm-Sr)    # Percent. of soil saturation
            if s<Sr:
                return 0.0
            elif pet>((Sm-Sr)*D):
                return Se*((Sm-Sr)*D)
            else:
                return Se*pet
           
        #______________________________________________________________________________#

        #_________________________Percolation function_________________________________#

        def perc(s,D):            
            Sg=(s-Sfc)/(Sm-Sfc)  # Percent. of gravitational water
            if s<=Sfc:
                return (0.0)
            elif s>Sfc and s<=Sm:
                if (Ks*Sg) > (D*(s-Sfc)):
                    return (D*(s-Sfc))
                else:
                    return (Ks*Sg)
            elif s>Sm:
                if Ks>(D*(Sm-Sfc)):
                    return (D*(Sm-Sfc))
                else:
                    return (Ks)
        #____________________________________________________________________________#

        #_________________________Ponding and surface runoff function________________#

        def pond(s,D):
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
        #__________________________________________________________________________#

    ##________________________________END OF FUNCTIONS DEFINITION__________________##


    ##________________________________FUNCTION CORE________________________________##

       
        #__________________Pe: Excess Precipitation________________________________#
        Pe=[]
        for p in P:
                if p<=(MAXIL):
                        Pe.append(0.0)
                else:
                        Pe.append(p-(MAXIL))

        #__________________FIRST STEP - INITIAL CONDITIONS WITH Si_________________#

        S_tmp = Si*D + Pe[0]
#        print "##########  time step 0  ############################"
#        print "S_tmp init" +str(S_tmp)

        # SUST and Qs
        SUST=[]               # Surface storage inicialization
        Qs=[]                 # Surface runoff inicialization
        pond_tmp = pond(S_tmp/D,D)
        SUST.append(pond_tmp[0])
        Qs.append(pond_tmp[1])

        # FS - Rp
        Rp=[]                   #Deep percolation inicialization
        Rp.append(perc(S_tmp/D,D))

        # FS - ETa
        ETa=[]                  #Actual evapotranspiration inicialization
        ETa.append(evp(S_tmp/D,PET[0]))

        # FS - S
        S=[]                    #Soil moisture storage inicialization
        S_tmp=S_tmp-Rp[0]-ETa[0]-(SUST[0]+Qs[0])

        if (S_tmp/D)<(Sr):
            
            print "-----------------------"      
            print "Warning in time step " + str(len(S)) + ": S<Sr!"
            print "SM " + str(S_tmp/D) + "   Sr " + str(Sr)
            print "-----------------------"      

            S_tmp = Si*D + Pe[0]

            # SUST and Qs
            SUST=[]               # Surface storage inicialization
            Qs=[]                 # Surface runoff inicialization
            pond_tmp = pond(S_tmp/D,D)
            SUST.append(pond_tmp[0])
            Qs.append(pond_tmp[1])
            S_tmp=S_tmp - (SUST[0]+Qs[0])

            # FS - Rp
            Rp=[]                   #Deep percolation inicialization
            Rp.append(perc(S_tmp/D,D))
            S_tmp=S_tmp-Rp[0]

            # FS - ETa
            ETa=[]                  #Actual evapotranspiration inicialization
            ETa.append(evp(S_tmp/D,PET[0]))
            S_tmp=S_tmp-ETa[0]

        S.append(S_tmp/D)


        #___________SECOND STEP - PROCESSING THE WHOLE DATA SET_________#

        for t in range(1,len(P)):               # t: current time step

            S_tmp = S[t-1]*D + Pe[t] + SUST[t-1]
#            print "######## time step " + str(t) + "##############################"
#            print "S_tmp init" +str(S_tmp)
            # SUST and Qs
            pond_tmp = pond(S_tmp/D,D)
            SUST.append(pond_tmp[0])
            Qs.append(pond_tmp[1])

            # SS - Rp
            Rp.append(perc(S_tmp/D,D))

            # SS - ETa
            ETa.append(evp(S_tmp/D,PET[t]))

            # SS - S
            S_tmp=S_tmp-Rp[t]-ETa[t]-(SUST[t]+Qs[t])

            if (S_tmp/D)<(Sr):
                print "-----------------------"                
                print "Warning in time step " + str(len(S)) + ": S<Sr!"
                print "SM " + str(S_tmp/D) + "  Sr " + str(Sr)
                print "Rp try1 " + str(Rp[t])
                print "ETa try1 " + str(ETa[t])
                print "S_tmp " + str(S_tmp)
                print "-----------------------"
                
                S_tmp = S[t-1]*D + Pe[t] + SUST[t-1]
                
                # SUST and Qs
                pond_tmp = pond(S_tmp/D,D)
                SUST[t] = pond_tmp[0]
                Qs[t] = pond_tmp[1]
                S_tmp = S_tmp - (SUST[t] + Qs[t])

                # FS - Rp
                Rp[t]=perc(S_tmp/D,D)
                S_tmp=S_tmp-Rp[t]

                # FS - ETa
                ETa[t] = evp(S_tmp/D,PET[t])
                S_tmp = S_tmp - ETa[t]

                print "Rp try2 " + str(Rp[t])
                print "ETa try2 " + str(ETa[t])
                print "S_tmp try2 " + str(S_tmp)
                print "-----------------------"                  

            S.append(S_tmp/D)

        return P, PET, Pe, SUST, Qs, ETa, S, Rp

        del P, PET, Pe, SUST, Qs, ETa, S, Rp

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
                    hi          Interception looses
                    h0          Max. Soil moiture storage 
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
            h_tmp= (h1[t] + h0*1000 )/1000
            h.append(h_tmp)

        return h

        del h_tmp, h1, h

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

        InputTest=0
        try:
            #_________________Open input file for date, P ad PET___________________________#
            if os.path.exists(input_path):
                input_ar=load(input_path, converters={0:datestr2num})
                DateInput = input_ar[:,0]
                for i in range(1,input_ar.shape[0]):
                    #__________________Check date consistency________________#
                    difDay=DateInput[i]-DateInput[i-1]
                    if (difDay !=1.0):
                        print 'DifDay = ' + str(difDay)
                        print 'Dates are not sequencial, read manual and check your daily time step!\nError in date ' + str(gmtime(DateInput[i]))
                        sys.exit()
                P = input_ar[:,1]
                PET = input_ar[:,2]
                #_________________Open measured piezometric levels___________________________#
                hmeas=[]
                if os.path.exists(hmeas_path):
                    hmeas_ar=load(hmeas_path, converters={0:datestr2num})
                    Datehmeas = hmeas_ar[:,0]
                    if Datehmeas[0]<DateInput[0]:
                        print 'WARNING, hmeas date is before P and PET date,\n data will not be plotted correctly'
                    if len(DateInput)<len(Datehmeas):
                        print 'WARNING, there is more hmeas data than P and PET data,\n data will not be plotted correctly'
                    hmeas_tmp = hmeas_ar[:,1]
                    j=0
                    for i in range(len(DateInput)):
                        if j<len(Datehmeas):
                            if DateInput[i]==Datehmeas[j]:
                                hmeas.append(hmeas_tmp[j])
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
                    Smeas_ar=load(Smeas_path, converters={0:datestr2num})
                    DateSmeas=Smeas_ar[:,0]
                    if DateSmeas[0]<DateInput[0]:
                        print 'WARNING, Smeas before P and PET,\n data will not be plotted correctly'
                    if len(DateInput)<len(DateSmeas):
                        print 'WARNING, there is more Smeas data than P and PET data,\n data will not be plotted correctly'
                    Smeas_tmp = Smeas_ar[:,1]
                    j=0
                    for i in range(len(DateInput)):
                        if j<len(DateSmeas):
                            if DateInput[i]==DateSmeas[j]:
                                Smeas.append(Smeas_tmp[j])
                                j=j+1
                            else:
                                Smeas.append(-999)
                        else:
                            Smeas.append(-999)
                else:
                    for i in range(len(DateInput)):
                        Smeas.append(-999)
            else:
                InputTest=-1
                DateInput, P , PET, hmeas, Smeas = 0,0,0,0,0          

        except (ValueError, TypeError, KeyboardInterrupt, IOError), e:
            InputTest=-2
            DateInput, P , PET, hmeas, Smeas = 0,0,0,0,0
            print e

        return InputTest, DateInput, P , PET, hmeas, Smeas

        del InputTest, DateInput, P , PET, hmeas, Smeas


##=========================================================================================##
##=============================| MAIN LOOP |===============================================##
##=========================================================================================##

def runReSo(strDateFormat, param, input_path, hmeas_path, Smeas_path):

    dataset=ReSo()

    InputTest, DateInput, P , PET, hmeas, Smeas = dataset.processInput(strDateFormat, input_path, hmeas_path, Smeas_path)    
    #__________________Check input file____________________________________________#
    if InputTest<0:
        DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas = 0,0,0,0,0,0,0,0,0,0,0,0,0
        return InputTest, DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas
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
    except ValueError:
        InputTest=-3
        DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas = 0,0,0,0,0,0,0,0,0,0,0,0,0
        return InputTest, DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas
        del InputTest, DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas
        sys.exit() 


    #__________________Check parameters______________________________________________#
    if Sm>Sfc>Sr and Sm>Si>=Sr:
        pass
    else:
        DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas = 0,0,0,0,0,0,0,0,0,0,0,0,0
        return InputTest, DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas
        del InputTest, DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas
        sys.exit() 

    # cal functions for reservoirs calculations
    P, PET, Pe, SUST, Qs, ETa, S, Rp = dataset.SOMOS(MAXIL, Sm, Sfc, Sr, Si, D, Ks, SUSTm, P, PET)
    R = dataset.LINRES(Rp, n, f)
    h = dataset.SATFLOW(R, hi, h0, RC, STO)

    


     
    return InputTest, DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas
    del InputTest, DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas


def callExportResults(DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, Output_path):
    dataset=ReSo()  
    dataset.ExportResults(DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas, Output_path)
    del DateInput, P, PET, Pe, SUST, Qs, ETa, S, Rp, R, h, hmeas, Smeas
