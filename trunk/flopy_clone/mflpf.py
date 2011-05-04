import numpy as np
from mbase import package

class mflpf(package):
    'Layer-property flow package class'
    def __init__(self, model, laytyp=1, layavg=0, chani=1.0, layvka=0, laywet=0, hdry=-1E+30, iwdflg=0, wetfct=0.1, iwetit=1, ihdwet=0, \
                 hk=1.0, hani=1.0, vka=1.0, ss=1e-5, sy=0.15, vkcb=0.0, wetdry=-0.01, extension='lpf', unitnumber=15):
        package.__init__(self, model, extension, 'LPF', unitnumber) # Call ancestor's init to set self.parent, extension, name and unit number
        self.url = 'lpf.htm'
        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        # First create arrays so that they have the correct size
        self.laytyp = np.empty(nlay, dtype='int32') # Specifies both the layer type (LAYCON) and the method of computing interblock conductance
        self.layavg = np.ones(nlay, dtype='int32') # Interblock transmissivity flag for each layer
        self.chani = np.ones(nlay) # Horizontal anisotropy flag for each layer
        self.layvka = np.ones(nlay, dtype='int32') # vertical hydraulic conductivity flag for each layer
        self.laywet = np.ones(nlay, dtype='int32') # wet dry flag for each layer
        #self.laycbd = np.ones(nlay, dtype='int32') # confining bed flag for each layer
        # Set values of all parameters
        self.assignarray( self.laytyp, laytyp )
        self.assignarray( self.layavg, layavg )
        self.assignarray( self.chani, chani )
        self.assignarray( self.layvka, layvka )
        self.assignarray( self.laywet, laywet )
        self.ilpfcb = 50 # Unit number for file with cell-by-cell flow terms
        self.nplpf = 0 # number of LPF parameters
        self.hdry = hdry # Head in cells that are converted to dry during a simulation
        self.iwdflg = iwdflg # Flag that determines if the wetting capability is active
        self.wetfct = wetfct # Factor that is included in the calculation of the head when a cell is converted from dry to wet
        self.iwetit = iwetit # Iteration interval for attempting to wet cells
        self.ihdwet = ihdwet # Flag that determines which equation is used to define the initial head at cells that become wet
        self.hk = np.empty((nrow, ncol, nlay))
        self.hani = np.empty((nrow, ncol, nlay))
        self.vka = np.empty((nrow, ncol, nlay))
        self.ss = np.empty((nrow, ncol, nlay))
        self.sy = np.empty((nrow, ncol, nlay))
        self.vkcb = np.empty((nrow, ncol, nlay))
        self.wetdry = np.empty((nrow, ncol, nlay))
        self.assignarray( self.hk, hk )
        self.assignarray( self.hani, hani )
        self.assignarray( self.vka, vka )
        self.assignarray( self.ss, ss )
        self.assignarray( self.sy, sy )
        self.assignarray( self.vkcb, vkcb )
        self.assignarray( self.wetdry, wetdry )
        self.parent.add_package(self)
    def __repr__( self ):
        return 'Layer-property flow package class'
    def write_file(self):
        nrow, ncol, nlay, nper = self.parent.nrow_ncol_nlay_nper
        # Open file for writing
        f_lpf = open(self.fn_path, 'w')
        # Item 0: IBCFCB, HDRY, NPLPF
        f_lpf.write('%10d%10.1e%10d\n' % (self.ilpfcb, self.hdry, self.nplpf))
        # LAYTYP array
        self.parent.write_array(f_lpf, self.laytyp, self.unit_number[0], False, 2, -40, 'LAYTYP(): Layer type of layers')
        # LAYAVG array
        self.parent.write_array(f_lpf, self.layavg, self.unit_number[0], False, 2, -40, 'LAYAVG(): Layer average of layers')
        # CHANI array
        self.parent.write_array(f_lpf, self.chani, self.unit_number[0], False, 2, -40, 'CHANI(): Horizontal anisotropy flag of layers')
        # LAYVKA array
        self.parent.write_array(f_lpf, self.layvka, self.unit_number[0], False, 2, -40, 'LAYVKA(): Vertical hydraulic conductivity flag of layers')
        # LAYWET array
        self.parent.write_array(f_lpf, self.laywet, self.unit_number[0], False, 2, -40, 'LAYWET(): Wet-Dry flag of layers')
        # Item 0: WETFCT, IWETIT, IHDWET
        iwetdry = self.laywet.sum()
        if iwetdry > 0:
        	f_lpf.write('%10f%10d%10d\n' % (self.wetfct, self.iwetit, self.ihdwet))
        transient = not self.parent.get_package('DIS').steady.all()
        for i in range(nlay):
        	comment = 'HK() = Horizontal hydraulic conductivity of layer ' + str(i + 1)
        	self.parent.write_array( f_lpf, self.hk[:,:,i], self.unit_number[0], True, 13, ncol, comment )
        	if self.chani[i] < 1:
	        	comment = 'HANI() = Ratio of horizontal hydraulic of columns to rows of layer ' + str(i + 1)
	        	self.parent.write_array( f_lpf, self.hani[:,:,i], self.unit_number[0], True, 13, ncol, comment )
        	comment = 'VKA() = Vertical hydraulic conductivity of layer ' + str(i + 1)
        	self.parent.write_array( f_lpf, self.vka[:,:,i], self.unit_number[0], True, 13, ncol, comment )
        	if transient == True:
        		comment = 'Ss() = Specific storage coefficient of layer ' + str(i + 1)
        		self.parent.write_array( f_lpf, self.ss[:,:,i], self.unit_number[0], True, 13, ncol, comment )
        		if self.laytyp[i] !=0:
        			comment = 'Sy() = Specific yield of layer ' + str(i + 1)
        			self.parent.write_array( f_lpf, self.sy[:,:,i], self.unit_number[0], True, 13, ncol, comment )
        	if self.parent.get_package('DIS').laycbd[i] > 0:
        		comment = 'VKCB() = Vertical hydraulic conductivity of quasi-three-dimensional confining bed of layer ' + str(i + 1)
        		self.parent.write_array( f_lpf, self.vkcb[:,:,i], self.unit_number[0], True, 13, ncol, comment )
        	if (self.laywet[i] != 0 and self.laytyp != 0):
        		comment = 'WETDRY() = Wetting threshold of layer ' + str(i + 1)
        		self.parent.write_array( f_lpf, self.wetdry[:,:,i], self.unit_number[0], True, 13, ncol, comment )
        f_lpf.close()

