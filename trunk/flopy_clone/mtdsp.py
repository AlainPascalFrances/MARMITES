from numpy import empty
from mbase import package

class mtdsp(package):
    'Dispersion package class'
    def __init__(self, model, al=0.01, trpt=0.1, trpv=0.01, dmcoef=1e-9, extension='dsp'):
        package.__init__(self, model, extension, 'DSP', 33) # Call ancestor's init to set self.parent, extension, name and unit number
        nrow, ncol, nlay, nper = self.parent.mf.nrow_ncol_nlay_nper
        # First create arrays so that they have the correct size
        self.al = empty((nrow, ncol, nlay))
        self.trpt = empty((nlay))
        self.trpv = empty((nlay))
        self.dmcoef = empty((nlay))
        # Set values of all parameters
        self.assignarray( self.al, al )
        self.assignarray( self.trpt, trpt )
        self.assignarray( self.trpv, trpv )
        self.assignarray( self.dmcoef, dmcoef )
        self.parent.add_package(self)
    def __repr__( self ):
        return 'Dispersion package class'
    def write_file(self):
        nrow, ncol, nlay, nper = self.parent.mf.nrow_ncol_nlay_nper
        # Open file for writing
        f_dsp = open(self.file_name[0], 'w')
        self.parent.write_array( f_dsp, self.al, self.unit_number[0], True, 13, -ncol, 'Longitudinal dispersivity for Layer')
        self.parent.write_array( f_dsp, self.trpt, self.unit_number[0], True, 13, -1, 'TRPT=(horizontal transverse dispersivity) / (Longitudinal dispersivity)')
        self.parent.write_array( f_dsp, self.trpv, self.unit_number[0], True, 13, -1, 'TRPV=(vertical transverse dispersivity) / (Longitudinal dispersivity)')
        self.parent.write_array( f_dsp, self.dmcoef, self.unit_number[0], True, 13, -1, 'Effective molecular diffusion coefficient')
        f_dsp.close()

