from numpy import empty
from mbase import package

class mswtvdf(package):
    '''VDF pacakge'''
    def __init__(self, model, mtdnconc=1, mfnadvfd=1, nswtcpl=1, iwtable=1, \
                 densemin=1.000, densemax=1.025, \
                 dnscrit=1e-2, \
                 denseref=1.000, denseslp=.025, \
                 firstdt=0.001,
                 indense=0,
                 dense=1.000, extension='vdf'):
        package.__init__(self, model, extension, 'VDF', 37) # Call ancestor's init to set self.parent, extension, name and unit number
        nrow, ncol, nlay, nper = self.parent.mf.nrow_ncol_nlay_nper
        self.mtdnconc = mtdnconc 
        self.mfnadvfd = mfnadvfd
        self.nswtcpl = nswtcpl
        self.iwtable = iwtable
        self.densemin = densemin 
        self.densemax = densemax
        self.dnscrit = dnscrit
        self.denseref = denseref
        self.denseslp = denseslp
        self.firstdt = firstdt
        self.indense = indense
        self.dense = empty((nrow, ncol, nlay))
        self.assignarray( self.dense, dense )
        self.parent.add_package(self)
    def __repr__( self ):
        return 'Variable density flow package class'
    def write_file(self):
        nrow, ncol, nlay, nper = self.parent.mf.nrow_ncol_nlay_nper
        f_vdf = open(self.file_name[0], 'w')
        f_vdf.write('%10i%10i%10i%10i\n' % (self.mtdnconc, self.mfnadvfd, self.nswtcpl, self.iwtable))
        f_vdf.write('%10.4f%10.4f\n' % (self.densemin, self.densemax))
        if (self.nswtcpl > 1):
            f_vdf.write('%10f\n' % (self.dnscrit))
        f_vdf.write('%10.4f%10.4f\n' % (self.denseref, self.denseslp))
        f_vdf.write('%10f\n' % (self.firstdt))
        if (self.mtdnconc == 0):
            f_vdf.write('%10i\n' % (self.indense))
            if (self.indense > 0):
                for i in range(nlay):
                    comment = 'Density of layer ' + str(i + 1)
                    self.parent.write_array( f_vdf, self.dense[:,:,i], self.unit_number[0], True, 13, ncol, comment )
        f_vdf.close()
                
