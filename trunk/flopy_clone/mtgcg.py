from mbase import package

class mtgcg(package):
    'Generalized Conjugate Gradient solver package class'
    def __init__(self, model, mxiter=1, iter1=50, isolve=3, ncrs=0, \
                 accl=1, cclose=1e-5, iprgcg=0, extension='gcg'):
        package.__init__(self, model, extension, 'GCG', 35) # Call ancestor's init to set self.parent, extension, name and unit number
        self.mxiter = mxiter
        self.iter1 = iter1
        self.isolve = isolve
        self.ncrs = ncrs
        self.accl = accl
        self.cclose = cclose
        self.iprgcg = iprgcg
        self.parent.add_package(self)
    def __repr__( self ):
        return 'Generalized Conjugate Gradient solver package class'
    def write_file(self):
        # Open file for writing
        f_gcg = open(self.file_name[0], 'w')
        f_gcg.write('%10d%10d%10d%10d\n' % (self.mxiter, self.iter1, self.isolve, self.ncrs))
        f_gcg.write('%10d%10f%10d\n' % (self.accl, self.cclose, self.iprgcg))
        f_gcg.close()

