from mbase import package

class mfpcg(package):
    '''Pcg pacakge
    Only programmed to work with the default values; may need work for other options'''
    def __init__(self, model, mxiter=50, iter1=30, npcond=1, \
                 hclose=1e-5, rclose=1e-5, relax=1.0, nbpol=0, iprpcg=0, mutpcg=3, damp=1.0, extension='pcg', unitnumber=27):
        package.__init__(self, model, extension, 'PCG', unitnumber) # Call ancestor's init to set self.parent, extension, name and unit number
        self.url = 'pcg.htm'
        self.mxiter = mxiter
        self.iter1 = iter1
        self.npcond = npcond
        self.hclose = hclose
        self.rclose = rclose
        self.relax = relax
        self.nbpol = nbpol
        self.iprpcg = iprpcg
        self.mutpcg = mutpcg
        self.damp = damp
        self.parent.add_package(self)
    def __repr__( self ):
        return 'Preconditioned conjugate gradient package class'
    def write_file(self):
        # Open file for writing
        f_pcg = open(self.fn_path, 'w')
        f_pcg.write('%10i%10i%10i\n' % (self.mxiter,self.iter1,self.npcond))
        f_pcg.write('%10f%10f%10f%10i%10i%10i%10f\n' % (self.hclose,self.rclose,self.relax,self.nbpol,self.iprpcg,self.mutpcg,self.damp))
        f_pcg.close()

