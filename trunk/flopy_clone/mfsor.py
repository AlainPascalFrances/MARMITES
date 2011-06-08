from mbase import package

class mfsor(package):
    '''Sor package
    '''
    def __init__(self, model, mxiter=200,
                 accl=1, hclose=1e-5, iprsor=0, extension='sor', unitnumber=26):
        package.__init__(self, model, extension, 'sor', unitnumber) # Call ancestor's init to set self.parent, extension, name and unit number
        self.url = 'sor.htm'
        self.mxiter = mxiter
        self.accl= accl
        self.hclose = hclose
        self.iprsor = iprsor
        self.parent.add_package(self)
    def __repr__( self ):
        return 'Slice-successive overrelaxation package class'
    def write_file(self):
        # Open file for writing
        f_sor = open(self.fn_path, 'w')
        f_sor.write('%10i\n' % (self.mxiter))
        f_sor.write('%10f%10f%10i\n' % (self.accl,self.hclose,self.iprsor))
        f_sor.close()

