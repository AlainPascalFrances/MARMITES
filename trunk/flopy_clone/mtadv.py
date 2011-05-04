from mbase import package

class mtadv(package):
    'MT3DMS ADV pacakge'
    def __init__(self, model, mixelm=3, percel=0.75, mxpart=800000, nadvfd=1, \
                 itrack=3, wd=0.5, \
                 dceps=1e-5, nplane=2, npl=10, nph=40, npmin=5, npmax=80, \
                 nlsink=0, npsink=15,
                 dchmoc=0.0001, extension='adv'):
        package.__init__(self, model, extension, 'ADV', 32) # Call ancestor's init to set self.parent, extension, name and unit number
        self.mixelm = mixelm
        self.percel = percel
        self.mxpart = mxpart
        self.nadvfd = nadvfd
        self.mixelm = mixelm
        self.itrack = itrack
        self.wd = wd
        self.dceps = dceps
        self.nplane = nplane
        self.npl = npl
        self.nph = nph
        self. npmin = npmin
        self.npmax = npmax
        self.interp = 1 # Command-line 'interp' might once be needed if MT3DMS is updated to include other interpolation method
        self.nlsink = nlsink
        self.npsink = npsink
        self.dchmoc = dchmoc
        self.parent.add_package(self)
    def __repr__( self ):
        return 'Advection package class'
    def write_file(self):
        f_adv = open(self.file_name[0], 'w')
        f_adv.write('%10i%10f%10i%10i\n' % (self.mixelm, self.percel, self.mxpart, self.nadvfd))
        if (self.mixelm > 0):
            f_adv.write('%10i%10f\n' % (self.itrack, self.wd))
        if ((self.mixelm == 1) or (self.mixelm == 3)):
            f_adv.write('%10f%10i%10i%10i%10i%10i\n' % (self.dceps, self.nplane, self.npl, self.nph, self. npmin, self.npmax))
        if ((self.mixelm == 2) or (self.mixelm == 3)):
            f_adv.write('%10i%10i%10i\n' % (self.interp, self.nlsink, self.npsink))
        if (self.mixelm == 3):
            f_adv.write('%10f\n' % (self.dchmoc))
        f_adv.close()
