from numpy import empty
from mbase import package

class mtssm(package):
    'Sink & Source Mixing package class'
    def __init__(self, model, cchd=1.0, crch=0.0, cpbc=1.0, extension='ssm'):
        package.__init__(self, model, extension, 'SSM', 34) # Call ancestor's init to set self.parent, extension, name and unit number
        nrow, ncol, nlay, nper = self.parent.mf.nrow_ncol_nlay_nper
        mfchd = self.parent.mf.get_package('CHD')
        mfrch = self.parent.mf.get_package('RCH')
        mfpbc = self.parent.mf.get_package('PBC')
        # The assignments below do not yet check if the dimensions of crch, cpbc, ..., etc. are
        # compatible with their corresponding MODFLOW structures! Better to fix...
        if (mfchd != None):
            self.cchd = []
            if (not isinstance(cchd, list)):
                cchd = [cchd]
            for a in cchd:
                b = empty((mfchd.layer_row_column_shead_ehead[0].shape[0]))
                self.assignarray(b , a )
                self.cchd = self.cchd + [b]
        if (mfrch != None):
            self.crch = []
            if (not isinstance(crch, list)):
                crch = [crch]
            for a in crch:
                crch_t = []
                if (not isinstance(a, list)):
                    a = [a]
                for b in a:
                    c = empty((nrow, ncol))
                    self.assignarray(c , b )
                    crch_t.append(c)
                self.crch.append(crch_t)
        if (mfpbc != None):
            self.cpbc = []
            if (not isinstance(cpbc, list)):
                cpbc = [cpbc]
            for a in cpbc:
                b = empty((mfpbc.layer_row_column_shead_ehead[0].shape[0]))
                self.assignarray(b , a )
                self.cpbc = self.cpbc + [b]
        self.parent.add_package(self)
    def __repr__( self ):
        return 'Sink & Source Mixing package class'
    def write_file(self):
        nrow, ncol, nlay, nper = self.parent.mf.nrow_ncol_nlay_nper
        mfchd = self.parent.mf.get_package('CHD')
        mfdrn = self.parent.mf.get_package('DRN')
        mfevt = self.parent.mf.get_package('EVT')
        mfghb = self.parent.mf.get_package('GHB')
        mfrch = self.parent.mf.get_package('RCH')
        mfpbc = self.parent.mf.get_package('PBC')
        mfriv = self.parent.mf.get_package('RIV')
        mfwel = self.parent.mf.get_package('WEL')
        # Open file for writing
        f_ssm = open(self.file_name[0], 'w')
        maxssm = 0
        if (mfwel == None):
            f_ssm.write('%2s' % ('F'))
        else:
            f_ssm.write('%2s' % ('T'))
            maxssm = maxssm + mfwel.ncells()
        if (mfdrn == None):
            f_ssm.write('%2s' % ('F'))
        else:
            f_ssm.write('%2s' % ('T'))
            maxssm = maxssm + mfdrn.ncells()
        if (mfrch == None):
            f_ssm.write('%2s' % ('F'))
        else:
            f_ssm.write('%2s' % ('T'))
            maxssm = maxssm + mfrch.ncells()
        if (mfevt == None):
            f_ssm.write('%2s' % ('F'))
        else:
            f_ssm.write('%2s' % ('T'))
            maxssm = maxssm + mfevt.ncells()
        if (mfriv == None):
            f_ssm.write('%2s' % ('F'))
        else:
            f_ssm.write('%2s' % ('T'))
            maxssm = maxssm + mfriv.ncells()
        if (mfghb == None):
            f_ssm.write('%2s' % ('F'))
        else:
            f_ssm.write('%2s' % ('T'))
            maxssm = maxssm + mfghb.ncells()
        # CHD package
        if (mfchd != None):
            maxssm = maxssm + mfchd.ncells()
        # PBC package
        if (mfpbc != None):
            maxssm = maxssm + mfpbc.ncells()
        f_ssm.write('%2s%2s%2s%2s\n' % ('F', 'F', 'F', 'F'))
        f_ssm.write('%10d\n' % (maxssm))
        for n in range(nper):
            # Recharge
            if (mfrch != None):
                if (n < len(self.crch)):
                    incrch = 1
                else:
                    incrch = -1
                f_ssm.write('%10i\n' % (incrch))
                if (n < len(self.crch)):
                    for s in range(len(self.crch[n])):
                        comment = 'Recharge concentration array of species %d for stress period %d' % (s + 1, n + 1)
                        self.parent.write_array( f_ssm, self.crch[n][s], self.unit_number[0], True, 13, -ncol, comment )
            # Evapotranspiration
            if (mfevt != None):
                print 'Not implemented yet!'
                break
            # WEL, DRN, EVT, RIV, GHB, PBC
            nss = 0
            need_nss = False
            if (mfchd != None):
                if (n < len(mfchd.layer_row_column_shead_ehead)):
                    ssmchd = mfchd.layer_row_column_shead_ehead[n]
                    need_nss = True
                else:
                    ssmchd = mfchd.layer_row_column_shead_ehead[-1]
                nss = nss + ssmchd.shape[0]
            if (mfpbc != None):
                if (n < len(mfpbc.layer_row_column_shead_ehead)):
                    ssmpbc = mfpbc.layer_row_column_shead_ehead[n]
                    need_nss = True
                else:
                    ssmpbc = mfpbc.layer_row_column_shead_ehead[-1]
                nss = nss + ssmpbc.shape[0]
            if (need_nss == True):
                f_ssm.write('%10i\n' % (nss))
                if (mfchd != None):
                    if (n < len(self.cchd)):
                        c = self.cchd[n]
                    else:
                        c = self.cchd[-1]
                    i = 0
                    for b in ssmchd:
                        f_ssm.write('%10i%10i%10i%10f%10i\n' % (b[0], b[1], b[2], c[i], 1) )
                        i = i + 1
                if (mfpbc != None):
                    if (n < len(self.cpbc)):
                        c = self.cpbc[n]
                    else:
                        c = self.cpbc[-1]
                    i = 0
                    for b in ssmpbc:
                        f_ssm.write('%10i%10i%10i%10f%10i\n' % (b[0], b[1], b[2], c[i], 51) )
                        i = i + 1
            else:
               f_ssm.write('%10i\n' % (nss))
        f_ssm.close()


