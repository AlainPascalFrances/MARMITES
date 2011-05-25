import numpy as np
import os as os
import subprocess as sp
import webbrowser as wb

# Global variables
iconst = 1 # Multiplier for individual array elements in integer and real arrays read by MODFLOW's U2DREL, U1DREL and U2DINT.
iprn = -1 # Printout flag. If >= 0 then array values read are printed in listing file.

class basemodel(object):
    'MODFLOW based models base class'
    def __init__(self, modelname = 'modflowtest', namefile_ext = 'nam', exe_name = 'mf2k.exe', model_ws = 'c:/temp'):
        self.__name = modelname
        self.namefile_ext = namefile_ext
        self.namefile = self.__name + '.' + self.namefile_ext
        self.packagelist = []
        self.heading = ''
        self.exe_name = exe_name
        self.model_ws= model_ws
        self.cl_params = ''
    def add_package(self, p):
        for pp in (self.packagelist):
            if (isinstance(p, type(pp))):
                pp = p
                return
        self.packagelist.append( p )
    def array2string(self, a, fmt_str, npl):
        '''Converts a 1D or 2D array into a string
        Input:
            a: array
            fmt_str: format string
            npl: number of numbers per line
        Output:
            s: string representation of the array'''

        aa = np.atleast_2d(a)
        nr, nc = np.shape(aa)[0:2]
        #print 'nr = %d, nc = %d\n' % (nr, nc)

        s = ''
        for r in range(nr):
            for c in range(nc):
                s = s + (fmt_str % a[r, c])
                if (((c + 1) % npl == 0) or (c == (nc - 1))):
                    s = s + '\n'

        return s
    def get_name_file_entries(self):
        s = ''
        for p in self.packagelist:
            for i in range(len(p.name)):
            	s = s + ('%s\t%3i\t%s\n' % (p.name[i], p.unit_number[i], p.file_name[i]))
        return s
    def get_package(self, name):
        for pp in (self.packagelist):
            if (pp.name[0] == name):
                return pp
        return None
    def run_model(self, pause = True, report = None):
        batch_file_name = os.path.join(self.model_ws, 'run.bat')
        error_message = ('Model executable %s not found!' % self.exe_name)
        assert os.path.exists(self.exe_name), error_message

        error_message = ('Name file %s not found!' % self.namefile)
        fn_path = os.path.join(self.model_ws, self.namefile)
        assert os.path.exists(fn_path), error_message

        # Create a batch file to call code so that window remains open in case of error messages
        f = open(batch_file_name, 'w')
        f.write('@ECHO Calling %s with %s\n' % (self.exe_name, self.namefile))
        f.write('%s %s %s\n' % (self.exe_name, self.namefile, self.cl_params))
        if (pause):
           f.write('@PAUSE\n')
        f.close()
        os.path.abspath = self.model_ws
        sp.call(batch_file_name, cwd=self.model_ws, stdout = report)
        os.remove(batch_file_name)
    def write_array(self, f, a, locat, write_fmtin, npi, npl, description=''):
        '''Writes an array of reals to a file.
        Input:
            f: handle to file
            a: array to be written
            locat: FORTRAN unit number
            write_fmtin: flag to indicate if format string must be written
            npi: number of characters per item
            npl: number of numbers per line
            description: string that is appended to the layer header

            Notes:
            npl < 0 indicates that this number of items per line must be enforced
            and no 'Constant' line will be written'''

        aa = np.atleast_3d(a)
        nr, nc, nl = np.shape(aa)[0:3]

        if (a.dtype == 'int32'):
            fmtin = '(%dI%d) ' % (abs(npl), npi) # FORTRAN format descriptor
            fmt_str = '%' + ('%d' % npi) + 'd'
        else:
            fmtin = '(%dG%d.0) ' % (abs(npl), npi) # FORTRAN format descriptor
            fmt_str = '%' + ('%d' % npi) + 'g' #' %13f'

        for l in range(nl):
            if ((npl > 0) and (aa[:, :, l].min() == aa[:, :, l].max())):
                if (nl > 1):
                    if a.dtype == 'int32':
                        f.write('%10s%10d %s %d\n' % ('Constant ', aa[0, 0, l], description, l + 1))
                    else:
                        f.write('%10s%10g %s %d\n' % ('Constant ', aa[0, 0, l], description, l + 1))
                else:
                    if a.dtype == 'int32':
                        f.write('%10s%10d %s\n' % ('Constant ', aa[0, 0, l], description))
                    else:
                        f.write('%10s%10g %s\n' % ('Constant ', aa[0, 0, l], description))
            else:
                if (write_fmtin):
                    if (nl > 1):
                        f.write('%10d%10d%20s%10d %s %d\n' % (locat, iconst, fmtin, iprn, description, l + 1))
                    else:
                        f.write('%10d%10d%20s%10d %s\n' % (locat, iconst, fmtin, iprn, description))

                f.write(self.array2string(aa[:, :, l], fmt_str, abs(npl)))
    def write_input(self):
        print self # Same as calling self.__repr__()
        print 'Writing packages:'
        for p in self.packagelist:
            p.write_file()
            print p.__repr__()
        self.write_name_file()
    def write_name_file(self):
        '''Every package needs its own writenamefile function'''
        print 'IMPLEMENTATION ERROR: writenamefile must be overloaded'
    def get_name(self):
        return self.__name
    def set_name(self, value):
        self.__name = value
        self.namefile = self.__name + '.' + self.namefile_ext
        for p in self.packagelist:
            for i in range(len(p.extension)):
                p.file_name[i] = self.__name + '.' + p.extension[i]
    name = property(get_name, set_name)

class package(object):
    'General package class'
    def __init__(self, parent, extension='glo', name='GLOBAL', unit_number=1):
        self.parent = parent # To be able to access the parent modflow object's attributes
        if (not isinstance(extension, list)):
            extension = [extension]
        self.extension = []
        self.file_name = []
        for e in extension:
            self.extension = self.extension + [e]
            self.file_name = self.file_name + [self.parent.name + '.' + e]
            self.fn_path = os.path.join(self.parent.model_ws,self.file_name[0])
        if (not isinstance(name, list)):
            name = [name]
        self.name = name
        if (not isinstance(unit_number, list)):
            unit_number = [unit_number]
        self.unit_number = unit_number
        self.url = 'http://water.usgs.gov/nrp/gwsoftware/modflow2005/Guide/index.html'
    def __repr__( self ):
        s = self.__doc__
        exclude_attributes = ['extension', 'heading', 'name', 'parent', 'url']
        for attr, value in sorted(self.__dict__.iteritems()):
            if not (attr in exclude_attributes):
                if (isinstance(value, list)):
                    if (len(value) == 1):
                        s = s + ' %s = %s (list)\n' % (attr, str(value[0]))
                    else:
                        s = s + ' %s (list, items = %d)\n' % (attr, len(value))
                elif (isinstance(value, np.ndarray)):
                    s = s + ' %s (array, shape = %s)\n' % (attr, value.shape.__str__()[1:-1] )
                else:
                    s = s + ' %s = %s (%s)\n' % (attr, str(value), str(type(value))[7:-2])
        return s
    def assignarray(self, destarray, srcarray):
        if (isinstance(srcarray, list)): srcarray = np.array(srcarray, dtype = destarray.dtype)
        if (np.isscalar(srcarray)):
            destarray[:] = srcarray
        else:
            assert destarray.dtype == srcarray.dtype, 'Error: Array types do not match'
            if srcarray.shape == destarray.shape:
                destarray[:] = srcarray
            elif ((srcarray.ndim == 1) and (destarray.ndim == 3) and (srcarray.size == destarray.shape[2])):
                # If srcarray is 1D and has a length equal to the number of layers
                # then assign the elements in srcarray as constant values to all
                # cells in the corresponding layers in destarray
                for l in range(srcarray.size):
                    destarray[ :, :, l ] = srcarray[l]
            elif ((srcarray.ndim == 2) and (destarray.ndim == 3) and (srcarray.shape == destarray[ :, :, 0].shape)):
                # If srcarray is 2D and the number of rows and columns matches
                # the number of rows and columns of destarray then assign the
                # values in srcarray to all the layers in destarray
                for l in range(destarray.shape[-1]):
                    destarray[ :, :, l ] = srcarray.copy()
            else:
                raise BaseException, 'Error: Cannot assign values to grid. Dimensions do not match.'
    def assign_layer_row_column_data(self, layer_row_column_data, ncols):
        if (layer_row_column_data != None):
            new_layer_row_column_data = []
            mxact = 0
            for a in layer_row_column_data:
                a = np.atleast_2d(a)
                nr, nc = a.shape
                assert nc == ncols, 'layer_row_column_data must have %d columns' % ncols
                mxact = max(mxact, nr)
                new_layer_row_column_data.append(a)
            return mxact, new_layer_row_column_data
        return
    def webdoc(self):
        if self.parent.version == 'mf2k':
            wb.open('http://water.usgs.gov/nrp/gwsoftware/modflow2000/Guide/' + self.url)
        else:
            wb.open('http://water.usgs.gov/nrp/gwsoftware/modflow2005/Guide/' + self.url)
    def write_file(self):
        '''Every package needs its own write_file function'''
        print 'IMPLEMENTATION ERROR: write_file must be overloaded'
    def write_layer_row_column_data(self, f, layer_row_column_data):
        for n in range(self.parent.get_package('DIS').nper):
            if (n < len(layer_row_column_data)):
                a = layer_row_column_data[n]
                itmp = a.shape[0]
                f.write('%10i%10i\n' % (itmp, self.np))
                for b in a:
                    f.write('%10i%10i%10i' % (b[0], b[1], b[2]) )
                    for c in b[3:]:
                        f.write('%13.6g' % c)
                    f.write('\n')
            else:
                itmp = -1
                f.write('%10i%10i\n' % (itmp, self.np))
