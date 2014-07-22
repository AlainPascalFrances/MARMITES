
import h5py, tempfile, os, sys, traceback
import numpy as np
import pylab as pylb

def plot(fn):
    #read the data on HDF5 file and plot them
    h5 = h5py.File(h5_fn)
    x = h5['x']
    th = h5['theory']
    sig = h5['signal']
    pylb.plot(x,th,'-',x,sig,'+')
    pylb.show()
    h5.close()

#create data
npts=500
x=np.linspace(0,10,npts)
theory=5*np.sin(2*np.pi*x/2.0)
noise=np.random.normal(0,np.std(theory)/5,npts)
sig=theory+noise

try:
    #store data into HDF5 file
    h5_fn = os.path.join(tempfile.gettempdir(),'_h5_MM.h5')
    h5 = h5py.File(h5_fn, 'w')
    x = 1+a
    h5.create_dataset(name = 'x', data = np.asarray(x))
    h5.create_dataset(name = 'theory', data = np.asarray(theory))
    h5.create_dataset(name = 'signal', data = np.asarray(sig))
    h5.create_dataset(name = 'noise', data = np.asarray(noise))
    h5.close()
    del x, theory, sig, noise
except StandardError, e:
    print '\nERROR!\nError description:'
    traceback.print_exc(file=sys.stdout)

plot(h5_fn)
print 'Done'
