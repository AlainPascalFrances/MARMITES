#! /usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
"""
Make signal with gaussian noise
"""
npts=500
x=np.linspace(0,10,npts)
theory=5*np.sin(2*np.pi*x/2.0)
noise=np.random.normal(0,np.std(theory)/5,npts)  #mean, std dev, num pts
sig=theory+noise

plt.plot(x,theory,'-',x,sig,'+')
plt.show()