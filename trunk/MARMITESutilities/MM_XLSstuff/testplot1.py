﻿import os, tempfile
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
print '\nBackend:\n' + mpl.get_backend() + '\n'

def plot_density(filename,i,t,psi_Na):
    plt.imshow(abs(psi_Na)**2,origin = 'lower')
    filename = os.path.join(tempfile.gettempdir(), filename + '_%04d.png'%i)
    plt.savefig(filename)
    plt.clf()

if __name__ == "__main__":
    x = np.linspace(-6e-6,6e-6,128,endpoint=False)
    y = np.linspace(-6e-6,6e-6,128,endpoint=False)
    X,Y = np.meshgrid(x,y)
    k = 1000000
    omega = 200
    times = np.linspace(0,100e-3,10,endpoint=False)
    fig = plt.figure(figsize=(8,6))
    for i,t in enumerate(times):
        psi_Na = np.sin(k*X-omega*t)
        plot_density('wavefunction',i,t,psi_Na)
        print i
    plt.close()