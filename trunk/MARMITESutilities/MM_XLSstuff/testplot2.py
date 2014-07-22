# BOF
import os, tempfile
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
print '\nBackend:\n' + mpl.get_backend() + '\n'

def plot_density(filename,i,t,psi_Na, image):
    if image is None:
        image = plt.imshow(abs(psi_Na)**2,origin = 'lower')
    else:
        image.set_data(abs(psi_Na)**2)
        fig.canvas.draw()
    filename = os.path.join(tempfile.gettempdir(), filename + '_%04d.png'%i)
    plt.savefig(filename)
    return image

if __name__ == "__main__":
    x = np.linspace(-6e-6,6e-6,128,endpoint=False)
    y = np.linspace(-6e-6,6e-6,128,endpoint=False)
    X,Y = np.meshgrid(x,y)
    k = 1000000
    omega = 200
    times = np.linspace(0,100e-3,100,endpoint=False)
    fig = plt.figure(figsize=(8,6))
    image = None
    for i,t in enumerate(times):
        psi_Na = np.sin(k*X-omega*t)
        image = plot_density('wavefunction',i,t,psi_Na, image)
        print i
    del image, fig
    plt.close()
# EOF