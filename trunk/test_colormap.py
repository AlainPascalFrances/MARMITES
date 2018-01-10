# -*- coding: utf-8 -*-
"""
Created on Thu Oct 02 09:08:00 2014

@author: apf
"""

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors as clr

cdict = {'red':  ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.8, 1.0),
                   (0.75,1.0, 1.0),
                   (1.0, 0.4, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.9, 0.9),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.4),
                   (0.25,1.0, 1.0),
                   (0.5, 1.0, 0.8),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }

#cdict['alpha'] = ((0.0, 1.0, 1.0),
#                #   (0.25,1.0, 1.0),
#                   (0.5, 0.3, 0.3),
#                #   (0.75,1.0, 1.0),
#                   (1.0, 1.0, 1.0))

plt.register_cmap(name='BlueRed3', data=cdict)
#plt.register_cmap(name='BlueRedAlpha', data=cdict)

x = np.arange(0, np.pi, 0.1)
y = np.arange(0, 2*np.pi, 0.1)
X, Y = np.meshgrid(x,y)
Z = np.cos(X) * np.sin(Y) * 10

# Make the figure:

plt.figure(figsize=(6,9))
plt.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)

# Make 4 subplots:

plt.subplot(1,1,1)
plt.imshow(Z, interpolation='nearest')
plt.colorbar()
plt.set_cmap('BlueRed3')

plt.show()
