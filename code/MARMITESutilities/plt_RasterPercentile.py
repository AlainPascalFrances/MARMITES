# -*- coding: utf-8 -*-
"""
Created on Mon Nov 03 21:20:15 2014

@author: frances.alain@gmail.com
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#import os

ncol = 6
nrow = 13
fmt = '%2.f'
ntick = 1
cmap = plt.cm.gist_rainbow_r

# Store some arrays for plotting
x = np.arange(0.5, ncol+1.5, 1)
y = np.arange(0.5, nrow+1.5, 1)
xg,yg = np.meshgrid(x,y)

x = np.arange(1, ncol+1, 1)
y = np.arange(1, nrow+1, 1)
xg1,yg1 = np.meshgrid(x,y)

Vtmp = np.asarray([[ 105.      ,  104.      ,  101.      ,  103.5     ,  106.      ,
         107.      ],
       [ 102.22843 ,  101.24186 ,   98.5     ,  101.25999 ,  102.96385 ,
         105.      ],
       [ 101.80586 ,  100.76963 ,   98.      ,   99.679477,  101.12545 ,
         102.5     ],
       [ 104.      ,  103.5     ,   97.5     ,   97.      ,   -999.9,
         -999.9 ],
       [ 105.      ,  103.      ,  104.      ,   96.5     ,   98.67086 ,
          -999.9],
       [ 100.59955 ,   99.368366,   95.5     ,   96.      ,   98.407378,
         101.      ],
       [  97.771869,   97.20416 ,   95.      ,   96.867147,   97.884638,
          98.939042],
       [  95.5     ,   96.100089,   94.5     ,   96.150685,   97.115566,
          97.115566],
       [  96.2914  ,   95.555256,   94.      ,   95.485017,   95.485017,
          95.485017],
       [  97.      ,   95.209708,   93.5     ,   94.924793,   96.336175,
          96.336175],
       [  95.605343,   94.488426,   93.      ,   94.362053,   95.490359,
          96.438411],
       [  94.338686,   93.821578,   92.5     ,   93.825062,   94.687194,
          95.262802],
       [  94.      ,   93.488017,   92.      ,   93.      ,   94.5     ,
          95.      ]], dtype = np.float)

Vtmp = np.ma.masked_values(Vtmp, -999.9, atol = 0.09)

ticks = np.percentile(Vtmp.compressed().flatten(),[0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0])
norm = mpl.colors.BoundaryNorm(ticks, cmap.N)

          
Vmin_tmp = np.min(Vtmp)
Vmax_tmp = np.max(Vtmp)
interval_diff = 1.0
files_tmp = []
fig = plt.figure()
figtitle = fig.suptitle('')
figtitle.set_text("Raster percentile color map")
plt.draw()
ax = fig.add_subplot(1,1,1, axisbg = 'silver')
ax.xaxis.set_ticks(np.arange(0,ncol+1,ntick))
ax.yaxis.set_ticks(np.arange(0,nrow+1,ntick))
plt.setp(ax.get_xticklabels(), fontsize=8)
plt.setp(ax.get_yticklabels(), fontsize=8)
plt.ylabel('row i', fontsize=10)
plt.xlabel('col j', fontsize=10)
ims=ax.pcolormesh(xg, yg, Vtmp, cmap = cmap, vmin = Vmin_tmp, vmax = Vmax_tmp, norm = norm)
ax.set_ylim(bottom = np.max(yg1), top = np.min(yg1))
ax.axis('scaled')
if max(x) > max(y):
    cax = fig.add_axes([.125, 0.035, 0.75, 0.025])
    CBorient = 'horizontal'
else:
    cax = fig.add_axes([0.035, 0.125, 0.025, 0.75])
    CBorient = 'vertical'
CB = fig.colorbar(ims, extend='both', ticks = ticks, cax = cax,  orientation = CBorient)
CB.set_label('elevation (m)', fontsize = 12)
if max(x) > max(y):
    cax.xaxis.set_label_position('top')
    plt.setp(CB.ax.get_xticklabels(), fontsize = 7)
else:
    cax.yaxis.set_label_position('left')
    plt.setp(CB.ax.get_yticklabels(), fontsize = 7)
#plt_export_fn = os.path.join(r'e:\temp', '__test.png')
#plt.savefig(plt_export_fn)
plt.show()