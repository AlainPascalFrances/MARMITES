import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

f = plt.figure()
f.suptitle('Calibration criteria', fontsize = 16)
lbl = ['','1/4','1/2','3/4','1','4/3','2','4']
fontsize = 8
ax = []
CS = []
CB = []

rmse = np.arange(49).reshape((7,7)).astype(float)
rmse = [[3.05,3.81,4.74,np.nan,np.nan,np.nan,12.99],
[2.87,3.18,3.92,4.62,5.49,7.05,10.49],
[np.nan,2.98,3.65,4.28,5.02,6.54,9.47],
[np.nan,2.86,3.38,4.07,4.72,6.10,9.14],
[np.nan,np.nan,3.17,3.77,4.42,5.67,8.72],
[np.nan,np.nan,2.94,3.27,3.91,4.84,7.93],
[3.16,np.nan,2.81,2.86,np.nan,np.nan,6.09]]

rsr = np.arange(49).reshape((7,7)).astype(float)
rsr = [[4.79,5.52,7.06,np.nan,np.nan,np.nan,22.62],
[4.76,4.60,5.54,6.69,8.34,11.19,17.54],
[np.nan,4.35,5.09,6.05,7.39,10.26,15.60],
[np.nan,4.28,4.68,5.70,6.76,9.44,14.95],
[np.nan,np.nan,4.39,5.25,6.24,8.64,14.14],
[np.nan,np.nan,4.25,4.50,5.45,7.02,12.67],
[5.47,np.nan,4.39,4.18,np.nan,np.nan,9.38]]

nse = np.arange(49).reshape((7,7)).astype(float)
nse = [[-29.90,-42.45,-81.93,np.nan,np.nan,np.nan,-857.70],
[-32.63,-27.77,-44.85,-74.94,-123.46,-227.68,-567.25],
[np.nan,-25.50,-36.42,-59.06,-97.60,-191.23,-450.87],
[np.nan,-25.72,-30.69,-49.58,-81.62,-162.72,-415.46],
[np.nan,np.nan,-27.51,-39.62,-65.92,-136.92,-371.29],
[np.nan,np.nan,-25.27,-28.94,-43.91,-88.85,-293.52],
[-40.75,np.nan,-28.08,-24.18,np.nan,np.nan,-157.45]]

r = [[0.69,0.71,0.71,np.nan,np.nan,np.nan,0.69],
[0.75,0.77,0.77,0.77,0.77,0.76,0.77],
[np.nan,0.76,0.76,0.76,0.76,0.76,0.76],
[np.nan,0.74,0.72,0.72,0.72,0.72,0.74],
[np.nan,np.nan,0.62,0.59,0.59,0.60,0.67],
[np.nan,np.nan,0.45,0.42,0.39,0.39,0.40],
[0.36,np.nan,0.17,0.10,np.nan,np.nan,0.18]]

for p, (sp, title, array) in enumerate(zip([221,222,223,224],['RMSE','RSR','NSE','r'],[rmse, rsr, nse, r])):
    ax.append(f.add_subplot(sp))
    ax[p].set_title(title)
    masked_array = np.ma.array(array, mask=np.isnan(array))
    if title == 'RMSE' or title == 'RSR':
        cmap = mpl.cm.jet
    else:
        cmap = mpl.cm.jet_r
    #cmap.set_bad('w',1.)
    CS.append(ax[p].imshow(masked_array, interpolation='nearest', cmap=cmap))
    ax[p].set_ylim(ax[p].get_ylim()[::-1])
    ax[p].set_xticklabels(lbl, fontsize = fontsize)
    ax[p].set_yticklabels(lbl, fontsize = fontsize)
    CB.append(plt.colorbar(CS[p]))
    plt.setp(CB[p].ax.get_yticklabels(), fontsize = fontsize)
    if sp == 221 or sp ==223:
        ax[p].set_ylabel('$S_{y,MRS}$ mult.')
    if sp == 223 or sp == 224:
        ax[p].set_xlabel('$K_{MRS}$ mult.')
    f.canvas.draw()

img_fn = os.path.join(os.path.expanduser('~'), 'Desktop\calib_criteria.png')
plt.savefig(img_fn)

print '\nDone!\nImage save in:\n%s' % img_fn