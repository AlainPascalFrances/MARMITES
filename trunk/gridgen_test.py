import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
flopypth = os.path.dirname(sys.executable) + r'\Lib\site-packages\flopy'
if flopypth not in sys.path:
    sys.path.append(flopypth)
import flopy
from flopy.utils.gridgen import Gridgen

print(sys.version)
print('numpy version: {}'.format(np.__version__))
print('matplotlib version: {}'.format(mpl.__version__))
print('flopy version: {}'.format(flopy.__version__))

Lx = 100.
Ly = 100.
nlay = 2
nrow = 51
ncol = 51
delr = Lx / ncol
delc = Ly / nrow
h0 = 10
h1 = 5
top = h0
botm = np.zeros((nlay, nrow, ncol), dtype=np.float32)
botm[1, :, :] = -10.

model_ws = r"E:\00code_ws\00_TESTS\gridgen" #os.path.join('.', 'data')
ms = flopy.modflow.Modflow(modelname='gridgen_test', model_ws=model_ws, exe_name=r'C:\00MODFLOW\mf6.0.3\bin\mf6.exe', rotation=-20.)
dis = flopy.modflow.ModflowDis(ms, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr,
                               delc=delc, top=top, botm=botm)

g = Gridgen(dis, model_ws=model_ws, exe_name=r'C:\00MODFLOW\gridgen.1.0.02\bin\gridgen_x64.exe')

# setup the active domain
adshp = os.path.join(model_ws, 'ad0')
adpoly = [[[(0, 0), (0, 60), (40, 80), (60, 0), (0, 0)]]]
# g.add_active_domain(adpoly, range(nlay))

# RIVER
river = [[[(-20, 10), (60, 60)]]]
g.add_refinement_features(river, 'line', 3, range(nlay))
rf1shp = os.path.join(model_ws, 'rf1')

# REFINE GRID
g.add_refinement_features(adpoly, 'polygon', 1, range(nlay))
rf2shp = os.path.join(model_ws, 'rf2')

#WELL
x = Lx * np.random.random(10)
y = Ly * np.random.random(10)
wells = list(zip(x, y))
g.add_refinement_features(wells, 'point', 3, range(nlay))
rf0shp = os.path.join(model_ws, 'rf0')

# Plot the Gridgen Input
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
mm = flopy.plot.ModelMap(model=ms)
mm.plot_grid()
flopy.plot.plot_shapefile(rf2shp, ax=ax, facecolor='lightgrey', edgecolor='none')
flopy.plot.plot_shapefile(rf1shp, ax=ax, linewidth=10, color = 'darkblue')
flopy.plot.plot_shapefile(rf0shp, ax=ax, facecolor='red', radius=1)
mpl.pyplot.savefig(os.path.join(model_ws, 'fig1.png'))

g.build(verbose=True)

# Plot the Grid
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
g.plot(ax, linewidth=0.5)
flopy.plot.plot_shapefile(rf2shp, ax=ax, facecolor='lightgrey', edgecolor='none', alpha=0.2)
flopy.plot.plot_shapefile(rf1shp, ax=ax, linewidth=10, alpha=0.2, color = 'darkblue')
flopy.plot.plot_shapefile(rf0shp, ax=ax, facecolor='red', radius=1, alpha=0.2)
mpl.pyplot.savefig(os.path.join(model_ws, 'fig2.png'))

mu = flopy.modflow.Modflow(model_ws=model_ws, modelname='mfusg')
disu = g.get_disu(mu)
disu.write_file()
# print(disu)

adpoly_intersect = g.intersect(adpoly, 'polygon', 0)
print(adpoly_intersect.dtype.names)
print(adpoly_intersect)
print(adpoly_intersect.nodenumber)

well_intersect = g.intersect(wells, 'point', 0)
print(well_intersect.dtype.names)
print(well_intersect)
print(well_intersect.nodenumber)

river_intersect = g.intersect(river, 'line', 0)
print(river_intersect.dtype.names)
# print(river_intersect)
# print(river_intersect.nodenumber)

# Plot Intersected Features
a = np.zeros((g.nodes), dtype=np.int)
a[adpoly_intersect.nodenumber] = 1
a[well_intersect.nodenumber] = 2
a[river_intersect.nodenumber] = 3
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
g.plot(ax, a=a, masked_values=[0], edgecolor='none', cmap='jet')
flopy.plot.plot_shapefile(rf2shp, ax=ax, facecolor='lightgrey', alpha=0.25)
mpl.pyplot.savefig(os.path.join(model_ws, 'fig3.png'))