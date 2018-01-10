#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      frances08512
#
# Created:     30-03-2011
# Copyright:   (c) frances08512 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

'''
    Brooks-Corey function
    UZF1 package for MODFLOW-2005
    TM6-A19 pag6
    http://pubs.usgs.gov/tm/2006/tm6a19
'''

lblspc = 0.05
mkscale = 0.5
fig = plt.figure(figsize=(11.7, 8.27))

# example of Table 1 UZF manual TM6-A19
Ks = 4.0E-6 # saturated vertical hydraulic conductivity (m/s)
Ts = 0.40 # saturated water content
Tr = 0.2 # residual water content, assimilated to specific retention
eps = [1.0, 2.0, 3.5, 5.0]
Sy = Ts - Tr # specific yield
dT = 0.01
T = np.arange(Tr,Ts+dT,dT)

fig.suptitle('Brooks-Corey function\n' + r'$\theta_{s}=%0.2f,\ \theta_{r}=%0.2f,\ K_{sat}=%.1E\ m.s^{-1}$' % (Ts,Tr,Ks))

ax1=fig.add_subplot(1,1,1)
plt.setp(ax1.get_xticklabels(), fontsize=8)
plt.setp(ax1.get_yticklabels(), fontsize=8)
for i, e in enumerate(eps):
    K_T = Ks*np.power((T-Tr)/(Ts-Tr),e)
    ax1.plot(T, K_T, '-', label = str(e))
plt.grid(True)
leg = ax1.legend(loc=0)
plt.ylabel(r'$K(\theta)\ (m.s^{-1})$')
plt.xlabel(r'$\theta\ (\%)$')
ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1E'))
ax1.set_ylim(ax1.get_ylim()[::-1])
plt.xlim((Tr,Ts))

##z0 = 5.0   # NOT CORRECT!
##T0 = Sy    # NOT CORRECT!
##for t in T:
##    v_T.append((eps*Ks/Sy)*np.math.pow((t-Tr)/Sy,eps-1))
##    z_T.append(z0*np.math.pow((t-Tr)/(T0-Tr),eps-1))
##ax2=fig.add_subplot(1,3,2)
##plt.setp(ax2.get_xticklabels(), fontsize=8)
##plt.setp(ax2.get_yticklabels(), fontsize=8)
##ax2.plot(T, v_T, '-', color='brown')
##plt.grid(True)
###plt.xlim((-0.05, 1.05))
###plt.ylim((-100,3000))
##plt.ylabel(r'$v(\theta)\ (m.s^{-1})$')
##plt.xlabel(r'$\theta\ (\%)$')
##ax2.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1E'))
##ax2.set_ylim(ax2.get_ylim()[::-1])
##
##ax3=fig.add_subplot(1,3,3)
##plt.setp(ax3.get_xticklabels(), fontsize=8)
##plt.setp(ax3.get_yticklabels(), fontsize=8)
##ax3.plot(T, z_T, '-', color='brown')
##plt.grid(True)
###plt.xlim((-0.05, 1.05))
###plt.ylim((-100,3000))
##plt.ylabel(r'$z(\theta)\ (m)$')
##plt.xlabel(r'$\theta\ (\%)$')
##ax3.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
##ax3.set_ylim(ax3.get_ylim()[::-1])

#plt.show()
img_path = 'BrooksCorey_function1.png'
plt.savefig(img_path,dpi=300)
print 'Plot printed:\n%s' % img_path