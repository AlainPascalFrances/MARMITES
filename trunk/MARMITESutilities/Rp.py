# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:25:53 2014

@author: apf
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
#import CreateColors

# INPUT
phi = 0.3413
theta_fc = 0.31
theta_fcl = 0.325
Ksat = 305.0
theta_real = [0.3413,0.331550151,0.32628381,0.322867804,0.320768384,0.318704548,0.316569544,0.315039458,0.313367039,0.311089702,0.310733868,0.310555951]
Rsoil_real = [303.87532, 89.62896686, 38.73018136, 12.56113989, 8.126176025, 2.232065299, 1.652781564, 0.937723985, 1.351587854, 0.182045506, 0.121716472, 0.103027722]

incr = 0.0001
lim_low = theta_fc-incr
lim_up  = phi+incr

fig = plt.figure(figsize=(11.7, 8.27))
ax1=fig.add_subplot(2,2,1)
plt.setp(ax1.get_xticklabels(), fontsize=8)
plt.setp(ax1.get_yticklabels(), fontsize=8) 
#plt.title(r'a) - $T_g/PT_g$ function of $\theta $', fontsize = 12)
theta = np.arange(theta_fcl,phi, incr)
theta
Rsoil = Ksat * (theta-theta_fcl)/(phi-theta_fcl)
ax1.plot([theta_fc, theta_fc], [0, Ksat], ':', color = 'black')
ax1.plot([theta_fcl, theta_fcl], [0, Ksat], ':', color = 'black')
ax1.plot([phi, phi], [0, Ksat], ':', color = 'black')
ax1.plot([0, phi], [Ksat, Ksat], ':', color = 'black')
ax1.plot(theta, Rsoil, '--', color = 'darkgray')
ax1.plot(theta, Rsoil, color = 'black')
ax1.plot(theta_real, Rsoil_real, '--o', color = 'gray', )
#    ax1.plot(theta, Esoil, color = 'black')
#ax1.annotate('$c = %.1f$' % c, (x_lbl, y_lbl), horizontalalignment='right', verticalalignment='bottom', fontsize = 10, color = 'black')
plt.grid(True, color = 'lightgray')
plt.xlim((0.30, 0.35))
plt.ylim((0.0,325))
plt.xlabel(r'$\theta$ (%)')#, labelpad=20)
plt.ylabel(r'$R_p$  ($mm.d^{-1}$)')
ax1.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
ax1.text(theta_fc+0.00075, 150, r'$\theta_{fc}$', horizontalalignment='left', verticalalignment='center')#, fontsize = 14)
ax1.text(theta_fcl+0.00075, 150, r"$\theta_{fc}'$", horizontalalignment='left', verticalalignment='center')#, fontsize = 14)
ax1.text(phi+0.00075, 150, r'$\phi$', horizontalalignment='left', verticalalignment='center')#, fontsize = 14)
ax1.text(0.315, Ksat+8, r'$Ksat$', horizontalalignment='left', verticalalignment='center')#, fontsize = 14)

#plt.subplots_adjust(left=0.075, bottom=0.05, right=0.95, top=0.95, wspace=0.25, hspace=0.3)    
#plt.show()
img_path = 'Rp_function'
plt.savefig(img_path,dpi=300)
print('Plot printed:\n%s' % img_path)