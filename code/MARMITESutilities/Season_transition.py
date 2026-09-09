#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      alf
#
# Created:     27-03-2011
# Copyright:   (c) alf 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import matplotlib.pyplot as plt
import os

J = []
for j in range(365):
    J.append(j+1)
J_vd = [150]
J_vw = [270]
TRANS_vdw = [20]
alfa_vd = [0.3]
alfa_vw = [0.5]

v = 0
Rns_v = []
alfa = []
for j in range(len(J)):
    if J[j] < J_vd[v] or J[j] > J_vw[v]:
        alfa.append(alfa_vw[v])
    else:
        if J[j] < J_vd[v] + TRANS_vdw[v]:
            alfa.append(alfa[j-1]-(alfa_vw[v]-alfa_vd[v])/(TRANS_vdw[v]+1))
        elif J[j] > J_vw[v] - TRANS_vdw[v]:
            alfa.append(alfa[j-1]+(alfa_vw[v]-alfa_vd[v])/(TRANS_vdw[v]+1))
        else:
            alfa.append(alfa_vd[v])
fig = plt.figure()
ax1=fig.add_subplot(111)
plt.plot(J,alfa)
plt.show()
