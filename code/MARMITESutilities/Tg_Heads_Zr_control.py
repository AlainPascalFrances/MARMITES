# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 11:01:50 2014

@author: apf
"""
import numpy as np

class cMF():
    def __init__(self):
        self.wel_yn = 1
        self.hnoflo = 9999.9
        self.outcropL = np.array([[1]])
        self.sy_actual = [np.array([[0.01]])]

nsl = 2
NVEG = 3
kT_min = [0.01, 0.5, 0.1]
kT_max = [0.02, 0.8, 0.9]
kT_n = [0.5, 2.0 , 2.0]
elevation = 100.0
HEADS_corr = 1000.0*(elevation - 2.0)
Zr_elev = 1000.0*np.array([0.2, 7.5, 15.0])
Zr_elev = HEADS_corr - Zr_elev
Ssoil_pc_tmp = [0.20,0.15]
Sr = [0.05, 0.075]
Sm = [0.3, 0.2]
PT =[1.0,1.5,1.75]
VEGarea = [25.0,25.0,25.0]
dgwt_corr = 10000.0-HEADS_corr
i = 0
j = 0
L = 0

cMF = cMF()

Tg_tmp_Zr = np.zeros([nsl, len(Zr_elev)])
Tg_tmp = 0.0
if cMF.wel_yn == 1:
    order = np.asarray(Zr_elev).argsort()
    for jj, (Zr_elev_, v, kT_min_, kT_max_, kT_n_) in enumerate(zip(np.asarray(Zr_elev)[order],np.asarray(range(NVEG))[order], np.asarray(kT_min)[order], np.asarray(kT_max)[order], np.asarray(kT_n)[order])):
        print '----'
        if HEADS_corr > Zr_elev_:
            for l in range(nsl):
                print 'veg%d, soil layer %d' %(v,l)
                if Ssoil_pc_tmp[l] != cMF.hnoflo:
                    if Ssoil_pc_tmp[l] >= Sr[l]:
                        if Ssoil_pc_tmp[l] <= Sm[l]:
                            kT = kT_min_+(kT_max_-kT_min_)*np.power(1-np.power(np.abs(Ssoil_pc_tmp[l]-Sm[l])/(Sm[l]-Sr[l]),kT_n_),1/kT_n_)
                        else:
                            kT = kT_max_
                    else:
                        kT = kT_min_
                    Tg_tmp_Zr[l,v] = PT[v]*kT
                    PT[v] -= Tg_tmp_Zr[l,v]
                    Tg_tmp1 = (Tg_tmp_Zr[l,v]*VEGarea[v]*0.01)
                    if Tg_tmp1 > 1E-6 and (HEADS_corr-Tg_tmp1/cMF.sy_actual[cMF.outcropL[i,j]-1][i][j]) < Zr_elev_:
                        print 'WARNING at cell [i=%d,j=%d,L=%d] with NVEG = %d\nTg = %.2f, HEADS_corr = %.2f, dgwt_corr= %.2f, Zr_elev = %.2f, Sy = %.3f' % (i+1,j+1,cMF.outcropL[i,j], v, Tg_tmp1, HEADS_corr, dgwt_corr, Zr_elev_, cMF.sy_actual[cMF.outcropL[i,j]-1][i][j])
                        Tg_tmp1 = (HEADS_corr - Zr_elev[v])*cMF.sy_actual[cMF.outcropL[i,j]-1][i][j]
                        print 'New Tg = %.2f mm' % Tg_tmp1
                    else:                        
                        print 'Tg veg%d = %.2f mm' %(v,Tg_tmp1)
                    Tg_tmp += Tg_tmp1
                    dgwt_corr += Tg_tmp1/cMF.sy_actual[cMF.outcropL[i,j]-1][i][j]
                    HEADS_corr -= Tg_tmp1/cMF.sy_actual[cMF.outcropL[i,j]-1][i][j]                                                   
        else:
            print 'veg%d: root too short!' % v
else:
    Tg_tmp = 0.0
print '----\nTg tot = %.2f mm' % Tg_tmp
print 'Done!'