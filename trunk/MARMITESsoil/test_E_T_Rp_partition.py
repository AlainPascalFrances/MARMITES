import matplotlib.pyplot as plt

def evp(s_tmp,Sm,Sr, pet, dt):
    '''
    Actual evapotranspiration function
    '''
    global evp_tmp
    if (s_tmp - Sr) < 1.0E-7:
        evp_tmp = 0.0
    elif ((s_tmp - Sr) > -1.0E-7) and (s_tmp - Sm) < 1.0E-7:
        Se = (s_tmp - Sr)/(Sm - Sr)
        if (pet*Se*dt - (s_tmp - Sr)) > -1.0E-7:
            evp_tmp = (s_tmp - Sr)/dt
        else:
            evp_tmp = pet*Se
    elif (s_tmp - Sr) > 1.0E-7:
        if (pet*dt-(Sm - Sr)) > -1.0E-7:
            evp_tmp = (Sm - Sr)/dt
        else:
            evp_tmp = pet
    return evp_tmp

def perc(s_tmp, Sm, Sfc, Ks, dt):
    '''
    Percolation function
    '''
    global rp_tmp
    if (s_tmp - Sfc) < 1.0E-7:
        rp_tmp = 0.0
    elif (s_tmp - Sfc) > -1.0E-7 and (s_tmp - Sm) < -1.0E-7:
        Sg = (s_tmp-Sfc)/(Sm-Sfc) # gravitational water
        if (Ks*Sg*dt) - (s_tmp-Sfc) > 1.0E-7:
            rp_tmp = (s_tmp-Sfc)/dt
        else:
            rp_tmp = Ks*Sg
    elif (s_tmp - Sm) > -1.0E-7:
        if (Ks*dt) - (Sm-Sfc) > -1.0E-7:
            rp_tmp = (Sm-Sfc)/dt
        else:
            rp_tmp = Ks
    return rp_tmp

Sm = 450.0
Sfc = 300.0
Sr = 50.0
E = []
T = []
R = []
S = [425.0]
t = 0
t_lst = [0]

while S[t] > (Sr + 0.1):
    print 'S =  %2f at t = %d' % (S[t],t)

    E.append(evp(S[t],Sm, Sr, pet = 3.0, dt = 1))
    S[t] -= E[t]
    T.append(evp(S[t],Sm, Sr, pet = 5.0, dt = 1))
    S[t] -= T[t]

    R.append(perc(S[t],Sm,Sfc, Ks = 10, dt = 1))
    S[t] -= R[t]

    S.append(S[t])
    t += 1
    t_lst.append(t)

R.append(perc(S[t],Sm,Sfc, Ks = 100, dt = 1))
S[t] -= R[t]
E.append(evp(S[t],Sm, Sr, pet = 3.0, dt = 1))
S[t] -= E[t]
T.append(evp(S[t],Sm, Sr, pet = 5.0, dt = 1))
S[t] -= T[t]

fig = plt.figure()    #(8.5,15), dpi=30)
ax1=fig.add_subplot(2,1,1)
ax1.plot(t_lst,S, label = 'S')
plt.legend(loc=0)
plt.grid(True)

ax2=fig.add_subplot(2,1,2)
ax2.plot(t_lst,R, label = 'R')
ax2.plot(t_lst,E, label = 'E')
ax2.plot(t_lst,T, label = 'T')
#ax2.plot(t_lst,Etmpnew, label = 'Enew')
#ax2.plot(t_lst,Ttmpnew, label = 'Tnew')
plt.legend(loc=0)
plt.grid(True)

plt.show()

print '\nDone!'


