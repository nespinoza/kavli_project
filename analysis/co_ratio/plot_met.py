import numpy as np
import utilities

r = np.linspace(1.,100.,10000)
T = 200.*((r)**(-0.62))
O_gas = np.zeros(len(r))
CO_ratio_gas = np.zeros(len(r))
CO_ratio_solids = np.zeros(len(r))
for i in range(len(T)):
    temp = T[i]
    Csolids = 0.0
    Osolids = 0.0
    Cgas = 0.0
    Ogas = 0.0
    if temp<=20:
        Csolids = Csolids + 1.5
        Osolids = Osolids + 1.5
    else:
        Cgas = Cgas + 1.5
        Ogas = Ogas + 1.5
    if temp<=47:
        Csolids = Csolids + 0.3
        Osolids = Osolids + 0.6
    else:
        Cgas = Cgas + 0.3
        Ogas = Ogas + 0.6
    if temp<=135:
        Osolids = Osolids + 0.9
    else:
        Ogas = Ogas + 0.9
    if temp<=500:
        Csolids = Csolids + 0.6
    else:
        Cgas = Cgas + 0.6
    if temp<=1500:
        Osolids = Osolids + 1.4
    else:
        Ogas = Ogas + 1.4
    if Ogas != 0:
        O_gas[i] = Ogas
        CO_ratio_gas[i] = Cgas/Ogas
    CO_ratio_solids[i] = Csolids/Osolids
import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.subplot('211')
idx = np.where(CO_ratio_gas!=0)[0]
plt.plot(r[idx],CO_ratio_gas[idx],'-')
plt.plot(r,CO_ratio_solids,'--')
plt.xscale('log')
plt.ylim(0.0,1.2)
#plt.xlabel('Radius (AU)')
plt.ylabel('C/O ratio')

epsilons = [0.5]
Zprotosolar = utilities.realZ()
plt.subplot('212')
idx = np.where(O_gas != 0.0)[0]
nidx = np.where(O_gas == 0.0)[0] 
for epsilon in epsilons:
    Z = np.zeros(len(CO_ratio_gas))
    Zgas = utilities.getZ(CO_ratio_gas,gas = True,no = O_gas)
    Zsolids = utilities.getZ(CO_ratio_solids, solids = True)
    print 'Zgas/Zsolar:    ',Zgas/Zprotosolar
    print 'Zsolids/Zsolar: ',Zsolids/Zprotosolar
    Z[idx] = ((Zgas[idx]+(epsilon*Zsolids[idx]))/(1.+epsilon))
    Z[nidx] = Zsolids[nidx]
    #print Zgas
    #print Zsolids
    #print ((Zgas+(epsilon*Zsolids))/(1.+epsilon))/Zprotosolar
    plt.plot(r,((Zgas+(epsilon*Zsolids))/(1.+epsilon))/Zprotosolar)
plt.xscale('log')
plt.ylabel('$Z_{pl}/Z_*$')
plt.xlabel('Radius (AU)')
plt.show()
