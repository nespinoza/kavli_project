import numpy as np
import disk_models
import matplotlib.pyplot as plt
plt.style.use('ggplot')

Myr = 1e6
r = 10**(np.linspace(-2,1,100))
times = [0.1,1.,3.,6.]
plt.subplot(311)
for t in times:
    T,P,Sigma = disk_models.ChambersModel(r,t*Myr)
    plt.plot(r,T,label=str(t)+' Myr')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Temperature (K)')
plt.legend()
plt.subplot(312)
for t in times:
    T,P,Sigma = disk_models.ChambersModel(r,t*Myr)
    plt.plot(r,P,label=None)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Pressure (bars)')
plt.subplot(313)
for t in times:
    T,P,Sigma = disk_models.ChambersModel(r,t*Myr)
    plt.plot(r,Sigma,label=None)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Surface Density (g/cm$^2$)')
plt.xlabel('Radius (AU)')
plt.show()
