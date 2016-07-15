import matplotlib.pyplot as plt
plt.style.use('ggplot')
import utils
import numpy as np
R = 8.31451

# Read data, format it for the fit
T_real,Cp_real,S_real,H_real = np.loadtxt('al2o_plus.dat',unpack=True)
DeltaH_0 = '12981.148' # H(298.15) - H(0)
DeltaH = 648970.248 # Heat of formation for AL2O+ according to CEA
H_real = H_real*1e3+DeltaH # H in Chase tables are on kJ/mol

# Perform constrained fit:
results,a1,a2,a3,a4,a5,a6,a7,b1,b2 = utils.get_coeffs(T_real,\
                                     Cp_real,H_real,S_real,\
                                     T_real[0],Cp_real[0],\
                                     H_real[0],S_real[0])

# Coefficients from CEA tables:
a_cea = [6.828925720E+04,-9.098504170E+02,8.896563010E+00,-7.772554380E-04,\
         -4.034655620E-07, 6.977315050E-10,-2.417262484E-13]

b_cea = [8.085007550E+04,-2.176216731E+01]

# Plot
T = np.linspace(np.min(T_real),np.max(T_real),100)

# CEA fits:
Cp_CEA = (a_cea[0]/T**2) + (a_cea[1]/T) + a_cea[2] + (a_cea[3]*T) + a_cea[4]*T**2+\
                a_cea[5]*T**3 + a_cea[6]*T**4

H_CEA = (-a_cea[0]/T) + (a_cea[1]*np.log(T)) + a_cea[2]*T + (a_cea[3]/2.)*T**2 + \
                (a_cea[4]/3.)*T**3+(a_cea[5]/4.)*T**4 + (a_cea[6]/5.)*T**5 + b_cea[0]

S_CEA = (-a_cea[0]/(2.*T**2)) - (a_cea[1]/T) + a_cea[2]*np.log(T) + (a_cea[3]*T) + \
                (a_cea[4]/2.)*T**2+(a_cea[5]/3.)*T**3 + (a_cea[6]/4.)*T**4 + b_cea[1]

plt.subplot('311')
plt.title('Al$_2$O$^+$')
plt.plot(T_real,Cp_real/R,'o',label = 'Data from Chase (1998)')
plt.plot(T,utils.modelCpR(T,results.params),label = 'Constrained fit (this work)')
plt.plot(T,Cp_CEA,linewidth=5,alpha=0.3,label = 'CEA table fit')
plt.ylabel('$C_p^o(T)/R$')
plt.xlim([200,1100])
plt.ylim([np.min(Cp_real/R)-0.1,np.max(Cp_real/R)+0.1])
plt.legend(loc='best')
plt.subplot('312')
plt.plot(T_real,S_real/R,'o')
plt.plot(T,S_CEA,linewidth=5,alpha=0.3)
plt.plot(T,utils.modelSR(T,results.params))
plt.ylabel('$S^o(T)/R$')
plt.xlim([200,1100])
plt.ylim([np.min(S_real/R)-1,np.max(S_real/R)+1])
plt.subplot('313')
plt.plot(T_real,H_real/R,'o')
plt.plot(T,utils.modelHR(T,results.params))
plt.plot(T,H_CEA,linewidth=5,alpha=0.3)
plt.ylabel('$H^o(T)/R$')
plt.xlabel('Temperature (K)')
plt.xlim([200,1100])
plt.ylim([np.min(H_real/R)-500,np.max(H_real/R)+500])
plt.show()

print 'CEA-like format:'

line1 = '    {0:4.3f}  {1:4.3f} 7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0'.format(np.min(T_real),np.max(T_real))
nspaces1 = 80-(len(line1)+len(DeltaH_0))
line1 = line1+''.join(nspaces1*[' '])+DeltaH_0
print line1
line2 = ''
for a in [a1,a2,a3,a4,a5]:
    if a>=0:
        line2 = line2+(' {0:.9e}'.format(a)).replace('e','D')
    else:
        line2 = line2+('{0:.9e}'.format(a)).replace('e','D')
print line2
line3 = ''
for a in [a6,a7,0.0,b1,b2]:
    if a>=0:
        line3 = line3+(' {0:.9e}'.format(a)).replace('e','D')
    else:
        line3 = line3+('{0:.9e}'.format(a)).replace('e','D')
print line3
