import matplotlib.pyplot as plt
plt.style.use('ggplot')
import utils
import numpy as np
R = 8.31451

###################################################################################
# Read data, format it for the fit (note it is assumed H comes in kJ/mol:
T_real,Cp_real,S_real,H_real = np.loadtxt('HS.dat',unpack=True,usecols=(0,1,2,4))

# Define...
DeltaH_0 = '9098.00' # H(298.15) - H(0), in J/mol
DeltaH = 142920.00   # Heat of formation in J/mol
###################################################################################

# Convert kJ/mol to J/mol; set H to match heat of formation at 298.15 K:
H_real = H_real*1e3+DeltaH

print '\n'
print '\t Results'
print '\t -------\n'
print '\t Ready-to-copy CEA format:'
utils.printCEAHeader()
# Perform constrained fit for the lower T (<1500 K):
idx_ref = np.where(T_real == 298.15)[0]
Tref = T_real[idx_ref][0]
Cpref = Cp_real[idx_ref][0]
Href = H_real[idx_ref][0]
Sref = S_real[idx_ref][0]
idx = np.where(T_real<=1500)[0]
results_lowT,a1_lowT,a2_lowT,a3_lowT,a4_lowT,a5_lowT,a6_lowT,a7_lowT,b1_lowT,b2_lowT = \
                                     utils.get_coeffs(T_real[idx],\
                                     Cp_real[idx],H_real[idx],S_real[idx],\
                                     Tref, Cpref, Href, Sref)

# Print in ready-to-copy CEA format:
utils.printCEA(np.min(T_real[idx]),np.max(T_real[idx]),a1_lowT,a2_lowT,a3_lowT,a4_lowT,\
               a5_lowT,a6_lowT,a7_lowT,b1_lowT,b2_lowT,DeltaH_0,DeltaH)

# Now do the same for mid-T (1500 K < T < 3000 K 
Tref = 1500.
Cpref = utils.modelCpR(1500.,results_lowT.params)*R
Href = utils.modelHR(1500.,results_lowT.params)*R
Sref = utils.modelSR(1500.,results_lowT.params)*R
idx = np.where((T_real>1500.)&(T_real<=3000.))[0]
results_midT,a1_midT,a2_midT,a3_midT,a4_midT,a5_midT,a6_midT,a7_midT,b1_midT,b2_midT = \
                                     utils.get_coeffs(T_real[idx],\
                                     Cp_real[idx],H_real[idx],S_real[idx],\
                                     Tref, Cpref, Href, Sref)

# Print in ready-to-copy CEA format:
utils.printCEA(np.min(T_real[idx-1]),np.max(T_real[idx]),a1_midT,a2_midT,a3_midT,a4_midT,\
               a5_midT,a6_midT,a7_midT,b1_midT,b2_midT,DeltaH_0,DeltaH)

# Finally, same thing for high-T (3000 K > T)
Tref = 3000.
Cpref = utils.modelCpR(3000.,results_midT.params)*R
Href = utils.modelHR(3000.,results_midT.params)*R
Sref = utils.modelSR(3000.,results_midT.params)*R
idx = np.where(T_real>3000.)[0]
results_highT,a1_highT,a2_highT,a3_highT,a4_highT,a5_highT,a6_highT,a7_highT,b1_highT,b2_highT = \
                                     utils.get_coeffs(T_real[idx],\
                                     Cp_real[idx],H_real[idx],S_real[idx],\
                                     Tref, Cpref, Href, Sref)
# Print in ready-to-copy CEA format:
utils.printCEA(np.min(T_real[idx-1]),np.max(T_real[idx]),a1_highT,a2_highT,a3_highT,a4_highT,\
               a5_highT,a6_highT,a7_highT,b1_highT,b2_highT,DeltaH_0,DeltaH)

# Plot the fits:
Tlow = np.linspace(np.min(T_real),1500.,100)
Tmid = np.linspace(1500.,3000.,100)
Thigh = np.linspace(3000.,np.max(T_real),100)

plt.subplot('311')
plt.title('Element fit')
plt.plot(T_real,Cp_real/R,'wo',label = 'Data')
plt.plot(Tlow,utils.modelCpR(Tlow,results_lowT.params),linewidth=2,label = 'Constrained fit, low T')
plt.plot(Tmid,utils.modelCpR(Tmid,results_midT.params),linewidth=2,label = 'Constrained fit, mid T')
plt.plot(Thigh,utils.modelCpR(Thigh,results_highT.params),linewidth=2,label = 'Constrained fit, high T')
plt.ylabel('$C_p^o(T)/R$')
plt.xlim([200,6100])
plt.ylim([np.min(Cp_real/R)-0.1,np.max(Cp_real/R)+0.1])
plt.legend(loc='best')
plt.subplot('312')
plt.plot(T_real,S_real/R,'wo')
plt.plot(Tlow,utils.modelSR(Tlow,results_lowT.params),linewidth=2)
plt.plot(Tmid,utils.modelSR(Tmid,results_midT.params),linewidth=2)
plt.plot(Thigh,utils.modelSR(Thigh,results_highT.params),linewidth=2)
plt.ylabel('$S^o(T)/R$')
plt.xlim([200,6100])
plt.ylim([np.min(S_real/R)-1,np.max(S_real/R)+1])
plt.subplot('313')
plt.plot(T_real,H_real/R,'wo')
plt.plot(Tlow,utils.modelHR(Tlow,results_lowT.params),linewidth=2)
plt.plot(Tmid,utils.modelHR(Tmid,results_midT.params),linewidth=2)
plt.plot(Thigh,utils.modelHR(Thigh,results_highT.params),linewidth=2)
plt.ylabel('$H^o(T)/R$')
plt.xlabel('Temperature (K)')
plt.xlim([200,6100])
plt.ylim([np.min(H_real/R)-500,np.max(H_real/R)+500])
plt.show()
