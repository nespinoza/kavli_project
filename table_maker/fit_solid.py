import matplotlib.pyplot as plt
plt.style.use('ggplot')
import utils
import numpy as np
R = 8.31451

###################################################################################
# Read data, format it for the fit (note it is assumed H comes in kJ/mol:
T_real,Cp_real,S_real,H_real = np.loadtxt('CaAlO.dat',unpack=True,usecols=(0,1,2,4))

# Define...
DeltaH_0 = '0.00'      # H(298.15) - H(0), in J/mol
DeltaH = -10703217.0   # Heat of formation in J/mol
###################################################################################

# Convert kJ/mol to J/mol; set H to match heat of formation at 298.15 K:
H_real = H_real*1e3+DeltaH

print '\n'
print '\t Results'
print '\t -------\n'
print '\t Ready-to-copy CEA format:'
utils.printCEAHeader()
# Perform constrained fit:
idx_ref = np.where(T_real == 298.15)[0]
Tref = T_real[idx_ref][0]
Cpref = Cp_real[idx_ref][0]
Href = H_real[idx_ref][0]
Sref = S_real[idx_ref][0]
T_transition = 2100.0
results,a1_1,a2_1,a3_1,a4_1,a5_1,a6_1,a7_1,b1_1,b2_1,a0,b1_2,b2_2 = \
                                     utils.get_coeffs_phase_transition(T_real,\
                                     Cp_real,H_real,S_real,\
                                     Tref, Cpref, Href, Sref,T_transition)

# Print in ready-to-copy CEA format:
utils.printCEA(np.min(T_real),T_transition,a1_1,a2_1,a3_1,a4_1,\
               a5_1,a6_1,a7_1,b1_1,b2_1,DeltaH_0,DeltaH)

utils.printCEAHeader()
utils.printCEA(T_transition,np.max(T_real),0.0,0.0,a0,0.0,\
               0.0,0.0,0.0,b1_2,b2_2,DeltaH_0,DeltaH)

# Plot the fits:
Tlow = np.linspace(np.min(T_real),T_transition,100)
Tmid = np.linspace(T_transition,np.max(T_real),100)

plt.subplot('311')
plt.title('Element fit')
plt.plot(T_real,Cp_real/R,'wo',label = 'Data')
plt.plot(Tlow,utils.modelCpR_PT(Tlow,results.params,T_transition),linewidth=2,label = 'Constrained fit, low T')
plt.plot(Tmid,utils.modelCpR_PT(Tmid,results.params,T_transition),linewidth=2,label = 'Constrained fit, mid T')
plt.ylabel('$C_p^o(T)/R$')
plt.xlim([200,3100])
plt.ylim([np.min(Cp_real/R)-0.1,np.max(Cp_real/R)+0.1])
plt.legend(loc='best')
plt.subplot('312')
plt.plot(T_real,S_real/R,'wo')
plt.plot(Tlow,utils.modelSR_PT(Tlow,results.params,T_transition),linewidth=2)
plt.plot(Tmid,utils.modelSR_PT(Tmid,results.params,T_transition),linewidth=2)
plt.ylabel('$S^o(T)/R$')
plt.xlim([200,3100])
plt.ylim([np.min(S_real/R)-1,np.max(S_real/R)+1])
plt.subplot('313')
plt.plot(T_real,H_real/R,'wo')
plt.plot(Tlow,utils.modelHR_PT(Tlow,results.params,T_transition),linewidth=2)
plt.plot(Tmid,utils.modelHR_PT(Tmid,results.params,T_transition),linewidth=2)
plt.ylabel('$H^o(T)/R$')
plt.xlabel('Temperature (K)')
plt.xlim([200,3100])
plt.ylim([np.min(H_real/R)-500,np.max(H_real/R)+500])
plt.show()
