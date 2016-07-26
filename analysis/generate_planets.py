# -*- coding: utf-8 -*-
amu = 1.660539040e-27 # kg
Mearth = 5.972e24     # kg
from scipy import interpolate
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import periodictable as pt
import numpy as np
import utils
import pickle
import os

if not os.path.exists('elemental_abundances_sim.pkl'):
    # First, generate or load the solar-system nebula 
    # abundances (normalized to the number of hydrogen atoms):
    if not os.path.exists('elemental_abundances_ssn.pkl'):
        abundances = {}
        Z,name,N = utils.read_abundances('ssabundances_4Gyr.dat')
        # Normalize number abundances with respect to 
        # hydrogen:
        for i in range(len(name)):
            abundances[name[i]] = N[i]
        Habundance = abundances['H']
        for i in range(len(name)):
            abundances[name[i]] = abundances[name[i]]/Habundance
        # Save:
        FILE = open('elemental_abundances_ssn.pkl','w')
        pickle.dump(abundances, FILE)
        FILE.close()

    else:
        FILE = open('elemental_abundances_ssn.pkl','r')
        abundances = pickle.load(FILE)
        FILE.close()
        name = abundances.keys()

    # Get data from the simulations:
    radius,times,results = utils.get_data()
    # Define dictionary to save the elemental abundances:
    elemental_abundances = {}

    elemental_abundances['r'] = radius
    # Generate elemental abundances (by number,
    # in units of hydrogen abundance) dictionary; 
    # this will save, at a given time, the abundance 
    # in gas and solids for each element at each radius:
    for i in range(len(times)):
        elemental_abundances[times[i]] = {}
        for element in name:
            elemental_abundances[times[i]][element+'(g)'] = np.zeros(len(radius))
            elemental_abundances[times[i]][element+'(s)'] = np.zeros(len(radius))
    # Now, at each time, extract elemental abundances in gas and 
    # solids:    
    for i in range(len(times)):
        t = times[i]
        for element in name:
            resulting_abundances = results[t]
            solids = 0.0
            gas = 0.0
            for species in resulting_abundances.keys():
                if utils.checkelement(element, species):
                    # Check how many atoms of the element we have in the species:    
                    fac = utils.howmany(element,species)
                    # Check if solid or gas:
                    if '(a)' in species or '(cr)' in species or '(b)' in species\
                       or '(c)' in species or '(I)' in species or "(I')" in species or\
                       'II' in species:
                        solids = solids + fac*resulting_abundances[species]
                    else:
                        gas = gas + fac*resulting_abundances[species]
            total_element = gas + solids
            # Transform to percentage of gas and solids (in number, normalized to H 
            # number), and save. If total_element is zero (i.e., element not considered in 
            # our model), then give it all to the gas:
            if type(total_element) is not float:
                idx_zeros = np.where(total_element==0)[0]
                if len(idx_zeros)>0:
                    print element
                elemental_abundances[t][element+'(g)'] = (gas/total_element)*abundances[element]
                elemental_abundances[t][element+'(s)'] = (solids/total_element)*abundances[element]
            else:
                elemental_abundances[t][element+'(g)'] = abundances[element]

    FILE = open('elemental_abundances_sim.pkl','w')
    pickle.dump(elemental_abundances, FILE)
    FILE.close()
else:
    FILE = open('elemental_abundances_sim.pkl','r')
    elemental_abundances = pickle.load(FILE)
    FILE.close()

FILE = open('elemental_abundances_ssn.pkl','r')
abundances = pickle.load(FILE)
FILE.close()
name = abundances.keys()

radius = elemental_abundances['r']
times_el = elemental_abundances.keys()
times_el.remove('r')
# Now load planet:
t_acc,a_acc,Mg_acc,Mt_acc,a_width_acc = utils.get_accretion_data('planets/accretion/planet005AU.sal') 
# Transform times to Myr, sort:
t_acc = t_acc/1e6
idx = np.argsort(t_acc)
t_acc = t_acc[idx]
a_acc = a_acc[idx]
Mg_acc = Mg_acc[idx]
Mt_acc = Mt_acc[idx]
a_width_acc = a_width_acc[idx]
# See when it starts accreting gas:
idx = np.where(Mg_acc>0.0)[0]
# Accrete material and keep track of the mass of each element in the 
# envelope. For this, create a dictionary that will save this info.
planet_mass = {}
for element in abundances.keys():
    planet_mass[element+'(g)'] = 0.0
    planet_mass[element+'(s)'] = 0.0
planet_mass_gas = 0.0
planet_mass_solids = 0.0
core_mass = 0.0
# Now go through each time-step:
for i in idx:
    c_t = t_acc[i]
    if core_mass == 0.0:
        core_mass = Mt_acc[i-1]
    # Get how much mass was accreted at this time:
    #print i,c_t
    Mt = Mt_acc[i]-Mt_acc[i-1]
    # How much gas was accreted:
    Mg = Mg_acc[i]-Mg_acc[i-1]
    #print 'Gas accreted:',Mg
    # Get how much solids were accreted:
    Ms = (Mt_acc[i]-Mt_acc[i-1])-Mg
    #print 'Solids accreted:',Ms
    # Interpolate in time, smooth in radius for each element:
    idx_radius = np.where(np.abs(radius-a_acc[i])<(a_width_acc[i]/2.))[0]
    X_g = np.zeros(len(times_el))
    X_s = np.zeros(len(times_el))
    # This keeps track of the denominator in order to get mass fractions
    totg = 0.0
    tots = 0.0
    c_abundances = {}
    for element in abundances.keys():
        if type(elemental_abundances[0.01][element+'(g)']) is not np.float64:
            # Go through each time-step in the chemical simulation, extract the 
            # mean gas and solid abundance of element (which is normalized to the 
            # total number of H atoms) inside 10R_hill of the planet...
            for j in range(len(times_el)):
                X_g[j] = np.mean(elemental_abundances[times_el[j]][element+'(g)'][idx_radius])
                X_s[j] = np.mean(elemental_abundances[times_el[j]][element+'(s)'][idx_radius])
            # Use those abundances with the times at which those abundances occur to create a 
            # function that interpolates the abundances at any time:
            ng = interpolate.interp1d(times_el,X_g)
            ns = interpolate.interp1d(times_el,X_s)
            # Get interpolated abundances at the current time:
            c_abundances[element+'(g)'] = ng(c_t)
            c_abundances[element+'(s)'] = ns(c_t)
            if element == 'CL':
                exec 'totg = totg + ng(c_t)*pt.Cl.mass*amu'
                exec 'tots = tots + ns(c_t)*pt.Cl.mass*amu'
            elif element == 'AL':
                exec 'totg = totg + ng(c_t)*pt.Al.mass*amu'
                exec 'tots = tots + ns(c_t)*pt.Al.mass*amu'
            else:
                exec 'totg = totg + ng(c_t)*pt.'+element+'.mass*amu'
                exec 'tots = tots + ns(c_t)*pt.'+element+'.mass*amu'
        else:
            c_abundances[element+'(g)'] = elemental_abundances[0.01][element+'(g)']
            c_abundances[element+'(s)'] = 0.0
            if element == 'CL':
                exec 'totg = totg + c_abundances["'+element+'(g)"]*pt.Cl.mass*amu'
            elif element == 'AL':
                exec 'totg = totg + c_abundances["'+element+'(g)"]*pt.Al.mass*amu'
            else:
                exec 'totg = totg + c_abundances["'+element+'(g)"]*pt.'+element+'.mass*amu'

    planet_mass_gas = planet_mass_gas + Mg
    if tots != 0.:
        planet_mass_solids = planet_mass_solids + Ms
    for element in abundances.keys():
        if element == 'CL':
            exec "planet_mass[element+'(g)'] = planet_mass[element+'(g)']+(c_abundances[element+'(g)']*pt."+\
                                               "Cl.mass*amu/totg)*Mg"
            if tots != 0.:
                exec "planet_mass[element+'(s)'] = planet_mass[element+'(s)']+(c_abundances[element+'(s)']*pt."+\
                                                   "Cl.mass*amu/tots)*Ms"
        elif element == 'AL':
            exec "planet_mass[element+'(g)'] = planet_mass[element+'(g)']+(c_abundances[element+'(g)']*pt."+\
                                               "Al.mass*amu/totg)*Mg"
            if tots != 0.:
                exec "planet_mass[element+'(s)'] = planet_mass[element+'(s)']+(c_abundances[element+'(s)']*pt."+\
                                                   "Al.mass*amu/tots)*Ms"
        else:
            exec "planet_mass[element+'(g)'] = planet_mass[element+'(g)']+(c_abundances[element+'(g)']*pt."+\
                                               element+".mass*amu/totg)*Mg"
            if tots != 0.:
                exec "planet_mass[element+'(s)'] = planet_mass[element+'(s)']+(c_abundances[element+'(s)']*pt."+\
                                                   element+".mass*amu/tots)*Ms"

#print '\n'
#print '\t ----------------------------------------'
#print 'Total planet mass:',core_mass + planet_mass_gas + planet_mass_solids,'Mearth (expected: ',Mt_acc[-1],'Mearth)'
#print '\nCore:'
#print '-----'
#print '  Core mass',core_mass,'Mearth'
#print '\nEnvelope:'
#print '---------'
#print '  Total gas mass:',planet_mass_gas,'Mearth'
#print '  Total solid mass:',planet_mass_solids,'Mearth\n'
#envelope_mass = planet_mass_gas + planet_mass_solids
#print '  Mass composition of the planetary envelope:'
ztot = 0.0
ztot_gas = 0.0
ztot_solids = 0.0
planet_mass[element] = {}
#print '    X: ',(planet_mass['H(g)'] + planet_mass['H(s)'])/(envelope_mass)
#print '    Y: ',(planet_mass['He(g)']+planet_mass['He(s)'])/(envelope_mass)
for element in abundances.keys():
    if element != 'H' and element != 'He':
        planet_mass[element] = planet_mass[element+'(g)'] + planet_mass[element+'(s)']
        ztot = ztot + planet_mass[element]
        ztot_gas = ztot_gas + planet_mass[element+'(g)']
        ztot_solids = ztot_solids + planet_mass[element+'(s)']
print '\n'
print '\t ----------------------------------------'
print 'Total planet mass:',core_mass + planet_mass_gas + planet_mass_solids,'Mearth (expected: ',Mt_acc[-1],'Mearth)'
print '\nCore:'
print '-----'
print '  Core mass',core_mass,'Mearth'
print '\nEnvelope:'
print '---------'
print '  Total gas mass:',planet_mass_gas,'Mearth'
print '  \t Mass composition:'
print '  \t X: ',planet_mass['H(g)']/planet_mass_gas
print '  \t Y: ',planet_mass['He(g)']/planet_mass_gas
print '  \t Z: ',ztot_gas/planet_mass_gas
print '  Total solid mass:',planet_mass_solids,'Mearth\n'
print '  \t Mass composition:'
print '  \t X: ',planet_mass['H(s)']/planet_mass_solids
print '  \t Y: ',planet_mass['He(s)']/planet_mass_solids
print '  \t Z: ',ztot_solids/planet_mass_solids
envelope_mass = planet_mass_gas + planet_mass_solids
print '  Total mass composition of the planetary envelope:'
planet_mass[element] = {}
print '    X: ',(planet_mass['H(g)'] + planet_mass['H(s)'])/(envelope_mass)
print '    Y: ',(planet_mass['He(g)']+planet_mass['He(s)'])/(envelope_mass)
print '    Z: ',ztot/envelope_mass,'(',(ztot/envelope_mass)/0.0141,'x solar)'
print '\n Total:'
print '\n ------'
print '  Total metallicity of the planet (core + envelope):'
Zenv_final = ztot/envelope_mass
epsilon = core_mass/envelope_mass
Zfin = (Zenv_final + epsilon)/(1.+epsilon)
print '    Z: ',Zfin,'(',(Zfin/0.0141),'x solar)'
print '\n Ratios:'
print '-------'
print '\n Na:'
solar = abundances['Na']/abundances['K']
print '    Na/K',((planet_mass['Na']/pt.Na.mass*amu)/(planet_mass['K']/pt.K.mass*amu))/solar
solar = abundances['Na']/abundances['Fe']
print '    Na/Fe',((planet_mass['Na']/pt.Na.mass*amu)/(planet_mass['Fe']/pt.Fe.mass*amu))/solar
solar = abundances['Na']/abundances['Ti']
print '    Na/Ti',((planet_mass['Na']/pt.Na.mass*amu)/(planet_mass['Ti']/pt.Ti.mass*amu))/solar
solar = abundances['Na']/abundances['V']
print '    Na/V',((planet_mass['Na']/pt.Na.mass*amu)/(planet_mass['V']/pt.V.mass*amu))/solar
solar = abundances['Na']/abundances['O']
print '    Na/O',((planet_mass['Na']/pt.Na.mass*amu)/(planet_mass['O']/pt.O.mass*amu))/solar
solar = abundances['Na']/abundances['C']
print '    Na/C',((planet_mass['Na']/pt.Na.mass*amu)/(planet_mass['C']/pt.C.mass*amu))/solar
solar = abundances['Na']/abundances['N']
print '    Na/N',((planet_mass['Na']/pt.Na.mass*amu)/(planet_mass['N']/pt.N.mass*amu))/solar
print '\n K:'
solar = abundances['K']/abundances['Fe']
print '    K/Fe',((planet_mass['K']/pt.K.mass*amu)/(planet_mass['Fe']/pt.Fe.mass*amu))/solar
solar = abundances['K']/abundances['Ti']
print '    K/Ti',((planet_mass['K']/pt.K.mass*amu)/(planet_mass['Ti']/pt.Ti.mass*amu))/solar
solar = abundances['K']/abundances['V']
print '    K/V',((planet_mass['K']/pt.K.mass*amu)/(planet_mass['V']/pt.V.mass*amu))/solar
solar = abundances['K']/abundances['O']
print '    K/O',((planet_mass['K']/pt.K.mass*amu)/(planet_mass['O']/pt.O.mass*amu))/solar
solar = abundances['K']/abundances['C']
print '    K/C',((planet_mass['K']/pt.K.mass*amu)/(planet_mass['C']/pt.C.mass*amu))/solar
solar = abundances['K']/abundances['N']
print '    K/N',((planet_mass['K']/pt.K.mass*amu)/(planet_mass['N']/pt.N.mass*amu))/solar
print '\n Fe:'
solar = abundances['Fe']/abundances['Ti']
print '    Fe/Ti',((planet_mass['Fe']/pt.Fe.mass*amu)/(planet_mass['Ti']/pt.Ti.mass*amu))/solar
solar = abundances['Fe']/abundances['V']
print '    Fe/V',((planet_mass['Fe']/pt.Fe.mass*amu)/(planet_mass['V']/pt.V.mass*amu))/solar
solar = abundances['Fe']/abundances['O']
print '    Fe/O',((planet_mass['Fe']/pt.Fe.mass*amu)/(planet_mass['O']/pt.O.mass*amu))/solar
solar = abundances['Fe']/abundances['C']
print '    Fe/C',((planet_mass['Fe']/pt.Fe.mass*amu)/(planet_mass['C']/pt.C.mass*amu))/solar
solar = abundances['Fe']/abundances['N']
print '    Fe/N',((planet_mass['Fe']/pt.Fe.mass*amu)/(planet_mass['N']/pt.N.mass*amu))/solar
print '\n Ti:'
solar = abundances['Ti']/abundances['V']
print '    Ti/V',((planet_mass['Ti']/pt.Ti.mass*amu)/(planet_mass['V']/pt.V.mass*amu))/solar
solar = abundances['Ti']/abundances['O']
print '    Ti/O',((planet_mass['Ti']/pt.Ti.mass*amu)/(planet_mass['O']/pt.O.mass*amu))/solar
solar = abundances['Ti']/abundances['C']
print '    Ti/C',((planet_mass['Ti']/pt.Ti.mass*amu)/(planet_mass['C']/pt.C.mass*amu))/solar
solar = abundances['Ti']/abundances['N']
print '    Ti/N',((planet_mass['Ti']/pt.Ti.mass*amu)/(planet_mass['N']/pt.N.mass*amu))/solar
print '\n V:'
solar = abundances['V']/abundances['O']
print '    V/O',((planet_mass['V']/pt.V.mass*amu)/(planet_mass['O']/pt.O.mass*amu))/solar
solar = abundances['V']/abundances['C']
print '    V/C',((planet_mass['V']/pt.V.mass*amu)/(planet_mass['C']/pt.C.mass*amu))/solar
solar = abundances['V']/abundances['N']
print '    V/N',((planet_mass['V']/pt.V.mass*amu)/(planet_mass['N']/pt.N.mass*amu))/solar
print '\n O:'
solar = abundances['O']/abundances['C']
print '    O/C',((planet_mass['O']/pt.O.mass*amu)/(planet_mass['C']/pt.C.mass*amu))/solar
solar = abundances['O']/abundances['N']
print '    O/N',((planet_mass['O']/pt.O.mass*amu)/(planet_mass['N']/pt.N.mass*amu))/solar
print '\n N:'
solar = abundances['N']/abundances['C']
print '    N/C',((planet_mass['N']/pt.N.mass*amu)/(planet_mass['C']/pt.C.mass*amu))/solar
