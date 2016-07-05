# -*- coding: utf-8 -*-
import numpy as np
############ Constants #############

# Seconds in a year:
year = 31557600. # s
# Astronomical unit:
AU = 1.496e13 # cm
# Solar mass:
M_sun = 1.98855e33 # grams
# Solar radius:
R_sun = 6.957e10 # cm
# Gravitational Constant:
G = 6.674e-8 # cm^3/g/s^2
# Hydrogen mass:
m_h = 1.6737236e-24 # grams
# Boltzmann's constant:
k_B = 1.3806504e-16 # ergs/K
# Steffan-Boltzmann's constant:
sigma = 5.6704e-5 # erg cm^-2 s^-1 K^-4

####################################
def ChambersModel(r, t, M0 = 0.1, M_s = 1., R_s = 3., T_s = 4200., \
                  s0 = 33., k0 = 3., alpha=1e-2, mu=2.4, gamma=1.4, T_e = 1380.):
    """
    Chambers (2009) disk model

    Given input r (radius, in AU, can be a vector) and 
    t (time, in seconds, can be a vector too), returns 
    the temperature (in K), pressure (in bars) and surface 
    density (in grams/cm^2) of the disk at each radius and time.


    Optional inputs:

        T_s:       Stellar temperature (K).

        M0:        Initial mass of the disk (in solar-mass units).

        R_s:       Stellar radius (in solar-radii units).

        M_s:       Stellar mass (in solar-mass units).

        s0:        Initial radius of the disk's outer edge (AU).

        k0:        Opacity (in cm^2/g).

        alpha:     Viscosity parameter.

        mu:        Mean molecular weight.

        gamma:     Adiabatic index.

        T_e:       Temperature of solid evaporation (K).

    """

    Omega_0 = np.sqrt((G*(M_s*M_sun))/((s0*AU)**3))

    Sigma_vis = (7.*(M0*M_sun))/(10.*np.pi*((s0*AU)**2))

    T_vis = ((27.*k0/(64.*sigma))**(1./3.))*\
        (((alpha*gamma*k_B)/(mu*m_h))**(1./3.))*\
        ((Sigma_vis)**(2./3.))*\
        ((Omega_0)**(1./3.))

    T_rad = ((4./7.)**(1./4.))*\
        (((T_s*k_B*(R_s*R_sun))/(G*(M_s*M_sun)*mu*m_h))**(1./7.))*\
        (((R_s*R_sun)/(s0*AU))**(3./7.))*\
        T_s

    tau_vis = (1./(16*np.pi))*\
          ((mu*m_h)/(alpha*gamma*k_B))*\
          ((Omega_0*(M0*M_sun))/(Sigma_vis*T_vis))

    Sigma_evap = Sigma_vis*((T_vis/T_e)**(14./19.))

    Sigma_rad = Sigma_vis*(T_vis/T_rad)

    # Transition radius between inner and mid-region:
    r_e = (s0)*((Sigma_evap/Sigma_vis)**(95./63.))*\
         (1. + ((t*year)/tau_vis))**(-19./36.)

    # Transition between mid-region and outer region:
    r_o = (s0)*((Sigma_rad/Sigma_vis)**(70./33.))*\
         (1. + ((t*year)/tau_vis))**(-133./132.)

    # Temperature in the outer region:
    T = np.zeros(len(r))
    Sigma = np.zeros(len(r))
    P = np.zeros(len(r))

    idx_in = np.where(r<r_e)[0]
    idx_mid = np.where((r>=r_e)&(r<r_o))[0]
    idx_out = np.where(r>=r_o)[0]

    T[idx_in] = (T_vis**(5./19.))*\
            (T_e**(14./19.))*\
            (r[idx_in]/s0)**(-9./38.)*\
            ((1. + ((t*year)/tau_vis))**(-1./8.))

    Sigma[idx_in] = Sigma_evap*\
                    ((r[idx_in]/s0)**(-24./19.))*\
                    ((1. + ((t*year)/tau_vis)))**(-17./16.)

    
    T[idx_mid] = T_vis*\
             ((r[idx_mid]/s0)**(-9./10.))*\
             ((1. + ((t*year)/tau_vis))**(-19./40.))

    Sigma[idx_mid] = (Sigma_vis)*\
                     ((r[idx_mid]/s0)**(-3./5.))*\
                     ((1. + ((t*year)/tau_vis))**(-57./80.))

    T[idx_out] = T_rad*\
             ((r[idx_out]/s0)**(-3./7.))

    Sigma[idx_out] = Sigma_rad*\
                     ((r[idx_out]/s0)**(-15./14.))*\
                     ((1. + ((t*year)/tau_vis))**(-19./16.))

    P = Sigma*np.sqrt((G*(M_s*M_sun)*k_B*T)/\
                    (2.*np.pi*mu*m_h*((r*AU)**3)))

    return T,P*1e-6,Sigma

def PowerLawProfile(r, q = 0.62, T0 = 200):
    """
    Simple power-law temperature profile for a disk

    Given input r (radius, in AU, can be a vector) it 
    returns temperature (in K) as a function of radius 
    of the form:
  
                      T = T0 * r**(-q)

    The values of q and T0 pre-defined are for the average 
    protoplanetary disk (Andrews & Williams, 2007).

    Optional inputs:

        q:       Power-law index of the profile.

        T0:      Temperature at 1 AU.

    """
    return T0*(r**(-q))
