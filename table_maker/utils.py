#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import lmfit
# Gas constant:
R = 8.31451

def getb1b2(exponents, coeffs, HF, S, T = 298.15):
    """

    Given the exponents and the coefficients of a fit to the 
    non-dimensional specific heat, Cp/R, of the form:

             Cp/R = A1*(T**a1) + A2*(T**a2) + ...

    where exponents = [a1,a2,...] and coeffs = [A1,A2,...], 
    along with the heat of formation (in Joules/mol) and the 
    entropy at the same reference temperature, this 
    returns the b1 and b2 coefficients needed for the thermo.inp 
    CEA file. 

    For more info see: 

     http://shepherd.caltech.edu/EDL/public/formats/nasa.html

    Input:

        exponents         a1,a2,... (exponents of the Cp/R fit).

        coefficients      A1,A2,... (coefficients of the Cp/R fit).

        HF                Heat of formation (J/mol)

        S                 Entropy at the reference temperature (T)

    Output:

        b1                Integration constant for the enthalpy.

        b2                Integration constant for the entropy.

    """

    HR = 0.0
    # First get the terms in H/R = int (Cp/R dT):
    for i in range(len(exponents)):
        if exponents[i] == -1:
            HR = HR + coeffs[i]*np.log(T)
        else:
            C = exponents[i]+1.
            HR = HR + coeffs[i]*(1./C)*(T**C)

    # Now HF = R*b1 + int Cp = R*b1 + R*int(Cp/R), 
    # so find for b1 and voila:
    b1 = (HF/R)-HR

    # Now get the terms in S/R = int (Cp/(R*T)) dT
    SR = 0.0
    for i in range(len(exponents)):
        if exponents[i] == 0.0:
            SR = SR + coeffs[i]*np.log(T)
        else:
            C = exponents[i]
            SR = SR + coeffs[i]*(1./C)*(T**C)
    b2 = (S/R)-SR
    return b1,b2

def modelCpR(T,params):
    return (params['a1'].value/T**2)+\
           (params['a2'].value/T)+\
           params['a3'].value+\
           params['a4'].value*T+\
           params['a5'].value*(T**2)+\
           params['a6'].value*(T**3)+\
           params['a7'].value*(T**4)

def modelHR(T,params):
    return -(params['a1'].value/T)+\
           (params['a2'].value*np.log(T))+\
           params['a3'].value*T+\
           (params['a4'].value/2.)*T**2+\
           (params['a5'].value/3.)*(T**3)+\
           (params['a6'].value/4.)*(T**4)+\
           (params['a7'].value/5.)*(T**5)+\
           params['b1'].value

def modelSR(T,params):
    return -(params['a1'].value/(2.*T**2))-\
           (params['a2'].value/T)+\
           params['a3'].value*np.log(T)+\
           params['a4'].value*T+\
           (params['a5'].value/2.)*(T**2)+\
           (params['a6'].value/3.)*(T**3)+\
           (params['a7'].value/4.)*(T**4)+\
           params['b2'].value

def get_coeffs(T, Cp, H, S, Tr, Cpr, Hr, Sr):
    """
    This function fits simultaneously Cp/R, H/R 
    and S/R given a temperature T, specific heat Cp, 
    enthalpy H, and entropy S table. The fit is of the form:

    Cp/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4

    H(T)/(R) = -a1*T^-1 + a2*ln(T) + a3*T + a4*T**2/2 + a5*T^3/3 +
               a6*T^4/4 + a7*T^5/5 + b1 

    S(T)/R = -a1*T^-2/2 - a2*T^-1 + a3*ln(T) + a4*T + a5*T^2/2 +
               a6*T^3/3 + a7*T^4/4 + b2

    The program also needs a reference temperature (Tr), 
    specific heat (Cpr), enthalpy (Hr) and entropy (Sr). 
    The fit is then constrained to match the reference 
    specific heat, enthalpy and entropy exactly.

    For more info see: http://www.grc.nasa.gov/WWW/CEAWeb/RP-1311.pdf

    Input:

        T          Temperatures to fit.
 
        Cp         Specific heats to fit (in J/K/mol).

        H          Enthalpies to fit (in J/mol; note that the absolute 
                   enthalpies must be entered, e.g., usually the enthalpy 
                   at 298.15 is the heat of formation).

        S          Entropy to fit (in J/K/mol).

        Tr         Reference temperature

        CpR        Reference specific heat

        Hr         Reference enthalpy.

        Sr         Reference entropy.

    Output

        result,a1,a2,a3,a4,a5,a6,a7,b1,b2

    Here, result is the lmfit result, the a_i the coefficients of the fitted functions and 
    b_i the integration constants for enthalpy (b1) and the entropy (b2).
    """

    Ndata = len(T)

    def residuals(params, X, Y):
        return np.sqrt((Y[:,0] - modelCpR(X,params))**2 + \
                       (Y[:,1] - modelHR(X,params))**2 + \
                       (Y[:,2] - modelSR(X,params))**2)

    # Initialize the fit. For this, fit Cp via common 
    # least-squares to get a1,a2,...,a7. Write the design
    # matrix, solve the system and voila:
    X = np.ones([Ndata,7]) 
    X[:,0] = (1./T**2) 
    X[:,1] = (1./T)
    X[:,3] = T
    X[:,4] = T**2
    X[:,5] = T**3
    X[:,6] = T**4
 
    ga1,ga2,ga3,ga4,ga5,ga6,ga7 = np.dot(np.dot(np.linalg.inv(np.dot(X.transpose(),X)),\
                                                X.transpose()),Cp/R)

    # Init lmfit
    prms = lmfit.Parameters()
    prms.add('a1',value = ga1)
    prms.add('a2',value = ga2)
    prms.add('a3',value = ga3)
    prms.add('a4',value = ga4)
    prms.add('a5',value = ga5)
    prms.add('a6',value = ga6)

    # Constrain for the fit to match Cpr/R:
    #prms.add('a7',value = ga7)
    prms.add('a7', expr='('+str(Cpr/R)+'-(a1/('+str(Tr)+'**2)+a2/('+str(Tr)+') + a3'+\
                        '+a4*'+str(Tr)+'+a5*'+str(Tr)+'**2 + a6*'+str(Tr)+'**3))/('+str(Tr)+'**4)')

    # Constrain for the fit to match Hr/R: 
    prms.add('b1', expr=str(Hr/R)+'+(a1/'+str(Tr)+')-a2*log('+str(Tr)+')'+\
                        '-a3*'+str(Tr)+'-(a4/2)*('+str(Tr)+'**2)'+\
                        '-(a5/3)*('+str(Tr)+'**3)'+'-(a6/4)*('+str(Tr)+'**4)'+\
                        '-(a7/5)*('+str(Tr)+'**5)')

    # Constrain for the fit to match Sr/R:
    #prms.add('b2',value = Sr/R - (-ga1/(2.*Tr**2)-(ga2/Tr)+ga3*np.log(Tr)+ga4*Tr+\
    #                             (ga5/2)*Tr**2 + (ga6/3.)*Tr**3 + (ga7/4.)*Tr**4))
    prms.add('b2', expr=str(Sr/R)+'+a1/(2.*'+str(Tr)+'**2)+a2/('+str(Tr)+')'+\
                        '-a3*log('+str(Tr)+')-(a4*'+str(Tr)+')'+\
                        '-(a5/2)*('+str(Tr)+'**2)'+'-(a6/3)*('+str(Tr)+'**3)'+\
                        '-(a7/4)*('+str(Tr)+'**4)')

    # Prepare data:
    Y = np.zeros([Ndata,3])
    Y[:,0] = Cp/R
    Y[:,1] = H/R
    Y[:,2] = S/R

    # Perform constrained fit:
    result = lmfit.minimize(residuals, prms, args=(T,Y))
    
    # If fit is successful, return fitted values. If not, cry:
    if result.success:
        return result, result.params['a1'].value, result.params['a2'].value, \
               result.params['a3'].value, result.params['a4'].value, \
               result.params['a5'].value, result.params['a6'].value,\
               result.params['a7'].value, result.params['b1'].value,\
               result.params['b2'].value
    else:
        print 'Fit failed. Too bad :(.'
        return result,0,0,0,0,0,0,0,0,0
