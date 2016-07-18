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

def printCEAHeader():
    print 'SN'+''.join(22*[' '])+'K. Lodders, B. Fegley, priv. comm.'
    print ' 3 NEF/16 '

def printCEA(Tmin,Tmax,a1,a2,a3,a4,a5,a6,a7,b1,b2,DeltaH_0,DeltaH):
    if Tmin<1000:
        line1 = '    {0:4.3f}  {1:4.3f} 7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0'.format(Tmin,Tmax)
    else:
        line1 = '   {0:4.3f}  {1:4.3f} 7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0'.format(Tmin,Tmax)
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

def modelCpR(T,params):
    if not (type(params) is list or type(params) is np.ndarray):
        return (params['a1'].value/T**2)+\
               (params['a2'].value/T)+\
               params['a3'].value+\
               params['a4'].value*T+\
               params['a5'].value*(T**2)+\
               params['a6'].value*(T**3)+\
               params['a7'].value*(T**4)
    else:
        return (params[0]/T**2)+\
               (params[1]/T)+\
               params[2]+\
               params[3]*T+\
               params[4]*(T**2)+\
               params[5]*(T**3)+\
               params[6]*(T**4)

def modelCpR_PT(T,params,T_transition,idx1 = None, idx2 = None, idx3 = None, model = '7+constant'):
    if (idx1 is None) and (type(T_transition) is not list):
       idx1 = np.where(T<=T_transition)[0]
       idx2 = np.where(T>T_transition)[0]
    elif idx1 is None:
       idx1 = np.where(T<=T_transition[0])[0]
       idx2 = np.where((T>T_transition[0])&(T<=T_transition[1]))[0]
       idx3 = np.where(T>T_transition[1])[0]

    result = np.ones(len(T))
    if not (type(params) is list or type(params) is np.ndarray):
        result[idx1] = (params['a1_1'].value/T[idx1]**2)+\
               (params['a2_1'].value/T[idx1])+\
               params['a3_1'].value+\
               params['a4_1'].value*T[idx1]+\
               params['a5_1'].value*(T[idx1]**2)+\
               params['a6_1'].value*(T[idx1]**3)+\
               params['a7_1'].value*(T[idx1]**4)
        if model == '7+constant':
            result[idx2] = result[idx2]*params['a0'].value
        elif model == '7+linear':
            result[idx2] = params['b0'].value*T[idx2] + params['a0'].value
        elif model == '7+constant+constant':
            result[idx2] = result[idx2]*params['a0'].value
            result[idx3] = result[idx3]*params['a0_2'].value
    else:
        result[idx1] = (params[0]/T[idx1]**2)+\
               (params[1]/T[idx1])+\
               params[2]+\
               params[3]*T[idx1]+\
               params[4]*(T[idx1]**2)+\
               params[5]*(T[idx1]**3)+\
               params[6]*(T[idx1]**4)
        if model == '7+constant':
            result[idx2] = result[idx2]*params[9]
        elif model == '7+linear':
            result[idx2] = params[10]*T[idx2] + params[9]
        elif model == '7+constant+constant':
            result[idx2] = result[idx2]*params[9]
            result[idx3] = result[idx3]*params[12]
    return result

def modelHR(T,params):
    if not (type(params) is list or type(params) is np.ndarray):
        return -(params['a1'].value/T)+\
           (params['a2'].value*np.log(T))+\
           params['a3'].value*T+\
           (params['a4'].value/2.)*T**2+\
           (params['a5'].value/3.)*(T**3)+\
           (params['a6'].value/4.)*(T**4)+\
           (params['a7'].value/5.)*(T**5)+\
           params['b1'].value
    else:
        return -(params[0]/T)+\
           (params[1]*np.log(T))+\
           params[2]*T+\
           (params[3]/2.)*T**2+\
           (params[4]/3.)*(T**3)+\
           (params[5]/4.)*(T**4)+\
           (params[6]/5.)*(T**5)+\
           params[7]

def modelHR_PT(T,params,T_transition,idx1 = None, idx2 = None, idx3 = None, model = '7+constant'):

    if (idx1 is None) and (type(T_transition) is not list):
       idx1 = np.where(T<=T_transition)[0]
       idx2 = np.where(T>T_transition)[0]
    elif idx1 is None:
       idx1 = np.where(T<=T_transition[0])[0]
       idx2 = np.where((T>T_transition[0])&(T<=T_transition[1]))[0]
       idx3 = np.where(T>T_transition[1])[0]

    result = np.ones(len(T))
    if not (type(params) is list or type(params) is np.ndarray):
        result[idx1] = -(params['a1_1'].value/T[idx1])+\
           (params['a2_1'].value*np.log(T[idx1]))+\
           params['a3_1'].value*T[idx1]+\
           (params['a4_1'].value/2.)*T[idx1]**2+\
           (params['a5_1'].value/3.)*(T[idx1]**3)+\
           (params['a6_1'].value/4.)*(T[idx1]**4)+\
           (params['a7_1'].value/5.)*(T[idx1]**5)+\
           params['b1_1'].value
        if model == '7+constant':
            result[idx2] = params['a0'].value*T[idx2] + params['b1_2'].value
        elif model == '7+linear':
            result[idx2] = params['a0'].value*T[idx2] + params['b0'].value*(T[idx2]**2)/2. + \
                           params['b1_2'].value
        elif model == '7+constant+constant':
            result[idx2] = params['a0'].value*T[idx2] + params['b1_2'].value
            result[idx3] = params['a0_2'].value*T[idx3] + params['b1_3'].value
    else:
        result[idx1] =  -(params[0]/T[idx1])+\
           (params[1]*np.log(T[idx1]))+\
           params[2]*T[idx1]+\
           (params[3]/2.)*T[idx1]**2+\
           (params[4]/3.)*(T[idx1]**3)+\
           (params[5]/4.)*(T[idx1]**4)+\
           (params[6]/5.)*(T[idx1]**5)+\
           params[7]
        if model == '7+constant':
            result[idx2] = params[9]*T[idx2] + params[10]
        elif model == '7+linear':
            result[idx2] = params[9]*T[idx2] + params[10]*(T[idx2]**2)/2. + params[11]
        elif model == '7+constant+constant':
            result[idx2] = params[9]*T[idx2] + params[10]
            result[idx3] = params[12]*T[idx3] + params[13]
    return result

def modelSR(T,params):
    if not (type(params) is list or type(params) is np.ndarray):
        return -(params['a1'].value/(2.*T**2))-\
           (params['a2'].value/T)+\
           params['a3'].value*np.log(T)+\
           params['a4'].value*T+\
           (params['a5'].value/2.)*(T**2)+\
           (params['a6'].value/3.)*(T**3)+\
           (params['a7'].value/4.)*(T**4)+\
           params['b2'].value

    else:
        return -(params[0]/(2.*T**2))-\
           (params[1]/T)+\
           params[2]*np.log(T)+\
           params[3]*T+\
           (params[4]/2.)*(T**2)+\
           (params[5]/3.)*(T**3)+\
           (params[6]/4.)*(T**4)+\
           params[8]

def modelSR_PT(T,params,T_transition,idx1 = None, idx2 = None, idx3 = None, model = '7+constant'):
    if (idx1 is None) and (type(T_transition) is not list):
       idx1 = np.where(T<=T_transition)[0]
       idx2 = np.where(T>T_transition)[0]
    elif idx1 is None:
       idx1 = np.where(T<=T_transition[0])[0]
       idx2 = np.where((T>T_transition[0])&(T<=T_transition[1]))[0]
       idx3 = np.where(T>T_transition[1])[0]

    result = np.ones(len(T))
    if not (type(params) is list or type(params) is np.ndarray):
        result[idx1] = -(params['a1_1'].value/(2.*T[idx1]**2))-\
           (params['a2_1'].value/T[idx1])+\
           params['a3_1'].value*np.log(T[idx1])+\
           params['a4_1'].value*T[idx1]+\
           (params['a5_1'].value/2.)*(T[idx1]**2)+\
           (params['a6_1'].value/3.)*(T[idx1]**3)+\
           (params['a7_1'].value/4.)*(T[idx1]**4)+\
           params['b2_1'].value
        if model == '7+constant':
            result[idx2] = params['a0'].value*np.log(T[idx2]) + params['b2_2'].value
        elif model == '7+linear':
            result[idx2] = params['a0'].value*np.log(T[idx2]) + params['b0'].value*T[idx2]+\
                           params['b2_2'].value
        elif model == '7+constant+constant':
            result[idx2] = params['a0'].value*np.log(T[idx2]) + params['b2_2'].value
            result[idx3] = params['a0_2'].value*np.log(T[idx3]) + params['b2_3'].value

    else:
        result[idx1] = -(params[0]/(2.*T[idx1]**2))-\
           (params[1]/T[idx1])+\
           params[2]*np.log(T[idx1])+\
           params[3]*T[idx1]+\
           (params[4]/2.)*(T[idx1]**2)+\
           (params[5]/3.)*(T[idx1]**3)+\
           (params[6]/4.)*(T[idx1]**4)+\
           params[8]
        if model == '7+constant':
            result[idx2] = params[9]*np.log(T[idx2]) + params[11]
        elif model == '7+linear':
            result[idx2] = params[9]*np.log(T[idx2]) + params[10]*T[idx2] + params[12]
        elif model == '7+constant+constant':
            result[idx2] = params[9]*np.log(T[idx2]) + params[11]
            result[idx3] = params[12]*np.log(T[idx3]) + params[14]
    return result

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

def get_coeffs_phase_transition(T, Cp, H, S, Tr, Cpr, Hr, Sr, T_transition, model = '7+constant'):
    """
    This function fits simultaneously Cp/R, H/R 
    and S/R given a temperature T, specific heat Cp, 
    enthalpy H, and entropy S table in the case when there 
    is a phase transition. The fit is similar to the 
    get_coeffs function, although with different coefficients 
    between each phase. For each phase, the fit is of the form:

    Cp/R = a1*T^-2 + a2*T^-1 + a3 + a4*T + a5*T^2 + a6*T^3 + a7*T^4

    H(T)/(R) = -a1*T^-1 + a2*ln(T) + a3*T + a4*T**2/2 + a5*T^3/3 +
               a6*T^4/4 + a7*T^5/5 + b1 

    S(T)/R = -a1*T^-2/2 - a2*T^-1 + a3*ln(T) + a4*T + a5*T^2/2 +
               a6*T^3/3 + a7*T^4/4 + b2

    The fit is constrained to (1) pass exactly at the reference temperature 
    (Tr), specific heat (Cpr), enthalpy (Hr) and entropy (Sr), which is usually 
    these thermodynamic quantities at 298.15 (in order to reproduce heats of formation 
    exactly) and (2) to give zero Gibbs free energy between the phases at temperature 
    T_transition. 

    For more info see: http://www.grc.nasa.gov/WWW/CEAWeb/RP-1311.pdf

    Input:

        T              Temperatures to fit.
 
        Cp             Specific heats to fit (in J/K/mol).

        H              Enthalpies to fit (in J/mol; note that the absolute 
                       enthalpies must be entered, e.g., usually the enthalpy 
                       at 298.15 is the heat of formation).

        S              Entropy to fit (in J/K/mol).

        Tr             Reference temperature

        CpR            Reference specific heat

        Hr             Reference enthalpy.

        Sr             Reference entropy.

        T_transition   Temperature at which phase transition occurs.

    Output

        result,[a1,a2,a3,a4,a5,a6,a7,b1,b2]^1,[a1,a2,a3,a4,a5,a6,a7,b1,b2]^2

    Here, result is the lmfit result, the a_i the coefficients of the fitted functions and 
    b_i the integration constants for enthalpy (b1) and the entropy (b2) for the phase before 
    T_Transition ([]^1) and after phase transition ([]^2).
    """

    if type(T_transition) is not list:
        idx1 = np.where(T<=T_transition)[0]
        idx2 = np.where(T>T_transition)[0]
        idx3 = None
        Ndata1 = len(idx1)
        Ndata2 = len(idx2)
        Ndata3 = None
    else:
        idx1 = np.where(T<=T_transition[0])[0]
        idx2 = np.where((T>T_transition[0])&(T<=T_transition[1]))[0]
        idx3 = np.where(T>T_transition[1])[0]
        Ndata1 = len(idx1)
        Ndata2 = len(idx2)
        Ndata3 = len(idx3)

    def residuals(params, X, Y):
        return np.sqrt((Y[:,0] - modelCpR_PT(X,params,T_transition,idx1=idx1,idx2=idx2, \
                                             idx3=idx3, model = model))**2 + \
                       (Y[:,1] - modelHR_PT(X,params,T_transition,idx1=idx1,idx2=idx2,\
                                             idx3=idx3, model = model))**2 + \
                       (Y[:,2] - modelSR_PT(X,params,T_transition,idx1=idx1,idx2=idx2,\
                                             idx3=idx3, model = model))**2)

    # Initialize the fit. For this, fit Cp via common 
    # least-squares to get a1,a2,...,a7 at each side of the 
    # Cp curve (before and after the transition). Write the design
    # matrix, solve the system and voila:
    X = np.ones([Ndata1,7])
    X[:,0] = (1./T[idx1]**2) 
    X[:,1] = (1./T[idx1])
    X[:,3] = T[idx1]
    X[:,4] = T[idx1]**2
    X[:,5] = T[idx1]**3
    X[:,6] = T[idx1]**4

    ga1_1,ga2_1,ga3_1,ga4_1,ga5_1,ga6_1,ga7_1 = np.dot(np.dot(np.linalg.inv(np.dot(X.transpose(),X)),\
                                                X.transpose()),Cp[idx1]/R)

    if model == '7+constant':
        X = np.ones([Ndata2,2])
        X[:,0] = np.log(T[idx2])

        ga0,b2_2 = np.dot(np.dot(np.linalg.inv(np.dot(X.transpose(),X)),\
                          X.transpose()),S[idx2]/R)
    elif model == '7+linear':
        X = np.ones([Ndata2,3])
        X[:,0] = np.log(T[idx2])
        X[:,1] = T[idx2]

        ga0,gb0,b2_2 = np.dot(np.dot(np.linalg.inv(np.dot(X.transpose(),X)),\
                          X.transpose()),S[idx2]/R)
    elif model == '7+constant+constant':
        X = np.ones([Ndata2,2])
        X[:,0] = np.log(T[idx2])

        ga0,b2_2 = np.dot(np.dot(np.linalg.inv(np.dot(X.transpose(),X)),\
                          X.transpose()),S[idx2]/R)

        X = np.ones([Ndata3,2])
        X[:,0] = np.log(T[idx3])

        ga0_2,b2_3 = np.dot(np.dot(np.linalg.inv(np.dot(X.transpose(),X)),\
                          X.transpose()),S[idx3]/R)

    # Init lmfit, left side:
    prms = lmfit.Parameters()
    prms.add('a1_1',value = ga1_1)
    prms.add('a2_1',value = ga2_1)
    prms.add('a3_1',value = ga3_1)
    prms.add('a4_1',value = ga4_1)
    prms.add('a5_1',value = ga5_1)
    prms.add('a6_1',value = ga6_1)

    # Constrain for the fit to match Cpr/R:
    #prms.add('a7',value = ga7)
    prms.add('a7_1', expr='('+str(Cpr/R)+'-(a1_1/('+str(Tr)+'**2)+a2_1/('+str(Tr)+') + a3_1'+\
                        '+a4_1*'+str(Tr)+'+a5_1*'+str(Tr)+'**2 + a6_1*'+str(Tr)+'**3))/('+str(Tr)+'**4)')

    # Constrain for the fit to match Hr/R: 
    prms.add('b1_1', expr=str(Hr/R)+'+(a1_1/'+str(Tr)+')-a2_1*log('+str(Tr)+')'+\
                        '-a3_1*'+str(Tr)+'-(a4_1/2)*('+str(Tr)+'**2)'+\
                        '-(a5_1/3)*('+str(Tr)+'**3)'+'-(a6_1/4)*('+str(Tr)+'**4)'+\
                        '-(a7_1/5)*('+str(Tr)+'**5)')

    # Constrain for the fit to match Sr/R:
    #prms.add('b2',value = Sr/R - (-ga1/(2.*Tr**2)-(ga2/Tr)+ga3*np.log(Tr)+ga4*Tr+\
    #                             (ga5/2)*Tr**2 + (ga6/3.)*Tr**3 + (ga7/4.)*Tr**4))
    prms.add('b2_1', expr=str(Sr/R)+'+a1_1/(2.*'+str(Tr)+'**2)+a2_1/('+str(Tr)+')'+\
                        '-a3_1*log('+str(Tr)+')-(a4_1*'+str(Tr)+')'+\
                        '-(a5_1/2)*('+str(Tr)+'**2)'+'-(a6_1/3)*('+str(Tr)+'**3)'+\
                        '-(a7_1/4)*('+str(Tr)+'**4)')

    # Now for the fit AFTER the phase transition:
    prms.add('a0',value = ga0)
    prms.add('b2_2',value = b2_2)

    if model == '7+constant':
        # Constrain fit to have zero difference in Gibbs free energy at the transition point:
        prms.add('b1_2',expr='-a1_1/(2*'+str(T_transition)+') + a2_1*(log('+str(T_transition)+')+1) + '+\
                         'a3_1*('+str(T_transition)+'-log('+str(T_transition)+')*'+str(T_transition)+')'+\
                         '-(a4_1/2.)*'+str(T_transition)+'**2-(a5_1/6.)*'+str(T_transition)+'**3 - '+\
                         '(a6_1/12.)*'+str(T_transition)+'**4 - (a7_1/20.)*'+str(T_transition)+'**5 + b1_1 '+\
                         '- b2_1*'+str(T_transition)+'-a0*('+str(T_transition)+'-'+str(T_transition)+'*log('+\
                         str(T_transition)+'))+'+str(T_transition)+'*b2_2')
    elif model == '7+linear':
        prms.add('b0',value = gb0)
        # Constrain fit to have zero difference in Gibbs free energy at the transition point:
        prms.add('b1_2',expr='-a1_1/(2*'+str(T_transition)+') + a2_1*(log('+str(T_transition)+')+1) + '+\
                         'a3_1*('+str(T_transition)+'-log('+str(T_transition)+')*'+str(T_transition)+')'+\
                         '-(a4_1/2.)*'+str(T_transition)+'**2-(a5_1/6.)*'+str(T_transition)+'**3 - '+\
                         '(a6_1/12.)*'+str(T_transition)+'**4 - (a7_1/20.)*'+str(T_transition)+'**5 + b1_1 '+\
                         '- b2_1*'+str(T_transition)+'-a0*('+str(T_transition)+'-'+str(T_transition)+'*log('+\
                         str(T_transition)+'))+'+str(T_transition)+'*b2_2 + b0*('+str(T_transition)+'**2/2.)')

    elif model == '7+constant+constant':
        # Constrain fit to have zero difference in Gibbs free energy at the first transition point:
        prms.add('b1_2',expr='-a1_1/(2*'+str(T_transition[0])+') + a2_1*(log('+str(T_transition[0])+')+1) + '+\
                         'a3_1*('+str(T_transition[0])+'-log('+str(T_transition[0])+')*'+str(T_transition[0])+')'+\
                         '-(a4_1/2.)*'+str(T_transition[0])+'**2-(a5_1/6.)*'+str(T_transition[0])+'**3 - '+\
                         '(a6_1/12.)*'+str(T_transition[0])+'**4 - (a7_1/20.)*'+str(T_transition[0])+'**5 + b1_1 '+\
                         '- b2_1*'+str(T_transition[0])+'-a0*('+str(T_transition[0])+'-'+str(T_transition[0])+'*log('+\
                         str(T_transition[0])+'))+'+str(T_transition[0])+'*b2_2')
 
        prms.add('a0_2',value = ga0_2)
        prms.add('b2_3',value = b2_3)

        # Constrain fit to have zero difference in Gibbs free energy at the second transition point:
        prms.add('b1_3',expr='a0*('+str(T_transition[1])+'-'+str(T_transition[1])+'*log('+str(T_transition[1])+\
                             ')) + b1_2 - '+str(T_transition[1])+'*b2_2 - a0_2*('+str(T_transition[1])+\
                             '-'+str(T_transition[1])+'*log('+str(T_transition[1])+')) + '+str(T_transition[1])+'*b2_3')
    # Prepare data:
    if len(model.split('+')) == 2:
        Y = np.zeros([Ndata1+Ndata2,3])
    elif len(model.split('+')) == 3:
        Y = np.zeros([Ndata1+Ndata2+Ndata3,3])
    Y[:,0] = Cp/R
    Y[:,1] = H/R
    Y[:,2] = S/R

    # Perform constrained fit:
    result = lmfit.minimize(residuals, prms, args=(T,Y))

    # If fit is successful, return fitted values. If not, cry:
    if result.success:
        if model == '7+constant':
            return result, result.params['a1_1'].value, result.params['a2_1'].value, \
               result.params['a3_1'].value, result.params['a4_1'].value, \
               result.params['a5_1'].value, result.params['a6_1'].value,\
               result.params['a7_1'].value, result.params['b1_1'].value,\
               result.params['b2_1'].value, result.params['a0'].value,\
               result.params['b1_2'].value, result.params['b2_2'].value
        elif model == '7+linear':
            return result, result.params['a1_1'].value, result.params['a2_1'].value, \
               result.params['a3_1'].value, result.params['a4_1'].value, \
               result.params['a5_1'].value, result.params['a6_1'].value,\
               result.params['a7_1'].value, result.params['b1_1'].value,\
               result.params['b2_1'].value, result.params['a0'].value,\
               result.params['b0'].value, result.params['b1_2'].value, result.params['b2_2'].value
        elif model == '7+constant+constant':
            return result, result.params['a1_1'].value, result.params['a2_1'].value, \
               result.params['a3_1'].value, result.params['a4_1'].value, \
               result.params['a5_1'].value, result.params['a6_1'].value,\
               result.params['a7_1'].value, result.params['b1_1'].value,\
               result.params['b2_1'].value, result.params['a0'].value,\
               result.params['b1_2'].value, result.params['b2_2'].value,\
               result.params['a0_2'].value, result.params['b1_3'].value,\
               result.params['b2_3'].value
    else:
        print 'Fit failed. Too bad :(.'
        return result,0,0,0,0,0,0,0,0,0
