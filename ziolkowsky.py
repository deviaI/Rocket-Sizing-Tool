# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 10:25:30 2022

@author: Devial
"""

import math
import numpy as np
import os.path
def f_reverse(n, m, m_pl, isp, **kwargs):
    """
    WORK IN PROGRESS NOT VALIDATED

    Method to return achievable delta v for a given Rocket Configuartion

    Inputs:
    n: Number of Stages
    m: Either Launch Mass of Rocket /wo Payload OR List of Stage Masses /wo Payload (such that the sum of the Stage Masses is the Launch Mass /wo Paylaod)
       If only Launch Mass is given, it is assumed that each stages mass is equal to half the previous Stages Mass
    m_pl: Mass of Payload
    isp: Engine Isp or List of Engine Isps
         If Only one Isp value is provided, it is assumed that Isp is identical for all stages
    
    Optional Input:
    m_f: List of propellant mass per stage
         If Not given, a strucure factor of 12% is assumed for each stage
    
    Returns:
        Achievable Delta V
    """
    _isp = []
    _m_s = []
    _m_f = []
    try:
        for i in range(0,n):
            _isp.append(isp[i])
    except TypeError:
        for i in range(0,n):
            _isp.append(isp)
    try:
        for i in range(0,n):
            _m_s.append(m[i])
    except TypeError:
        fac = 0
        for i in range(0,n):
            fac = fac + math.pow(0.5, i)
        _m_s.append((m)/fac)
        for i in range(1,n):
            _m_s.append(_m_s[i-1]*0.5)
    if "m_f" in kwargs:
        if type(kwargs["m_f"]) != list:
            raise TypeError("m_f must be list")
        for i in range(0,n):
            _m_f = kwargs["m_f"][i]
    else:
        for i in range(0,n):
            _m_f.append(_m_s[i] * 0.88)
    _m_s.append(m_pl)

    delta_v = 0
    for i in range(0,n):
        delta_v = delta_v + _isp[i]*9.81*math.log(sum(_m_s, i)/(sum(_m_s, i) - _m_f[i]), math.e)
    return delta_v

def f( mu, isp, m_pl, delv, limit):
    """
    Method for returning the launch mass of an SSTO based on:
    Engine Isp,  Structure Factor, Payload Mass and target Delta V
    
    Inputs:
        mu: Structure Factor 
        isp: Engine Isp in Seconds
        m_pl: Payload Mass in kg
        delv: Target Delta V
        limit: Upper mass limit for divergence cut off

    Returns:
        Launch Mass of Rocket in tons (if converged)

        1e99 (if diverged)
    """
    m_f = 250
    m_0 = 250
    delv = 9000
    m_01 = 250
    v_star = isp*9.81
    while True:
        m_01 = m_f * math.exp(delv/v_star)
        m_f = m_pl + mu * m_01
        if abs((m_01-m_0)/m_01) < 0.001:
            #print("converged for mu=", mu, "and isp =", isp)
            break
        if m_01 > limit:
            m_01 = 1e102
            #print("diverged for mu=", mu, "and isp =", isp)
            break
        m_0 = m_01
    return m_01/1000

def f_twoStage(mu, isp_1, isp_2, m_pl, delv, limit, **kwargs):
    """
    Method for returning the launch mass of an SSTO based on:
    Engine Isp,  Structure Factor, Payload Mass and target Delta V

    Inputs:
        mu: Structure Factor 
        isp_1: First stage engine Isp in Seconds
        isp_2: Second stage engine Isp in Seconds
        m_pl: Payload Mass in kg
        delv: Target Delta V
        limit: Upper mass limit for divergence cut off

    Optional Inputs:
        debug: debug Flag // Default: Off
        size_fac: Relative Mass of second stage compared to first stage // Default: 0.5

    Returns:
        Launch Mass of Rocket in tons(if converged)

        1e99 (if diverged)
    """
    m_0 = 2*m_pl
    v_star_1 = isp_1*9.81
    v_star_2 = isp_2*9.81
    if "size_fac" in kwargs:
        size_fac = kwargs["size_fac"]
    else:
        size_fac = 1/2
    #for i in range(0,20):
    while True:
        m_s1 = (m_0 - m_pl)/(1+size_fac)
        m_s2 = size_fac * m_s1
        m_01 = m_0
        m_02 = m_s2 + m_pl
        m_f1 = m_01 * mu + m_02
        m_f2 = m_02 * mu + m_pl
        sanity_check = ((m_s1 + m_s2 >= (m_0 - m_pl)*0.995) or ((m_s1 + m_s2 <= (m_0 - m_pl)*1.005)) ) and (m_f2 > m_pl)
        if not(sanity_check):
            print("sanity check failed")
            return -1
        delv_1 = v_star_1 * math.log(m_01/m_f1, math.e)
        delv_2 = v_star_2 * math.log(m_02/m_f2, math.e)
        delv_tot = delv_2 + delv_1
        m_0_ = (1+(delv-delv_tot)/delv) * m_0
        if "debug" in kwargs:
            print("Total delta v:" +str(delv_tot))
            print("Convergence Factor Launch masss " +str((1+(delv-delv_tot)/delv)))
            print("m_0_:" +str(m_0_))
            print("m_0:" +str(m_0))
            print("m_s1:" +str(m_s1))
            print("m_s2:" +str(m_s2))
            print("m_01:" +str(m_01))
            print("m_02:" +str(m_02))
            print("m_f1:" +str(m_f1))
            print("m_f2:" +str(m_f2))
            print("//")
        #print(m_0_)
        #print("//")
        if abs((m_0_-m_0)/m_0) < 0.0001:
            #print("converged for mu=", mu, "and isp =", isp)
            break
        if m_0_ > limit:
            m_0_ = 1e102
            #print("diverged for mu=", mu, "and isp =", isp)
            break
        m_0 = m_0_
    m_0 = m_0_
    return m_0/1000

def f_threeStage(mu, isp_1, isp_2, isp_3, m_pl, delv, limit, **kwargs):
    """
    Method for returning the launch mass of an SSTO based on:
    Engine Isp,  Structure Factor, Payload Mass and target Delta V
    Inputs:
        mu: Structure Factor 
        isp_1: First stage engine Isp in Seconds
        isp_2: Second stage engine Isp in Seconds
        isp_3: Third stage engine Isp in Seconds
        m_pl: Payload Mass in kg
        delv: Target Delta V
        limit: Upper mass limit for divergence cut off
    Optional Inputs:
        debug: debug Flag
    Returns:
        Launch Mass of Rocket in tons(if converged)
        1e99 (if diverged)
    """
    m_0 = 2*m_pl
    v_star_1 = isp_1*9.81
    v_star_2 = isp_2*9.81
    v_star_3 = isp_3*9.81
    while True:
        m_s1 = 4/7 * (m_0 - m_pl)
        m_s2 = 1/2 * m_s1
        m_s3 = 1/2 * m_s2
        m_01 = m_0
        m_02 = m_s2 + m_s3 + m_pl
        m_03 = m_s3 + m_pl
        m_f1 = m_0 * mu + m_02
        m_f2 = m_02 * mu + m_03
        m_f3 = m_03 * mu + m_pl
        delv_1 = v_star_1 * math.log(m_01/m_f1, math.e)
        delv_2 = v_star_2 * math.log(m_02/m_f2, math.e)
        delv_3 = v_star_3 * math.log(m_03/m_f3, math.e)
        delv_tot = delv_2 + delv_1 + delv_3
        sanity_check = ((m_s1 + m_s2 +m_s3 >= (m_0 - m_pl)*0.995) or ((m_s1 + m_s2 +m_s3 <= (m_0 - m_pl)*1.005)) ) and (m_f3 > m_pl)
        if not(sanity_check):
            print("sanity check failed")
            return -1
        m_0_ = (1+(delv-delv_tot)/delv) * m_0
        if "debug" in kwargs:
            print("Total delta v:" +str(delv_tot))
            print("Convergence Factor Launch masss " +str((1+(delv-delv_tot)/delv)))
            print("m_0_:" +str(m_0_))
            print("m_0:" +str(m_0))
            print("m_s1:" +str(m_s1))
            print("m_s2:" +str(m_s2))
            print("m_s3:" +str(m_s3))
            print("m_01:" +str(m_01))
            print("m_02:" +str(m_02))
            print("m_03:" +str(m_03))
            print("m_f1:" +str(m_f1))
            print("m_f2:" +str(m_f2))
            print("m_f3:" +str(m_f3))
            print("//")
        #print(m_0_)
        #print("//")
        if abs((m_0_-m_0)/m_0) < 0.0001:
            #print("converged for mu=", mu, "and isp =", isp)
            break
        if m_0_ > limit:
            m_0_ = 1e102
            #print("diverged for mu=", mu, "and isp =", isp)
            break
        m_0 = m_0_
    m_0 = m_0_
    return m_0/1000

def optimise_2Stage(m_0, isp1, isp2, mu, m_pl, **kwargs):
    """
    Method for finding the optimal mass ratio between first and second stage for 
    a given Launch Mass and Stage Isps
   Inputs:
        mu: Structure Factor 
        isp_1: First stage engine Isp in Seconds
        isp_2: Second stage engine Isp in Seconds
        m_pl: Payload Mass in kg
        delv: Target Delta V
        limit: Upper mass limit for divergence cut off
    Optional Inputs:
        debug: debug Flag // Default: Off
    Returns:
        Dictionary with entries:
            "Achieved Delta V": Delta V at optimal relative stage mass
            "Optimal Stage Sizing Factor": Optimal mass of second stage as a factor of first stage mass
    """
    delv_tot_max = 0
    size_fac_max = 0
    v_star_1 = isp1*9.81
    v_star_2 = isp2*9.81
    size_fac = np.linspace(0.01, 0.99, 99)
    for i in range(0,99):
        m_s1 = (m_0 - m_pl)/(1+size_fac[i])
        m_s2 = size_fac[i] * m_s1
        m_01 = m_0
        m_02 = m_s2 + m_pl
        m_f1 = m_01 * mu + m_02
        m_f2 = m_02 * mu + m_pl
        sanity_check = ((m_s1 + m_s2 >= (m_0 - m_pl)*0.995) or ((m_s1 + m_s2 <= (m_0 - m_pl)*1.005)) ) and (m_f2 > m_pl)
        if not(sanity_check):
            print("sanity check failed")
            return -1
        delv_1 = v_star_1 * math.log(m_01/m_f1, math.e)
        delv_2 = v_star_2 * math.log(m_02/m_f2, math.e)
        delv_tot = delv_2 + delv_1
        if "debug" in kwargs:
            print("Total delta v:" +str(delv_tot))
            print("m_0:" +str(m_0))
            print("m_s1:" +str(m_s1))
            print("m_s2:" +str(m_s2))
            print("m_01:" +str(m_01))
            print("m_02:" +str(m_02))
            print("m_f1:" +str(m_f1))
            print("m_f2:" +str(m_f2))
            print("//")
        if delv_tot > delv_tot_max:
            delv_tot_max = delv_tot
            size_fac_max = size_fac[i]
    results = {}
    results["Achieved Delta V"] = delv_tot_max
    results["Optimal Stage Sizing Factor"] = size_fac_max
    return results

            
def TwoDAlt(Stages, Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit, **kwargs):
    """
    Parent Method for generating a 2D Plotable csv Data set of Launch Mass against Isp or Launch Mass against Structure Factor

    Inputs:
        Stages: Number of Stages (1,2 or 3)
        Mode: Flag determining Mode ("Isp" for constant Isp, "Mu" for constant Mu)
        FixVal: Value of fixed Parameter (Isp in "Isp" Mode or Mu in "Mu" Mode)
        size_X_Axis: Number of Calculation Points along X-Axis
        X_Axis_LL:  Lower Limit of X-Axis
        X_Axis_UL: Upper Limit of X-Axis
    Optional Inputs:
        isp_2: Isp of second stage engine (only usable for Mode="Isp" and Stages >=2) // Default: FixVal

        isp_3: Isp of third stage engine (only usable for Mode="Isp" and Stages =3) // Default: isp_2
    Returns:
        None 
    """
    if Stages == 1:
        TwoDAlt_SSTO(Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit)
    elif Stages == 2:
        if "isp_2" in kwargs:
            isp2 = kwargs["isp_2"]
        else:
            isp2 = FixVal 
        TwoDAlt_2Stage(Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit,isp_2 = isp2)
    else:
        if "isp_2" in kwargs:
            isp2 = kwargs["isp_2"]
        else:
            isp2 = FixVal
        if "isp_3" in kwargs:
            isp3 = kwargs["isp_3"]
        else:
            isp3 = isp2
        TwoDAlt_3Stage(Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit, isp_2 = isp2, isp_3 = isp3)
               
def TwoDAlt_SSTO(Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit):
    """
    Method for generating a 2D Plotable csv Data set of Launch Mass against Isp or Launch Mass against Structure Factor
    for an SSTO
    Inputs:
        Mode: Flag determining Mode ("Isp" for constant Isp, "Mu" for constant Mu)
        FixVal: Value of fixed Parameter (Isp in "Isp" Mode or Mu in "Mu" Mode)
        size_X_Axis: Number of Calculation Points along X-Axis
        X_Axis_LL:  Lower Limit of X-Axis
        X_Axis_UL: Upper Limit of X-Axis

    Returns:
        None 
    """

    if Mode == "Isp":
        isp = FixVal
        mu = np.linspace(x_Axis_UL, X_Axis_LL, size_X_Axis)
        Y = np.zeros((size_X_Axis, 2))
        fName = os.path.join(os.path.dirname(__file__), "data",  "SSTO_ISP=" + str(round(isp, 2))+ ".csv")
        for k in range(0, size_X_Axis):
            Y[k, 0] = f(mu[k], isp, m_pl, delv, limit)
            Y[k, 1] = mu[k]
        np.savetxt(fName, Y, delimiter=",")
        with open(fName,'a') as fd:
            fd.write("ISP = " + str(isp) + " Configuration: SSTO")
    else:
        isp = np.linspace(X_Axis_LL, x_Axis_UL, size_X_Axis)
        mu = FixVal
        Y = np.zeros((size_X_Axis,2))
        fName = os.path.join(os.path.dirname(__file__), "data", "SSTO_Mu=" + str(round(mu, 5))+ ".csv")
        for k in range(0, size_X_Axis):
            Y[k] = f(mu, isp[k], m_pl, delv, limit)
            Y[k, 1] = isp[k]
        np.savetxt(fName, Y, delimiter=",")
        with open(fName,'a') as fd:
            fd.write("Mu = " + str(mu) + " Configuration: SSTO")

        
def TwoDAlt_2Stage(Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit, **kwargs):
    """
    Method for generating a 2D Plotable csv Data set of Launch Mass against Isp or Launch Mass against Structure Factor
    for a 2 Stage Rocket
    Inputs:
        Mode: Flag determining Mode ("Isp" for constant Isp, "Mu" for constant Mu)
        FixVal: Value of fixed Parameter (Isp in "Isp" Mode or Mu in "Mu" Mode)
        size_X_Axis: Number of Calculation Points along X-Axis
        X_Axis_LL:  Lower Limit of X-Axis
        X_Axis_UL: Upper Limit of X-Axis
    Optional Inputs:
        isp_2: Isp of second stage engine (Only usable if Mode is set to "Isp") // Default: FixVal
    Returns:
        None 
    """

    if Mode == "Isp":
        isp = FixVal
        if "isp_2" in kwargs:
            isp2 = kwargs["isp_2"]
        else:
            isp2 = isp
        mu = np.linspace(x_Axis_UL, X_Axis_LL, size_X_Axis)
        Y = np.zeros((size_X_Axis, 2))
        fName = os.path.join(os.path.dirname(__file__), "data", "2Stage_ISP1=" + str(round(isp, 2)) + ";ISP2=" + str(round(isp2, 2))+ ".csv")
        for k in range(0, size_X_Axis):
            Y[k, 0] = f_twoStage(mu[k], isp, isp2, m_pl, delv, limit)
            Y[k, 1] = mu[k]
        np.savetxt(fName, Y, delimiter=",")
        with open(fName,'a') as fd:
            fd.write("Configuration: 2-Stage, mass Stage 1 = 2 x mass Stage 2" + " Isp Stage 1:" + str(isp) + "s Isp Stage 2: " + str(isp2) + "s")
    else:
        isp = np.linspace(X_Axis_LL, x_Axis_UL, size_X_Axis)
        mu = FixVal
        Y = np.zeros((size_X_Axis,2))
        fName = os.path.join(os.path.dirname(__file__), "data", "2Stage_Mu=" + str(round(mu, 5))+ ".csv")
        for k in range(0, size_X_Axis):
            Y[k] = f_twoStage(mu, isp[k], 1e6, isp[k])
            Y[k, 1] = isp[k]
        np.savetxt(fName, Y, delimiter=",")
        with open(fName,'a') as fd:
            fd.write(" Configuration: 2-Stage, mass Stage 1 = 2 x mass Stage 2, Mu = " + str(mu) )

def TwoDAlt_3Stage(Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit, **kwargs):
    """
    Method for generating a 2D Plotable csv Data set of Launch Mass against Isp or Launch Mass against Structure Factor
    for a 3 Stage Rocket
    Inputs:
        Mode: Flag determining Mode ("Isp" for constant Isp, "Mu" for constant Mu)
        FixVal: Value of fixed Parameter (Isp in "Isp" Mode or Mu in "Mu" Mode)
        size_X_Axis: Number of Calculation Points along X-Axis
        X_Axis_LL:  Lower Limit of X-Axis
        X_Axis_UL: Upper Limit of X-Axis
    Optional Inputs:
        isp_2: Isp of second stage engine (Only usable if Mode is set to "Isp") // Default: FixVal
        isp_3: Isp of third stage engine (Only usable if Mode is set to "Isp") // Default: isp_2
    Returns:
        None 
    """
    if Mode == "Isp":
        isp = FixVal
        if "isp_2" in kwargs:
            isp2 = kwargs["isp_2"]
        else:
            isp2 = isp
        if "isp_3" in kwargs:
            isp3 = kwargs["isp_3"]
        else:
            isp3 = isp2
        mu = np.linspace(x_Axis_UL, X_Axis_LL, size_X_Axis)
        Y = np.zeros((size_X_Axis, 2))
        fName = os.path.join(os.path.dirname(__file__), "data", "3 Stage_ISP1=" + str(round(isp, 2)) + ";ISP2=" + str(round(isp2, 2)) + ";ISP3=" + str(round(isp3, 2)) + ".csv")
        for k in range(0, size_X_Axis):
            Y[k, 0] = f_threeStage(mu[k], isp, isp2, isp3, m_pl, delv, limit, )
            Y[k, 1] = mu[k]
        np.savetxt(fName, Y, delimiter=",")
        with open(fName,'a') as fd:
            fd.write("Configuration: 3-Stage, mass Stage 1 = 2 x mass Stage 2 = 2 x mass Stage 3" + " Isp Stage 1:" + str(isp) + "s Isp Stage 2: " + str(isp2) + "s Isp Stage 3: " + str(isp3) + "s")
    else:
        isp = np.linspace(X_Axis_LL, x_Axis_UL, size_X_Axis)
        mu = FixVal
        Y = np.zeros((size_X_Axis,2))
        fName = os.path.join(os.path.dirname(__file__), "data", "3Stage_Mu=" + str(round(mu, 5))+ ".csv")
        for k in range(0, size_X_Axis):
            Y[k] = f_twoStage(mu, isp[k], 1e6, isp[k], isp[k])
            Y[k, 1] = isp[k]
        np.savetxt(fName, Y, delimiter=",")
        with open(fName,'a') as fd:
            fd.write(" Configuration: 3-Stage, mass Stage 1 = 2 x mass Stage 2 = 2 x mass Stage 3, Mu = " + str(mu) )

def Optimised2Stage(Isp1, Isp2, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit):
    """
    Method for generating a 2D Plotable csv Data set of Launch Mass against Isp or Launch Mass against Structure Factor
    for a 2 Stage Rocket with optimised relative mass between stages
    Inputs:
       Isp1: Isp of first stage engine
       Isp2: Isp of second stage engine
       size_x_Axis: Number of calculation steps for mu range
       X_Axis_LL: Lower limit of Mu
       x_Axis_UL: Upper limit of Mu
    Returns:
        None 
    """

    mu = np.linspace(x_Axis_UL, X_Axis_LL, size_X_Axis)
    Y = np.zeros((size_X_Axis, 3))
    fName = os.path.join(os.path.dirname(__file__), "data", "2-Stage Opt.;ISP1=" + str(Isp1) + "s;ISP2=" +str(Isp2) +"s"+ ".csv") 
    for k in range(0, size_X_Axis):
        m_0 = 1e-3
        m_0_ = f_twoStage(mu[k], Isp1, Isp2, m_pl, delv, limit)
        while abs((m_0_-m_0)/m_0) > 0.0001:
            m_0 = m_0_
            results = optimise_2Stage(m_0*1000, Isp1, Isp2, mu[k], m_pl)
            m_0_ = f_twoStage(mu[k], Isp1, Isp2, m_pl, delv, limit, size_fac = results["Optimal Stage Sizing Factor"])
        Y[k,0] = m_0_
        Y[k,1] = mu[k]
        Y[k,2] = results["Optimal Stage Sizing Factor"]
    np.savetxt(fName, Y, delimiter=",")
    with open(fName,'a') as fd:
        fd.write(" Configuration: 2-Stage optimised , Relative Stage Sizing = third column, Isp1=" + str(Isp1) + ";Isp2=" + str(Isp2))
    