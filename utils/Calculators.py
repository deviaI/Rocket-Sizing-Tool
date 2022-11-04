# -*- coding: utf-8 -*-
"""
Created on  2022/10/26 09:23:69

@author: Devial
"""

import math
import numpy as np


class Calculator(object):


    def Tsiolkowsky(self, isp, m0, mf):
        """
        Tsiolkowsky Equation

        Inputs:
            Isp: Isp
            m0: Mass at start of burn
            mf: Mass at end of burn

        Returns:
            delta v
        """
        return isp*9.81*math.log(m0/mf, math.e)

    def MassSplit(self, n, m_0, m_pl, mu, size_fac):
        """
        Computes list of fuel masses and total masses for each individual stage

        Inputs:
            n: number of stages
            m_0: total launch mass
            m_pl: payload mass
            mu: Structural Factor
            mu: relative stage to stage mass
        
        Returns:
            list of stage masses (including payload as final stage), such that the sum of stage masses is the total launch mass
            list of fuel masses (one element less than stage masses)
        """
        m_s = []
        m_f = []
        fac = 0
        for k in range(0,n):
            fac = fac + math.pow(size_fac, k)
        m_s.append((m_0-m_pl)/fac)
        for k in range(1,n):
            m_s.append(m_s[k-1]*size_fac)
        m_s.append(m_pl)
        for k in range(0,n):
            m_f.append(m_s[k] * (1-mu))
        sanity_check = sum(m_s) < m_0*1.01 and sum(m_s) > m_0*0.99
        if not(sanity_check):
            print("sanity check failed")
            return -1
        return m_s, m_f

    def calcDelV(self, n, m, m_pl, isp, **kwargs):
        """
        Method to return achievable delta v for a given Rocket Configuartion

        Inputs:
        n: Number of Stages
        m: Total Launch Mass of Rocket 
           OR
           List of Stage Masses /wo Payload (such that the sum of the Stage Masses is the Launch Mass /wo Paylaod)
           ///If only launch mass is given, an optimal stage to stage mass ratio is determined
        m_pl: Mass of Payload
        isp: Engine Isp or List of Engine Isps
            If Only one Isp value is provided, it is assumed that Isp is identical for all stages
        
        Optional Input:
        m_f: List of propellant mass per stage
            //If Not given, a structure factor of 12% is assumed. 
            //Can only be specified if m is a list of stage masses
        
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
        if "m_f" in kwargs:
            if type(kwargs["m_f"]) != list:
                raise TypeError("m_f must be list")
            if type(m) != list:
                raise TypeError("m must be list if m_f is given")
            for i in range(0,n):
                _m_s.append(m[i])
            for i in range(0,n):
                 _m_f.append(kwargs["m_f"][i])
        else:
            if type(m) != list:
                _m_s, _m_f = self.MassSplit(n, m, m_pl, 0.12, 0.5)
            else:
                for i in range(0,n):
                    _m_s.append(m[i])
                    _m_f = _m_s[i] * 0.88
        delta_v = 0
        for i in range(0,n):
            delta_v += self.Tsiolkowsky(_isp[i], sum(_m_s[i:n]), sum(_m_s[i:n]) - _m_f[i])
        return delta_v


    def calcPoint(self, n, isp, m_pl, mu = 0.12, delv = 9000, limit = 1e6, **kwargs):
        """
        Method for returning the theoretical launch mass of a Rocket based on:
        Engine Isp,  Structure Factor, Payload Mass and target Delta V

        Inputs:
            n:   Number of Stages
            Isp: Isp of all Engines
                 OR
                 List of Isps for each Engine
            m_pl: Payload Mass in kg
            mu: Structure Factor // Default = 0.12
            delv: Target Delta V // Default = 90000 m/s
            limit: Upper mass limit for divergence cut off // Default = 1e6

        Optional Inputs:
            size_fac: relative stage to stage mass
                      If not given, it is calculated to be optimal
        Returns:
            Launch Mass of Rocket in kg(if converged)
            OR
            1e99 (if diverged)

            Relative Stage Mass Factor (will be identical to input, if given as input)
        """

        m_0 = 2*m_pl
        if "size_fac" in kwargs:
            size_fac = kwargs["size_fac"]
        while True:   
            delv_tot = 0
            if not("size_fac" in kwargs):
                size_fac = self.optimiseMassRatio(n, m_0, isp, m_pl, mu)["Optimal Stage Sizing Factor"]
            m_s, m_f = self.MassSplit(n, m_0, m_pl, mu, size_fac)
            for k in range(0,n):
                if type(isp) == list:
                    delv_tot += self.Tsiolkowsky(isp[k], sum(m_s[k:n+1]), sum(m_s[k:n+1]) - m_f[k])
                else:
                    delv_tot += self.Tsiolkowsky(isp, sum(m_s[k:n+1]), sum(m_s[k:n+1]) - m_f[k])
            m_0_ = (1+(delv-delv_tot)/delv) * m_0
            if abs((m_0_-m_0)/m_0) < 0.0001:
                break
            if m_0_ > limit:
                m_0_ = 1e99
                break
            m_0 = m_0_
        m_0 = m_0_
        return m_0, size_fac

    
    def calcRange(self, n, RangeVariable,  isp, m_pl, mu, delv = 9000, limit = 1e6, numSteps=100, **kwargs):
        """
        Method for returning the theoretical launch mass of a Rocket based on:
        Engine Isp,  Structure Factor, Payload Mass and target Delta V

        Inputs:
            n:   Number of Stages
            Isp: Isp of all Engines
                 OR
                 List of Isps for each Engine
            m_pl: Payload Mass in kg
            mu: Structure Factor // Default = 0.12
            delv: Target Delta V // Default = 90000 m/s
            limit: Upper mass limit for divergence cut off // Default = 1e6

        Optional Inputs:
            size_fac: Relative Mass of second stage compared to first stage 
                      If not given, it is calculated to be optimal

        Returns:
            np.array with column [0] containing mass and column [1] containing the corresponding Range Variable (either isp or mu)
        """
        Y = np.zeros((3, numSteps))
        if RangeVariable == "isp":
            if not(type(isp) == list):
                raise TypeError("RangeVariable must be list")
            isp = np.linspace(isp[0], isp[-1], numSteps)
            for i in range(0, numSteps):
                if "size_fac" in kwargs:
                    Y[0, i], Y[2, i] = self.calcPoint(n, isp[i], m_pl, mu, delv, limit, kwargs["size_fac"])
                else:
                    Y[0, i], Y[2, i] = self.calcPoint(n, isp[i], m_pl, mu, delv, limit)
                Y[1, i] = isp[i]
        elif RangeVariable == "mu":
            if not(type(mu) == list):
                raise TypeError("RangeVariable must be list")
            mu = np.linspace(mu[0], mu[-1], numSteps)
            for i in range(0, numSteps):
                if "size_fac" in kwargs:
                    Y[0, i], Y[2, i] = self.calcPoint(n, isp, m_pl, mu[i], delv, limit, kwargs["size_fac"])
                else:
                    Y[0, i], Y[2, i] = self.calcPoint(n, isp, m_pl, mu[i], delv, limit)
                Y[1, i] = mu[i]
        else:
            raise Exception("Unknown RangeVariable")
        return Y


    def optimiseMassRatio(self, n, m_0, isp, m_pl, mu = 0.12):
        """
        Method for finding the optimal mass ratio between stages for 
        a given Launch Mass and Stage Isps
    Inputs:
            n:  Number of Stages
            mu: Structure Factor//Default = 12%
            isp1: Isp of all Engines 
                  OR 
                  List of Isps for each Engine
            m_pl: Payload Mass in kg
    Returns:
            Dictionary with entries:
                "Achieved Delta V": Delta V at optimal relative stage mass
                "Optimal Stage Sizing Factor": Optimal mass of second/third... stage as a factor of first/second... stage mass
        """
        delv_tot_max = 0
        size_fac_max = 0
        for i in range(1,99):
            delv_tot = 0
            m_s, m_f = self.MassSplit(n, m_0, m_pl, mu, i*1e-2)
            sanity_check = sum(m_s) < m_0*1.01 and sum(m_s) > m_0*0.99
            if not(sanity_check):
                print("sanity check failed")
                print(m_s)
                print(i)
                return -1
            for k in range(0,n):
                if type(isp) == list:
                    delv_tot += self.Tsiolkowsky(isp[k], sum(m_s[k:n+1]), sum(m_s[k:n+1]) - m_f[k])
                else:
                    delv_tot += self.Tsiolkowsky(isp, sum(m_s[k:n+1]), sum(m_s[k:n+1]) - m_f[k])
            if delv_tot > delv_tot_max:
                delv_tot_max = delv_tot
                size_fac_max = i*1e-2
        results = {}
        results["Achieved Delta V"] = delv_tot_max
        results["Optimal Stage Sizing Factor"] = size_fac_max
        return results


    #//////////////////////////////////////
    #OBSOLETE
    #/////////////////////////////////////0
    def f(self, mu, isp, m_pl, delv, limit):
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

    def f_twoStage(self, mu, isp_1, isp_2, m_pl, delv, limit, **kwargs):
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
            delv_1 = self.Tsiolkowsky(isp_1, m_01, m_f1)
            delv_2 = self.Tsiolkowsky(isp_2, m_02, m_f2)
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


    def f_threeStage(self, mu, isp_1, isp_2, isp_3, m_pl, delv, limit, **kwargs):
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
            delv_1 = self.Tsiolkowsky(isp_1, m_01, m_f1)
            delv_2 = self.Tsiolkowsky(isp_2, m_02, m_f2)
            delv_3 = self.Tsiolkowsky(isp_3, m_03, m_f3)
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

    def TwoDAlt(self, Stages, Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit, **kwargs):
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
            return self.TwoDAlt_SSTO(Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit)
        elif Stages == 2:
            if "isp_2" in kwargs:
                isp2 = kwargs["isp_2"]
            else:
                isp2 = FixVal 
            return self.TwoDAlt_2Stage(Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit,isp_2 = isp2)
        else:
            if "isp_2" in kwargs:
                isp2 = kwargs["isp_2"]
            else:
                isp2 = FixVal
            if "isp_3" in kwargs:
                isp3 = kwargs["isp_3"]
            else:
                isp3 = isp2
            return self.TwoDAlt_3Stage(Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit, isp_2 = isp2, isp_3 = isp3)
                
    def TwoDAlt_SSTO(self, Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit):
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

            for k in range(0, size_X_Axis):
                Y[k, 0] = self.f(mu[k], isp, m_pl, delv, limit)
                Y[k, 1] = mu[k]
        else:
            isp = np.linspace(X_Axis_LL, x_Axis_UL, size_X_Axis)
            mu = FixVal
            Y = np.zeros((size_X_Axis,2))
    
            for k in range(0, size_X_Axis):
                Y[k] = self.f(mu, isp[k], m_pl, delv, limit)
                Y[k, 1] = isp[k]
        return Y


            
    def TwoDAlt_2Stage(self, Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit, **kwargs):
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
            for k in range(0, size_X_Axis):
                Y[k, 0] = self.f_twoStage(mu[k], isp, isp2, m_pl, delv, limit)
                Y[k, 1] = mu[k]
        else:
            isp = np.linspace(X_Axis_LL, x_Axis_UL, size_X_Axis)
            mu = FixVal
            Y = np.zeros((size_X_Axis,2))
            for k in range(0, size_X_Axis):
                Y[k] = self.f_twoStage(mu, isp[k], 1e6, isp[k])
                Y[k, 1] = isp[k]
        return Y

    def TwoDAlt_3Stage(self, Mode, FixVal, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit, **kwargs):
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
            for k in range(0, size_X_Axis):
                Y[k, 0] = self.f_threeStage(mu[k], isp, isp2, isp3, m_pl, delv, limit, )
                Y[k, 1] = mu[k]
        else:
            isp = np.linspace(X_Axis_LL, x_Axis_UL, size_X_Axis)
            mu = FixVal
            Y = np.zeros((size_X_Axis,2))
            for k in range(0, size_X_Axis):
                Y[k] = self.f_twoStage(mu, isp[k], 1e6, isp[k], isp[k])
                Y[k, 1] = isp[k]
        return Y

    def Optimised2Stage(self, Isp1, Isp2, size_X_Axis, X_Axis_LL, x_Axis_UL, m_pl, delv, limit):
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
        for k in range(0, size_X_Axis):
            m_0 = 1e-3
            m_0_ = self.f_twoStage(mu[k], Isp1, Isp2, m_pl, delv, limit)
            while abs((m_0_-m_0)/m_0) > 0.0001:
                m_0 = m_0_
                results = self.optimise_2Stage(m_0*1000, Isp1, Isp2, mu[k], m_pl)
                m_0_ = self.f_twoStage(mu[k], Isp1, Isp2, m_pl, delv, limit, size_fac = results["Optimal Stage Sizing Factor"])
            Y[k,0] = m_0_
            Y[k,1] = mu[k]
            Y[k,2] = results["Optimal Stage Sizing Factor"]
        return Y
        