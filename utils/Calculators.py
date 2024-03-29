# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 10:25:30 2022

@author: Devial
"""

import math
import numpy as np
import utils.ThermalSizing as TS
import utils.TM as TM
class Calculator(object):

    def __init__(self):
        """
        Calculator Module
        All Inputs/Outputs in SI units
        """

    def Tsiolkowsky(self, isp, m0, mf):
        """
        Simple Tsiolkowsky Calculator

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
            mu: Structural Factor or list of structural factors per stage
            size_fac: relative stage to stage mass
        
        Returns:
            list of stage masses (including payload as final stage), such that the sum of stage masses is the total launch mass
            list of fuel masses (one element less than stage masses)
        """
        m_s = []
        m_f = []
        try:
            if len(mu) == 1:
                mu = mu[0]
            elif len(mu) != n:
                raise ValueError("Mu must either be value, or equal in length to number of Stages")
        except:
            mu = mu
        fac = 0
        for k in range(0,n):
            fac = fac + math.pow(size_fac, k)
        m_s.append((m_0-m_pl)/fac)
        for k in range(1,n):
            m_s.append(m_s[k-1]*size_fac)
        m_s.append(m_pl)
        for k in range(0,n):
            try:
                m_f.append(m_s[k] * (1-mu[k]))
            except:
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
        m: Launch Mass of Rocket /wo Payload
           OR
           List of Stage Masses /wo Payload (such that the sum of the Stage Masses is the Launch Mass /wo Paylaod)
           ///If only launch mass is given, an optimal stage to stage mass ratio is determined
        m_pl: Mass of Payload
        isp: Engine Isp or List of Engine Isps
            If Only one Isp value is provided, it is assumed that Isp is identical for all stages
        
        Optional Input:
        m_f: List of propellant mass per stage
            //If Not given, a structure factor of 12% is assumed. 
            //Can only be specified if len(m) == len(m_f) == n (can be value or array like for n = 1, must be array like for n>1)
        
        Returns:
            Achievable Delta V
        """
        _isp = []
        _m_s = []
        _m_f = []
        _mu = []
        try:
            for i in range(0,n):
                _isp.append(isp[i])
        except TypeError:
            for i in range(0,n):
                _isp.append(isp)
        if "m_f" in kwargs:
            m_f = kwargs["m_f"]
            if n > 1:
                try:
                    if len(m) == len(m_f):
                        for i in range(0, n):
                            _m_s.append(m[i])
                            _m_f.append(m_f[i])
                            _mu.append(1 - m_f[i]/m[i])
                except TypeError:
                    raise ValueError("m and m_f must both be arraay like with length n")
            else:
                try:
                    _m_s.append(m[0])
                except:
                    _m_s.append(m)
                try:
                    _m_f.append(m_f[0])
                except:
                    _m_f.append(m_f)
                _mu.append(1 - _m_f[0]/_m_s[0])
        else:
            try:
                for i in range(0,n):
                    _m_s.append(m[i])
                    _m_f.append(_m_s[i] * 0.88)
                    _mu.append(0.12)
            except TypeError:
                _mu.append(0.12)
                return self.optimiseMassRatio(n,m + m_pl, isp, m_pl, _mu)["Achieved Delta V"]
        _m_s.append(m_pl)
        delta_v = 0
        for i in range(0,n):
            delta_v += self.Tsiolkowsky(_isp[i], sum(_m_s[i:n+1]), sum(_m_s[i:n+1]) - _m_f[i])
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
        max_dv = 0
        for i in range(0,n):
            try:
                max_dv += self.Tsiolkowsky(isp[i], 1, mu)
            except:
                max_dv += self.Tsiolkowsky(isp, 1, mu)
        if max_dv <= delv:
            raise Exception("Target Delta V unachievable")
        m_0 = 2*m_pl
        if "size_fac" in kwargs:
            size_fac = kwargs["size_fac"]
        while True:   
            delv_tot = 0
            if not("size_fac" in kwargs):
                size_fac = self.optimiseMassRatio(n, m_0, isp, m_pl, mu)["Optimal Stage Sizing Factor"]
            m_s, m_f = self.MassSplit(n, m_0, m_pl, mu, size_fac)
            for k in range(0,n):
                try:
                    delv_tot += self.Tsiolkowsky(isp[k], sum(m_s[k:n+1]), sum(m_s[k:n+1]) - m_f[k])
                except:
                    delv_tot += self.Tsiolkowsky(isp, sum(m_s[k:n+1]), sum(m_s[k:n+1]) - m_f[k])
            m_0_ = (1+(delv-delv_tot)/delv) * m_0
            if abs((m_0_-m_0)/m_0) < 0.0001:
                break
            if m_0_ > limit:
                m_0_ = 1e99
                break
            m_0 = m_0_
        m_0 = m_0_
        return m_0, size_fac, m_f

    def calcRange(self, n, RangeVariable,  isp, m_pl, mu, delv = 9000, limit = 1e6, numSteps=100, **kwargs):
        """
        Method for returning the theoretical launch mass of a Rocket based on:
        Engine Isp,  Structure Factor, Payload Mass and target Delta V

        Inputs:
            n:   Number of Stages
            RangeVariable: "Isp" or "Mu"
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
            np.array with column [0] containing mass, column [1] containing the corresponding Range Variable (either isp or mu) and column 3 containing the rel. stage sizing
        """
        Y = np.zeros((3, numSteps))
        if RangeVariable == "Isp":
            try:
                isp = np.linspace(isp[0], isp[-1], numSteps)
            except TypeError:
                raise TypeError("Range Variable must be array like")
            for i in range(0, numSteps):
                if "size_fac" in kwargs:
                    Y[0, i], Y[2, i] = self.calcPoint(n, isp[i], m_pl, mu, delv, limit, kwargs["size_fac"])[0:2]
                else:
                    Y[0, i], Y[2, i] = self.calcPoint(n, isp[i], m_pl, mu, delv, limit)[0:2]
                Y[1, i] = isp[i]
        elif RangeVariable == "Mu":
            try:
                mu = np.linspace(mu[0], mu[-1], numSteps)
            except:
                raise TypeError("RangeVariable must be list")
            for i in range(0, numSteps):
                if "size_fac" in kwargs:
                    Y[0, i], Y[2, i] = self.calcPoint(n, isp, m_pl, mu[i], delv, limit, kwargs["size_fac"])[0:2]
                else:
                    Y[0, i], Y[2, i] = self.calcPoint(n, isp, m_pl, mu[i], delv, limit)[0:2]
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
        try:
            if len(isp) == 1:
                isp = isp[0]
        except:
            pass
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
                try:
                    delv_tot += self.Tsiolkowsky(isp[k], sum(m_s[k:n+1]), sum(m_s[k:n+1]) - m_f[k])
                except:
                    delv_tot += self.Tsiolkowsky(isp, sum(m_s[k:n+1]), sum(m_s[k:n+1]) - m_f[k])
            if delv_tot > delv_tot_max:
                delv_tot_max = delv_tot
                size_fac_max = i*1e-2
        results = {}
        results["Achieved Delta V"] = delv_tot_max
        results["Optimal Stage Sizing Factor"] = size_fac_max
        return results

    def calcReqFuel(self, n, isp, m_pl, mu = 0.12, delv = 9000, limit = 1e6, **kwargs):
        """
        Method for returning the required fuel of a Rocket based on:
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
            mix_rat: Mass mixture ratio O/F 
            dens: density of propelant [dens Oxidzer, dens Fuel]
            fueltype: fuel type. If handed, mix_rat and dens are taken from typical values
                    Recognised Values: HydroLox, KeroLox, MethaLox
                    // if both fueltype and mix_rat are handed, the specified mixture ratio will be used
        Returns:
            Total Propellant mass
            Oxidzer Mass if "fueltype" or "mix_rat" is given. Else -1
            Fuel Mass if "fueltype" or "mix_rat" is given. Else -1
            Oxidzer Volume (m^3) if "fueltype" or "mix_rat" and "dens" is given. Else -1
            Fuel Volume (m^3) if "fueltype" or "mix_rat" and "dens" is given. Else -1
            

        """
        fuelSpecs = {}
        fuelSpecs["density"] = {}
        fuelSpecs["Mix Ratio"] = {}
        fuelSpecs["density"]["HydroLox"] = (1140 , 71)
        fuelSpecs["density"]["KeroLox"] = (1140 , 800)
        fuelSpecs["density"]["MethaLox"] = (1140 , 423)
        fuelSpecs["Mix Ratio"]["HydroLox"] = 6
        fuelSpecs["Mix Ratio"]["KeroLox"] = 2.3
        fuelSpecs["Mix Ratio"]["MethaLox"] = 3.5
        dels = []
        for val in kwargs:
            if kwargs[val] == -1:
                dels.append(val)
        for val in dels:
            del kwargs[val] 
        if "size_fac" in kwargs:
            m_p = self.calcPoint(n, isp, m_pl, mu, delv, limit, size_fac = kwargs["size_fac"])[-1]
        else:
            m_p = self.calcPoint(n, isp, m_pl, mu, delv, limit)[-1]
        m_p = sum(m_p)
        if "fueltype" in kwargs:
            try:
                dens = fuelSpecs["density"][kwargs["fueltype"]]
                mix_rat = fuelSpecs["Mix Ratio"][kwargs["fueltype"]]
            except:
                raise ValueError("Unrecognised Fuel Type. Options are: 'HydroLox', 'KeroLox' or 'MethaLox'")
        elif "mix_rat" not in kwargs:
            return m_p, -1, -1, -1, -1
        if "mix_rat" in kwargs:
            mix_rat = kwargs["mix_rat"]
        if "dens" in kwargs:
            dens = kwargs["dens"]
        m_o =  (m_p * mix_rat)/(1+mix_rat)
        m_f =  m_p/(1+mix_rat)
        if "dens" not in kwargs and "fueltype" not in kwargs:
            return m_p, m_o, m_f, -1, -1
        else:
            v_o =  (m_o / dens[0])*1000
            v_f = (m_f / dens[1])*1000
            return m_p, m_o, m_f, v_o, v_f

    def calcAscent_Ideal_DV(self, LEOalt):
        """
        Calculate ideal ascent delta v

        Inputs:
            LEOalt: Altitude of target LEO
        Returns:
            Ideal ascent delta v
        """
        r = 6378000 + LEOalt
        mu_earth = 3.986e14
        dv_ascent = 2*9.81*(6378000-6378000**2/(6378000+LEOalt))
        dv_orbit = mu_earth*(1/r)
        dv_total = np.sqrt(dv_ascent+dv_orbit)
        return dv_total

    def gen_ascent_path_preview(self, C1, C2,h_grav_turn = 1000,  h_cutoff = 200000):
        """
        Generate a preview of the ascent profile
        h(x) = (C2 * x)^(C1)
        Inputs:
            C1: Constant
            C2: Constant
            g
        Returns:
            h: Array
            x: Array
        """
        x = 1e-6 #Function not defined at x = 0, thus start "infinitessimally" close to 0
        dx = 10
        gamma = np.pi/2
        while abs(gamma - np.deg2rad(87)) > 1e-5:
            while gamma > np.deg2rad(87):
                gamma = np.arctan(C1*(C2**(C1))*x**(C1-1))
                x+=dx
            x -= 2*dx
            dx *= 0.5
            gamma = np.arctan(C1*(C2**(C1))*(x**(C1-1)))
        h_offset = h_grav_turn - (C2 * (x))**(C1) 
        x_step = 0
        dx = 1e-3
        x_array = []
        h_array = []
        x_array.append(0)
        x_array.append(0)
        h_array.append(0)
        h_array.append(h_grav_turn)
        h_step = h_grav_turn
        while h_step <= h_cutoff:
            if h_step >= C2*4 and h_step <= C2*4+1000:
                temp = 0
            x_step += dx
            h_step = (C2 * (x_step + x))**(C1) + h_offset
            h_array.append(h_step)
            x_array.append(x_step)
            dx *= 500/(h_array[-1] - h_array[-2])
        return h_array, x_array

    def calcAscent(self, C1, T, beta, C2, m0, cDrag, A_front, propburn, dt = 1e-3, grav_turn = [1, 1000], h_cutoff = 100000, mf = 0, steer_rate = 0.017452006980803, throttle_rate = 0.5, a_lim = 100.0):
        """
        Calculate the Ascent Portion of a Launch Vehicle, based on a target trajectory of the form h(x) = (0.25 LEOalt (x))^(c), c = 0..1 
        The LV will ascend vertically for 1km, proceed to initiate a gravity turn by 3° and then start attempting to follow the target launch trajectory
        The LV will continue ascending until reaching either the cutoff altitude or the cutoff mass
        The LV data can be provided in array like format, to allow for mid ascent staging

        Inputs:
            LEOalt: Altitude of Target LEO, affects the target flight path (higher LEOalt = steeper ascent)
            T: Maximal Thrust, or array like of Maximal Thrusts per stage
            beta: latitude of launch site
            C: constant "c" in target launch trajectory
            m0: Start mass, or array like of start masses per stage
            cDrag: Drag coefficient, or array like of drag coefficients per stage. Can be provided as value, if equal for all stages
            A_front: Frontal LV area, or array like of frontal LV area per stage. Can be provided as value, if equal for all stages
            propburn: mass flow rate, or array like of mass flow rates per stage
            dt: time step. Default = 1e-3
            h_cutoff: cutoff altitude. Default = 100000
            mf: cutoff mass. Default = 0
            steer_rate: maximum TVC deflection rate. Default = 1° per second
            throttle_rate: maximum Throttle up/down rate. Default = 0.5 (50% per second)
            a_lim: maximum accelration. Default = 100 m/s^2 (approx. 10g)
            grav_turn: Initiation of gravity turn, [1, altitude] OR [0, speed]. Default = [1, 1000]
        Returns:
            dictionary of ascent data for each time step with entries:
                "t": Time
                "h": Height
                "x": Downrange
                "D": Drag
                "alpha": Steering Angle
                "v": Speed
                "a": Accelration
                "T": Thrust
                "m": Mass
                "q": Dynamic pressure
                "rot_loss": Rotational Loss
                "grav_loss": Gravitational Loss
                "steer_loss": Steering Loss
                "drag_loss": Drag Loss
                "tot_loss": Total Loss
        """ 
        #x = downrange distance
        #Assume trajectory with some vertical ascent, then perform gravity turn to gamma = 87°
        #Then following a fligth path of h(x) = (C2 (x + d))^(C1) + h_offset
        #With d being set such that the flight path at the end of the gravity turn smoothly transitions into profile
        #By Ensuring that gamma = arctan(h'(x_EndOfGravTUrn)) = 87°
        #X is offset by x_87 given by the condition: arctan(h'(x_87)) = 87°
        #Thus for h(x) = (C2 * (x+x_87))^(C1), arctan(h'(0)) = 87°
        #At the End of the Gravity Turn the Offset is then modified by the current downrange such that d = x_87 - x_EndOfGravTurn
        #Thus h(x) = (C2 * (x+x_87-x_EndOfGravTurn))^(C1) thus arctan(h'(x_EndOfGravTurn)) = 87° and a smooth transition of the flight profile is ensured

        n=0
        H = 6700
        pass_flag = False
        gamma = np.pi/2
        x = 1e-6 #Function not defined at x = 0, thus start "infinitessimally" close to 0
        dx = 10
        #Bisection Determination of x_87
        while abs(gamma - np.deg2rad(87)) > 1e-5:
            while gamma > np.deg2rad(87):
                gamma = np.arctan(C1*(C2**(C1))*x**(C1-1))
                x+=dx
            x -= 2*dx
            dx *= 0.5
            gamma = np.arctan(C1*(C2**(C1))*x**(C1-1))
        x_offset = x
        beta = np.deg2rad(beta)
        rot_loss = - 0.464*np.cos(beta)
        x_step = 0
        if grav_turn[0] == 0:
            h_gt = 0
            v_gt = grav_turn[1]
        else:
            h_gt = grav_turn[1]
            v_gt = 0
        try:
            m_step = m0[n]
        except:
            m_step = m0
        try:
            m_cutoff = mf[n]
        except:
            m_cutoff = mf
        try:
            m_dot = propburn[n]
        except:
            m_dot = propburn
        try:
            cD = cDrag[n]
        except:
            cD = cDrag
        try:
            A = A_front[n]
        except:
            A = A_front
        try:
            T_step = T[n]
            T_Max = T[n]
        except:
            T_step = T
            T_Max = T
        if m_dot > 0:
            m_dot *= -1

        v_step = 0
        gamma_step = np.pi/2
        h_step = 0
        steer_loss = 0
        drag_loss = 0
        grav_loss = 0
        D_step = 0
        gamma_tar = np.pi/2
        alpha_step = 0
        t = 0

        ascent_data = {}
        ascent_data["h"]  = [0]
        ascent_data["x"] = [0]
        ascent_data["D"] = [0]
        ascent_data["alpha"] = [0]
        ascent_data["v"] = [0]
        ascent_data["a"] = [0]
        ascent_data["T"] = [T_step]
        ascent_data["m"] = [m_step]
        ascent_data["t"] = [0]
        ascent_data["rot_loss"] = [0]
        ascent_data["grav_loss"] = [0]
        ascent_data["steer_loss"] = [0]
        ascent_data["tot_loss"] = [0]
        ascent_data["drag_loss"] = [0]
        ascent_data["q"] = [0]
        
        while h_step <= h_cutoff:
            #Debugging
            if h_step == 0:
                temp =0 
            if h_step > 1000 and h_step < 1001:
                temp = 0
            if h_step > 2000 and h_step < 2001:
                temp = 0
            rho_step = 1.225 * np.exp(-h_step/H)
            if h_step > 1000 and gamma_step > 1.5184364:
                gamma_tar = 1.5184364 #87°
            elif h_step > 1000:
                temp = 0
            alpha_step = self.calcAlpha(gamma_step, gamma_tar, T_step, m_step, v_step, h_step + 6378000, dt, alpha_step, steer_rate)
            q = rho_step*0.5*v_step**2
            D_step = q*cD*A
            steer_loss += 2*T_step*np.sin(alpha_step*0.5)**2/m_step * dt
            drag_loss += D_step/m_step * dt
            grav_loss += 9.81*np.sin(gamma_step) * dt
            v_step += ((T_step*np.cos(alpha_step)-D_step)/m_step - 9.81*np.sin(gamma_step))*dt
            a_step = ((T_step*np.cos(alpha_step)-D_step)/m_step - 9.81*np.sin(gamma_step))
            if a_step >= a_lim:
                T_step -= T_Max*throttle_rate*dt
            else:
                T_step += T_Max*throttle_rate*dt
            if T_step > T_Max:
                T_step = T_Max
            h_step += v_step * np.sin(gamma_step) * dt
            x_step += abs(v_step * np.cos(gamma_step) * dt) 
            if not (h_step > h_gt and v_step > v_gt):
                x_step = 0
                gamma_step = np.pi/2
                gamma_tar = np.pi/2
            else:
                #initiate Gravity Turn
                gamma_step += (T_step * np.sin(alpha_step)/(m_step*v_step) - 9.81/v_step*np.cos(gamma_step) + v_step/(6378000+h_step)*np.cos(gamma_step)) * dt
                #At End of Gravity Turn, determine offsets for smooth transition into profile
                if gamma_step <= np.deg2rad(87) and not pass_flag:
                    x_offset -= x_step
                    #Ensure only done once
                    pass_flag = True
                gamma_tar = np.arctan(C1*(C2**(C1))*(x_step + x_offset)**(C1-1))
            m_step += m_dot*dt

            t += dt
            if m_step < m_cutoff:
                n += 1
                try:
                    m_step = m0[n]
                    T_step = 0
                    T_Max = T[n]
                    m_dot = propburn[n]
                    try:
                        A = A_front[n]
                        cD = cDrag[n]
                    except:
                        pass
                    print("Stage Seperation at" +  str(h_step))
                except:
                    print("Out of Stage Fuel//End of Ascent")
                    h_cutoff = 0

            ascent_data["h"].append(h_step)
            ascent_data["x"].append(x_step)
            ascent_data["D"].append(D_step)
            ascent_data["v"].append(v_step)
            ascent_data["a"].append(a_step)
            ascent_data["alpha"].append(alpha_step)
            ascent_data["T"].append(T_step)
            ascent_data["m"].append(m_step)
            ascent_data["q"].append(q)
            ascent_data["t"].append(t)
            ascent_data["rot_loss"].append(rot_loss)
            ascent_data["grav_loss"].append(drag_loss)
            ascent_data["steer_loss"].append(steer_loss)
            ascent_data["drag_loss"].append(drag_loss)
            ascent_data["tot_loss"].append(rot_loss+grav_loss+steer_loss)
        return  ascent_data

    def calcAlpha(self, gamma_0, gamma_tar, T, m, v, r, dt, alpha_0, steer_rate):
        #iteratively determine the steering angle in order to best follow the profile
        if alpha_0 >= 0.436300: 
            return 0.34904
        elif alpha_0 <= -0.34904:
            return -0.436300
        try:
            alpha = np.arcsin(m*v/T *((gamma_tar-gamma_0)/dt + np.cos(gamma_0)*(9.81/v - v/r)))
        except RuntimeWarning:
            if gamma_tar - gamma_0 > 0:
                alpha = alpha_0 + steer_rate * dt
            else:
                alpha = alpha_0 - steer_rate * dt
            return alpha
        except ZeroDivisionError:
            return alpha_0
        if alpha > alpha_0 + steer_rate*dt:
            alpha = alpha_0 + steer_rate * dt
        elif alpha < alpha_0 - steer_rate*dt:
            alpha = alpha_0 - steer_rate * dt
        elif math.isnan(alpha):
            if gamma_tar - gamma_0 > 0:
                alpha = alpha_0 + steer_rate * dt
            else:
                alpha = alpha_0 - steer_rate * dt
            return alpha
        return alpha

    def calcBoosterCont_OptimalCore(self, n, n_booster, isp_booster, isp_core, m_pl, MFR, mu_booster = 0.12, mu_core = 0.12, delv = 9000, limit = 1e6, **kwargs):
        """
        Unstable
        Size a rocket with boosters
        Inputs:
            n: Number of Stages of core rocket
            n_booster: Number of boosters
            isp_booster: Efficiency of boosters
            isp_core: isp of core rocket 
                        OR
                      list of isps for core rocket
            m_pl: Mass of Payload
            MFR: Mass flow ratio (m_dot_core / m_dot_booster)
            mu_booster: Structure Factor booster (default = 0.12)
            mu_core: Strucutre factor core rocket (default = 0.12)
            delv: Target delta v (default = 9000 m/s)
            limit: Divergence limit (default = 1e6 kg)
        optional inputs:
            "m_f": List of fuel masses for core rocket 
            "range": Range of mass fractions (m_booster/m_core_stage_1) [imin, imax]
        Outputs:
            mininum total launch mass
            minimum booster mass
            minimum core rocket mass
            Mass Fraction Booster/Core Stage 1
        """
        if "range" in kwargs:
            i_min = kwargs["range"][0]
            i_max = kwargs["range"][1]
        else:
            i_max = int(100 * 1/MFR)
            i_min = 1
        delv_tar = delv
        m_tot_min = limit
        isp_weighted = (isp_core * MFR + isp_booster*1/MFR * n_booster)/(n_booster*1/MFR+1*MFR)
        for i in range(i_min, i_max+1):
            while True:
                m0_core, factor, m_f = self.calcPoint(n, isp_core, m_pl, mu = mu_core, delv = delv_tar, limit = limit)
                m0_core = self.MassSplit(n, m0_core, m_pl, mu_core, factor)[0]
                m0_core_1 = m0_core[0]
                m0_booster = m0_core_1 * i * 1e-2
                dv_tot = self.calcDelV_par(n, m0_booster, mu_booster, m0_core[0:n], n_booster, m_pl, MFR, isp_weighted, isp_core, m_f = m_f)
                if dv_tot <= 1.001 * delv and dv_tot >= 0.999 * delv:
                    break
                delv_tar -= dv_tot - delv
            m_tot = sum(m0_core) + n_booster*m0_booster
            if m_tot < m_tot_min:
                m_tot_min = m_tot
                m0_booster_min = m0_booster
                m0_core_min = m0_core[0:n+1]
                i_min = i
        return m_tot_min, m0_booster_min, m0_core_min, i_min
    
    def calcBoosterCont_FixedCore(self, n_booster, isp_booster, isp_core, MFR, m0_core, m_pl, mu_booster = 0.12, mu_core = 0.12, delv = 9000, **kwargs):
        """
        Size boosters for a fixed core rocket
        Inputs:
            n_booster: Number of boosters
            isp_booster: Efficiency of boosters
            isp_core: isp of core rocket 
                        OR
                      list of isps for core rocket
            m_pl: Mass of Payload
            MFR: Mass flow ratio (m_dot_core / m_dot_booster)
            m0_coore: List of Stage masses of core rocket (excl. payload)
            mu_booster: Structure Factor booster (default = 0.12)
            mu_core: Strucutre factor core rocket (default = 0.12)
            delv: Target delta v (default = 9000 m/s)
            limit: Divergence limit (default = 1e6 kg)
        Optional Inputs:
            "m_f": List of fuel masses for core rocket 
            "range": Range of mass fractions (m_booster/m_core_stage_1) [imin, imax]
        Outputs:
            mininum total launch mass
            minimum booster mass
            minimum core rocket mass
            Mass Fraction Booster/Core Stage 1
        """
        #MFR = m_dot_core / m_dot_booster
        n = len(m0_core)
        if "range" in kwargs:
            i_min = kwargs["range"][0]
            i_max = min(kwargs["range"][1], int(100 * 1/MFR))
        else:
            i_max = int(100 * 1/MFR)
            i_min = 1
        if "m_f" in kwargs:
            mp_core = kwargs["m_f"][0:n]
        else:
            mp_core = [val*(1-mu_core) for val in m0_core]
        try:
            isp_core_0 = isp_core[0]
        except TypeError:
            isp_core_0 = isp_core
        delv_core = self.calcDelV(n, m0_core, m_pl, isp_core, m_f = mp_core)
        if delv_core >= delv:
            return 0
        m0_core_1 = m0_core[0]
        isp_weighted = (isp_core_0 * MFR + isp_booster*1/MFR * n_booster)/(n_booster*1/MFR+1*MFR)
        dv_tot = delv_core
        for i in range(i_min, i_max+1):
            m0_booster = m0_core_1 * i * 1e-2
            dv_tot_ = self.calcDelV_par(n, m0_booster, mu_booster, m0_core[0:n], n_booster, m_pl, MFR, isp_weighted, isp_core, m_f = mp_core)
            if dv_tot_ < dv_tot:
                break
            dv_tot = dv_tot_
            if dv_tot >= delv:
                break
        return m0_booster, i, dv_tot

    def calcBoosterDisc_FixedCore(self, isp_core, isp_booster, m0_core, m0_booster,m_pl, MFR, mu_core = 0.12, mu_booster = 0.12 , delv = 9000, booster_align = [2,3,4,5,6,7,8], **kwargs):
        """
        For fixed core rocket and fixed booster size, determine configuration with minimum number of boosters to reach target dv
        Inputs:
            n_booster: Number of boosters
            isp_core: isp of core rocket 
                        OR
                      list of isps for core rocket
            isp_booster: Efficiency of boosters
            m0_coore: List of Stage masses of core rocket (excl. payload)
            m0_booster: Booster mass
            m_pl: Mass of Payload
            MFR: Mass flow ratio (m_dot_core / m_dot_booster)
            mu_core: Strucutre factor core rocket (default = 0.12)
            mu_booster: Structure Factor booster (default = 0.12)
            delv: Target delta v (default = 9000 m/s)
            booster_allign: List of available booster alignments (default: [2,3,4,5,6,7,8])
        Optional Inputs:
            "m_f": List of fuel masses for core rocket 
        Outputs:
            mininum total launch mass
            minimum booster mass
            minimum core rocket mass
            Mass Fraction Booster/Core Stage 1
        """
        n = len(m0_core)
        if "m_f" in kwargs:
            mp_core = kwargs["m_f"][0:n]
        else:
            mp_core = [val*(1-mu_core) for val in m0_core]
        m0_core_1 = m0_core[0]
        try:
            isp_core_0 = isp_core[0]
        except TypeError:
            isp_core_0 = isp_core
        if int(100 * 1/MFR)*m0_core_1 < m0_booster:
            raise Exception("Given MFR and Booster Mass results in Booster Burn Time exeeding Core Stage 1 Burn Time")
        for i in booster_align:
            isp_weighted = (isp_core_0 * MFR + isp_booster*1/MFR * i)/(i*1/MFR+1*MFR)
            dv_tot = self.calcDelV_par(n, m0_booster, mu_booster, m0_core[0:n], i, m_pl, MFR, isp_weighted, isp_core, m_f = mp_core)
            if dv_tot >= delv:
                return i, dv_tot
        return i, dv_tot
    
    def calcDelV_par(self, n, m0_booster, mu_booster, m0_core, n_booster, m_pl, MFR, isp_weighted, isp_core, mu_core = 0.12,  **kwargs):
            """
            Calculate deltaV of a parallel staged Rocket
            Inputs:
                n: Number of stages
                m0_booster: Booster mass+
                mu_booster: structure fator of booster
                m0_core: List of Stage masses of core rocket (excl. payload)
                n_booster: Number of boosters
                m_pl: payload
                MFR: Mass Flow Ratio (m_dot_core/m_dot_booster)
                isp_weighted: Weighted Isp of parallel stage
                isp_core: isp of core rocket
                            OR
                          List of stage Isps of core rocket
                mu_core: strucure factor of core rocket
            Optional Inputs:
                "m_f": List of fuel masses for core rocket 
            Outputs:
                mininum total launch mass
                minimum booster mass
                minimum core rocket mass
                Mass Fraction Booster/Core Stage 1
            """
            if "m_f" in kwargs:
                mp_core = kwargs["m_f"][0:n]
            else:
                mp_core = [val*(1-mu_core) for val in m0_core]
            mf_booster = m0_booster * mu_booster
            m0_LV = sum(m0_core) + m0_booster * n_booster + m_pl
            mf_core = sum(m0_core) + m_pl - (m0_booster - mf_booster)*MFR 
            mf_LV = mf_booster*n_booster + mf_core
            dv_BoosterStage = self.Tsiolkowsky(isp_weighted, m0_LV, mf_LV)
            m =  m0_core[0:n]
            m[0] -= (m0_booster - mf_booster)*MFR 
            mp_core[0] -=  (m0_booster - mf_booster)*MFR 
            n_ = n
            if mp_core[0] <= 0:
                try:
                    isp_core = isp_core[1:n]
                except TypeError:
                    pass
                m = m[1:n]
                mp_core = mp_core[1:n]
                n_=n-1
            dv_CoreStages = self.calcDelV(n_, m, m_pl, isp_core, m_f = mp_core)
            return dv_CoreStages + dv_BoosterStage


    def InjCalc(self, dens, K, dPrel, pc, mdot, A):
        inja = [[],[]]
        ello = []
        injv = [[],[]]
        numbel = []
        injdp = [[],[]]
        dP_fu_0 = dPrel[0]*pc*100000
        dP_ox_0 = dPrel[1]*pc*100000
        v_fu_0 = mdot[0]/(A*dens[0])
        v_ox_0 = mdot[1]/(A*dens[1])
        for N in range(50,500):
            A_inj_fu = mdot[0]/(N*(np.sqrt(dP_fu_0*2*dens[0]/K)+ dens[0] * v_fu_0))
            A_inj_ox = mdot[1]/(N*(np.sqrt(dP_ox_0*2*dens[1]/K)+ dens[1] * v_ox_0))
            elem_load = sum(mdot)/N
            v_fu = mdot[0]/(A_inj_fu*N*dens[0])
            v_ox = mdot[1]/(A_inj_ox*N*dens[1])
            mu_ox = 195.3e-6
            mu_fu = 117.2e-6
            ks = 30e-6
            l = 0.03
            inja[0].append(A_inj_fu)
            inja[1].append(A_inj_ox)
            ello.append(elem_load)
            injv[0].append(v_fu)
            injv[1].append(v_ox)
            numbel.append(N)
        return inja, ello, injv, numbel

    def solveLam(self, ks, d, Re):
        A = -2*np.log10((ks/d)/3.7 + 12/Re)
        B = -2*np.log10((ks/d)/3.7 + 2.51*A/Re)
        C = -2*np.log10((ks/d)/3.7 + 2.51*B/Re)
        return (A-(B-A)**2/(C-2*B+A))**-2
    def CoolingCalc(self):
        TS
    def TMCalc(self, p_pb, of_pb, TET, cp, gamma):
        TM.run(p_pb, of_pb, TET, cp, gamma)

    #//////////////////////////////////////
    #//////////////////////////////////////
    #//////////////////////////////////////
    #//////////////////////////////////////
    #//////////////////////////////////////
    #//////////////////////////////////////
    #OBSOLETE
    #//////////////////////////////////////
    #//////////////////////////////////////
    #//////////////////////////////////////
    #//////////////////////////////////////
    #//////////////////////////////////////
    #//////////////////////////////////////
    #//////////////////////////////////////

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