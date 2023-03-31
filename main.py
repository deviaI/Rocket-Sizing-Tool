#!/usr/bin/env python
'''
# -*- coding: utf-8 -*-
Created on  2022/10/26 09:23:69
@author  Devial
'''
from utils.Calculators import Calculator
from utils.Exporter import Exporter
from utils.Plotters import Plotter
from utils.GUI import GUI
from utils.Lookup import LookUpTable
import os.path
import sys
import csv

def main():
    gui = False
    exe_path = None
    Calc = Calculator()
    if getattr(sys, 'frozen', False):
        exe_path = os.path.dirname(sys.executable)
        Exp = Exporter(os.path.join(exe_path, "data"))
        Plt = Plotter(os.path.join(exe_path, "data"))
        try:
            if not os.path.isfile(os.path.join(exe_path, "files", ".datadir")):
                with open(os.path.join(exe_path, "files", ".datadir"), "w") as f:
                    f.write(os.path.join(exe_path, "data"))
                    f.close()
        except:
            pass
        gui = GUI(plotter = Plt,exporter= Exp, calculator = Calc, exeDir = exe_path)
        gui.run()
        return 0
    if gui == True:
        base = os.path.dirname(__file__)
        Exp = Exporter(os.path.join(base, "data"))
        Plt = Plotter(os.path.join(base, "data"))
        try:
            if not os.path.isfile(os.path.join(base, "utils", "files", ".datadir")):
                with open(os.path.join(base, "utils", "files", ".datadir"), "w") as f:
                    f.write(os.path.join(base, "data"))
                    f.close()
        except:
            pass
        gui = GUI(plotter = Plt,exporter= Exp, calculator = Calc, exeDir = exe_path)
        gui.run()
    else:
        base = os.path.dirname(__file__)
        Plt = Plotter(os.path.join(base, "data"))
        Exp = Exporter(os.path.join(base, "data"))

        

        """
        ARCHIVED CALCULATIONS
        rez = Calc.calcRange(1, "Mu", 400, 200, [0.03,0.07], 9300)
        m = rez[0]
        isp = rez[1]
        Plt.plot2D(isp, m, xLab = "Mu [-]", yLab = "Mass [kg]")
                #DVs: 
        #SL2: 13.28
        #GTO: 11.8
        #Venus Escape: 12.8
        #LEO: 9.3
        #SSO: 9.5
        n = 3
        pl = 15000
        dv_target = 9200
        isp = [320, 350, 440]
        mu_core = 0.12  	
        mu_booster = 0.1
        MFR = 1.5
        m = [103466.80868422716, 45525.39582105995, 20031.17416126638, pl]
        m_fuel = [91050.7916421199, 40062.34832253276, 17627.433261914415]
        m_booster = 68288.09373158992
        print(Calc.calcAscent_Ideal_DV(393000))
        print(Calc.calcDelV(n, m[0:n], pl, isp[0:n], m_f = m_fuel[0:n]))
        print(Calc.calcBoosterDisc_FixedCore(isp[0:n], 320, m[0:n], m_booster, pl, MFR, mu_core, mu_booster, dv_target, m_f = m_fuel[0:n], booster_align=[2,4]))


        base = os.path.dirname(__file__)
        Plt = Plotter(os.path.join(base, "data"))
        Exp = Exporter(os.path.join(base, "data"))
        # TMLU = LookUpTable()
        # TMLU.showKeyList()
        # Eff = TMLU.returnTable(0)
        # Speed = TMLU.returnTable(1, "AL2219-T87")
        # print(TMLU.SpeedLookup("AL2219-T87", 550))
        # Plt.plot2D(Speed[0], Speed[1], [200, 1310], [100,550])
        
        print(Calc.calcAscent_Ideal_DV(400000))
        rez = Calc.calcAscent(0.4, 2100000, 0, 2500000, 151000, 0.75, 4.5, 600, mf = 60000)
        time = rez["t"]
        downrange = rez["x"]
        alt = rez["h"]
        vel = rez["v"]
        q = rez["q"]
        T = rez["T"]
        a = rez["a"]
        m = rez["m"]
        D = rez["D"]
        alpha = rez["alpha"]
        tot_loss = rez["tot_loss"]
        for i in range(0, len(a)):
            a[i] /= 9.81
        Plt.plot2D(dataX = time, dataY = alt, yLab="Altitude[m]", xLab="Time", filename= "alt_t", savefile=1, show=0)
        Plt.plot2D(dataX = downrange, dataY = alt, yLab="Altitude[m]", xLab="Downrange", filename= "alt_x", savefile=1, show=0)
        Plt.plot2D(dataX = time, dataY = vel, yLab="Velocity[m/s]", xLab="Time", filename= "v_t", savefile=1, show=0)
        Plt.plot2D(dataX = time, dataY = q, yLab="Dynamic Pressure[Pa]", xLab="Time", filename= "q_t", savefile=1, show=0)
        Plt.plot2D(dataX = time, dataY = T, yLab="Thrust[N]", xLab="Time", filename= "T_t", savefile=1, show=0)
        Plt.plot2D(dataX = time, dataY = a, yLab="Gs[-]", xLab="Time", filename= "a_t", savefile=1, show=0)
        Plt.plot2D(dataX = time, dataY = D, yLab="Drag[N]", xLab="Time", filename= "D_t", savefile=1, show=0)
        Plt.plot2D(dataX = time, dataY = m, yLab="Mass[kg]", xLab="Time", filename= "m_t", savefile=1, show=0)
        Plt.plot2D(dataX = time, dataY = alpha, yLab="Steering Angle[rad]", xLab="Time", filename= "alpha_t", savefile=1, show=0)
        Plt.plot2D(dataX = time, dataY = tot_loss, yLab="Ascent Losses[m/s]", xLab="Time", filename= "loss_t", savefile=1, show=0)
        Exp.ExportData(data=alt, fType = ".csv", fName = "Altitude")
        Exp.ExportData(data=time, fType = ".csv", fName = "Time")

        # h, x = Calc.gen_ascent_path_preview(0.4, 2500000, 500, 150000)
        # x = [val*1e-3 for val in x]
        # karman_line = [100000 for k in x]
        # Plt.plot2D(x, h, dataY_List = [karman_line], xlim = [-max(x)*0.1, max(x)], xLab = "Downrange [km]", yLab = "Altitude [m]", data_Labels = ["Ascent Profile" ,"Karman Line"])


        # [inja, ello, injv, numel] = Calc.InjCalc([7.146, 1141], 1.5, [0.05, 0.1], 222, [26.32, 73.68], 0.0404)
        # Plt.plot2D(numel, inja[0], yLab = "injector Area [m^2]", xLab = "number injectors", dataY_List=[inja[1]], data_Labels=["Fuel", "Oxidizer"])
        # Plt.plot2D(numel, injv[0], yLab = "injection velocity[m/s]",  xLab = "number injectors", dataY_List=[injv[1]], data_Labels=["Fuel", "Oxidizer"])
        # Plt.plot2D(numel, ello, yLab = "element loading (kg/s)", xLab = "number injectors")
        # mdotO/mdotF = 2.8
        # mdotO/mdot-mdotO = 2.8
        # mdot -mdotO(2.8)=mdotO
        # mdot2.8 = mdotO + mdotO2.8
        #mdot 2.8 = mdotO*(1+2.8
        # mdotO = mdot*2.8/(1+2.8)
        # print(Calc.calcAscent_Ideal_DV(250000))
        # m, factor, fuel = Calc.calcPoint(n, isp, pl, mu_core, dv_core, limit = 1e9)
        # print(m)
        # m, fuel = Calc.MassSplit(n, m, pl, mu_core, factor)
        # print(m)
        # print(fuel)
        # print(Calc.calcDelV(2, m[0:2], 3600, [320, 350]))
        # print(Calc.calcBoosterCont_FixedCore(4, 320, isp, MFR, m[0:n], pl, 0.1, mu_core, dv_tot, m_f = fuel))
        #print(Calc.calcBoosterCont_FixedCore(6, 360, 360, 2, m[0:2], 250, delv = 12000, m_f = fuel))
        #print(Calc.calcBoosterCont_OptimalCore(2, 2, 250, 360, 250, 1))
        #print(Calc.calcPoint(2, 360, 250))



        
        n = 3
        pl = 4000
        dv_core = 11000
        dv_tot = 12800
        isp = [320, 350, 440]
        mu_core = 0.12  	
        mu_booster = 0.1
        MFR = 1.5
        m = [103466.80868422716, 45525.39582105995, 20031.17416126638, 4000]
        m_fuel = [91050.7916421199, 40062.34832253276, 17627.433261914415]
        m_booster = 68288.09373158992
        print(Calc.calcAscent_Ideal_DV(393000))
        print(Calc.calcDelV(n, m[0:n], pl, isp[0:n], m_f = m_fuel[0:n]))
        print(Calc.calcBoosterDisc_FixedCore(isp[0:n], 320, m[0:n], m_booster, pl, MFR, mu_core, mu_booster, dv_tot, m_f = m_fuel[0:n]))
        # print(Calc.calcAscent_Ideal_DV(250000))
        m, factor, fuel = Calc.calcPoint(n, isp, pl, mu_core, dv_core, limit = 1e9)
        # print(m)
        m, fuel = Calc.MassSplit(n, m, pl, mu_core, factor)
        print(m)
        # # print(fuel)
        # # print(Calc.calcDelV(2, m[0:2], 3600, [320, 350]))
        # print(Calc.calcBoosterCont_FixedCore(4, 320, isp, MFR, m[0:n], pl, 0.1, mu_core, dv_tot, m_f = fuel))
        #print(Calc.calcBoosterCont_FixedCore(6, 360, 360, 2, m[0:2], 250, delv = 12000, m_f = fuel))
        #print(Calc.calcBoosterCont_OptimalCore(2, 2, 250, 360, 250, 1))
        #print(Calc.calcPoint(2, 360, 250))
        # Calc.CoolingCalc()
        # Calc.TMCalc(30000000, 0.38, 1169, 9221, 1.1361)
        """
        
if __name__ == '__main__':
    main()