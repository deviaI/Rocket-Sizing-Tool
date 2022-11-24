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
import os.path
import sys

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
        n = 2
        pl = 4000
        dv_core = 11000
        dv_tot = 9500
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
if __name__ == '__main__':
    
    main()