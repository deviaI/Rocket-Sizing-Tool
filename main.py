#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on  2022/11/11 17:13:19
@author  Devial
'''

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
        m, factor, fuel = Calc.calcPoint(2, 360, 250, delv = 10000)
        m, fuel = Calc.MassSplit(2, m, 250, 0.12, factor)
        #print(Calc.calcBoosterDisc_FixedCore(360, 360, m[0:2], 2450, 250, 2, delv = 12000, m_f = fuel))
        # print(Calc.calcBoosterCont_FixedCore(6, 360, 360, 2, m[0:2], 250, delv = 12000, m_f = fuel))
        print(Calc.calcBoosterCont_OptimalCore(2, 2, 250, 360, 250, 1))
        #print(Calc.calcPoint(2, 360, 250))
if __name__ == '__main__':
    
    main()
