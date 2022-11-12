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
    gui = True
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
        data = Calc.calcAscent(400000, 224300, 0, 0.5, 12500, 0.7, 4.5239, -85)
        Exp = Exporter(os.path.join(base, "data"))
        Plt = Plotter(os.path.join(base, "data"))
        Plt.plot2D(data["x"], data["h"], "downrange", "height")
        Plt.plot2D(data["x"], data["v"], "downrange", "speed")
        Plt.plot2D(data["x"], data["a"], "downrange", "acceleration")
        Plt.plot2D(data["x"], data["T"], "downrange", "Thrust")
        Plt.plot2D(data["x"], data["alpha"], "downrange", "alpha")
if __name__ == '__main__':
    
    main()
