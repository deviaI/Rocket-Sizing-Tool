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

    base = os.path.dirname(__file__)    
    Calc = Calculator()
    Exp = Exporter(os.path.join(base, "data"))
    Plt = Plotter(os.path.join(base, "data"))
    if getattr(sys, 'frozen', False):
        exe_path = os.path.dirname(sys.executable)
    else:
        exe_path = None
    gui = GUI(plotter = Plt,exporter= Exp, calculator = Calc, exeDir = exe_path)
    gui.run()

if __name__ == '__main__':
    
    main()
