#!/usr/bin/env python
'''
# -*- coding: utf-8 -*-
Created on  2022/10/26 09:23:69
@author  Devial
'''

# NON GUI MAIN

from utils.Calculators import Calculator
from utils.Exporter import Exporter
from utils.Plotters import Plotter
import os.path

def main():

    base = os.path.dirname(__file__)    
    Calc = Calculator()
    Exp = Exporter(os.path.join(base, "data"))
    Plt = Plotter(os.path.join(base, "data"))
    
    
    losses, data = Calc.calcAscent(400000, 224300, 0, 0.5, 12500, 0.7, 4.52, -85)
    print(losses)
    Plt.plot2D(data["x"], data["h"], "downrange", "height", savefile=1)
    Plt.plot2D(data["x"], data["alpha"], "downrange", "alpha", savefile=1)
    Plt.plot2D(data["x"], data["v"], "downrange", "speed", savefile=1)
    Plt.plot2D(data["x"], data["a"], "downrange", "accel.", savefile=1)
    return 0

if __name__ == '__main__':
    
    main()
