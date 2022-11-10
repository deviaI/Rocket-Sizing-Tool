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
    
    
    losses, data = Calc.calcAscent(400000, 224300, 0, 0.5, 12500, 0.7, 4.52, -85)
    print(losses)
    Plt.plot2D(data["x"], data["h"], "downrange", "heigh")
    Plt.plot2D(data["x"], data["alpha"], "downrange", "alpha")
    return 0
    if getattr(sys, 'frozen', False):
        exe_path = os.path.dirname(sys.executable)
    else:
        exe_path = None
    gui = GUI(plotter = Plt,exporter= Exp, calculator = Calc, exeDir = exe_path)
    gui.run()




    
    # #Example Calculations 

    # print(Calc.optimiseMassRatio(2, 12500, [311,343], 250, 0.1))

    # result = Calc.calcDelV(2, 12500, 250, [311, 343])  #Calculate Delta V of 12.5 ton 2 Stage rocket with Isps of 311 and 343
    # print(result)

    # result = Calc.calcPoint(2, [311, 343], 250, 0.1, delv = 9000) #Calculate theoretical launch mass to lift 250kg payload to LEO with a 
    #                                                                #2 Stage Rocket with Isp 311 and 343 engines, structural factor of 10%
    #                                                                # and optimal stage to stage mass ratio
    # print("Calculating Range")
    # results = Calc.calcRange(2, "mu", [311, 343], 250, [0.05, 0.17]) #Calculate theoretical launch mass for Rocket given above across a range of mu from 0.05 to 0.17
    # print("Done!")
    # Plt.plot2D(results[1], results[0], xLab = "Mu", yLab = "Launch Mass")

    # fName = "Isp1=311;Isp2=343;Mu=0.05:0.17"
    # Exp.ExportData(results, fName, "csv")  #Save results in csv file
    # Description = "Isp 1 = 311s;Isp 2 = 343s, stage to stage mass ratio: Optimized;Column 1: Mass, Column 2: Mu, Column 3: stage to stage mass ratio"
    # Exp.AppendData(Description, fName + ".csv") #Append Data description to File



if __name__ == '__main__':
    
    main()
