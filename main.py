#!/usr/bin/env python
'''
# -*- coding: utf-8 -*-
Created on  2022/10/26 09:23:69
@author  Devial
'''
from utils.Calculators import Calculator
from utils.Exporter import Exporter
from utils.Plotters import Plotter
<<<<<<< HEAD
import numpy as np
=======
from utils.GUI import GUI
>>>>>>> GUI
import os.path
import math

def main():

    base = os.path.dirname(__file__)
    Calc = Calculator()
    Exp = Exporter(os.path.join(base, "data"))
    Plt = Plotter(os.path.join(base, "data"))
    gui = GUI(plotter = Plt,exporter= Exp, calculator = Calc)
    gui.run()  


    
    # #Example Calculations

    # print(Calc.optimiseMassRatio(2, 12500, [311,343], 250, 0.1))

    # result = Calc.calcDelV(2, 12500, 250, [311, 343])  #Calculate Delta V of 12.5 ton 2 Stage rocket with Isps of 311 and 343
    # print(result)

<<<<<<< HEAD
    result = Calc.calcPoint(2, [311, 343], 250, 0.1, delv = 9000) #Calculate theoretical launch mass to lift 250kg payload to LEO with a 
                                                                   #2 Stage Rocket with Isp 311 and 343 engines, structural factor of 10%
                                                                   # and optimal stage to stage mass ratio
    # print("Calculating Range")
    # resultsA = Calc.calcRange(2, "mu", [280, 311], 250, [0.05, 0.17]) #Calculate theoretical launch mass for Rocket given above across a range of mu from 0.05 to 0.17
    # resultsB = Calc.calcRange(2, "mu", [311, 343], 250, [0.05, 0.17]) #Calculate theoretical launch mass for Rocket given above across a range of mu from 0.05 to 0.17
    # resultsC = Calc.calcRange(2, "mu", [343, 410], 250, [0.05, 0.17]) #Calculate theoretical launch mass for Rocket given above across a range of mu from 0.05 to 0.17
    # print("Done!")
    # Plt.plot2D(resultsA[1], resultsA[0], xLab = "Mu", yLab = "Launch Mass", dataY_List=(resultsC[0], resultsB[0]), data_Labels = ("280/311", "311/43", "343410"))
    n=0
    len = 5
    rez = np.zeros((len, 50))
    isps = np.linspace(280, 410, len)
    for isp in isps:
        print("Progress: " + str(n) +  " out of " + str(len))
        results = Calc.calcRange(1, "mu", isp, 250, [0.05, 0.11], numSteps= 50) #Calculate theoretical launch mass for Rocket given above across a range of mu from 0.05 to 0.17
        rez[len-n-1]  = results[0]
        n+=1
    for i in range(0,len):
        for k in range(0,50):
            if rez[i][k] >= 1e90:
                rez[i][k] = 1e6
            rez[i][k] = math.log(rez[i][k], 10)
    isps = np.flip(isps)
    Plt.colourMap(results[1], isps, rez)
    Plt.contour(results[1], isps, rez)
=======
    # result = Calc.calcPoint(2, [311, 343], 250, 0.1, delv = 9000) #Calculate theoretical launch mass to lift 250kg payload to LEO with a 
    #                                                                #2 Stage Rocket with Isp 311 and 343 engines, structural factor of 10%
    #                                                                # and optimal stage to stage mass ratio
    # print("Calculating Range")
    # results = Calc.calcRange(2, "mu", [311, 343], 250, [0.05, 0.17]) #Calculate theoretical launch mass for Rocket given above across a range of mu from 0.05 to 0.17
    # print("Done!")
    # Plt.plot2D(results[1], results[0], xLab = "Mu", yLab = "Launch Mass")
>>>>>>> GUI

    # fName = "Isp1=311;Isp2=343;Mu=0.05:0.17"
    # Exp.ExportData(results, fName, "csv")  #Save results in csv file
    # Description = "Isp 1 = 311s;Isp 2 = 343s, stage to stage mass ratio: Optimized;Column 1: Mass, Column 2: Mu, Column 3: stage to stage mass ratio"
    # Exp.AppendData(Description, fName + ".csv") #Append Data description to File



if __name__ == '__main__':
    
    main()
