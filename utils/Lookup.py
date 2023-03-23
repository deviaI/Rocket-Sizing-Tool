#!/usr/bin/env python
'''
# -*- coding: utf-8 -*-
Created on  2023/03/21 09:23:69
@author  Devial
'''

import numpy as np


class LookUpTable():

    def __init__(self):
        self.genTables()

    def genTables(self):
        self.EffTable = np.array([[0.8, 0.85, 0.95, 1.05, 1.15, 1.3, 1.55, 1.7, 2, 2.1, 2.5, 3, 4, 5, 6, 7, 8, 9, 10],
                                  [0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,  0.77, 0.79, 0.8, 0.825, 0.833, 0.85, 0.85, 0.845, 0.837, 0.83, 0.825, 0.816]])
        self.SpeedTable = {}
        temps = np.linspace(320, 1310, 19)
        temps = np.array(temps)
        self.SpeedTable["AL2219-T87"] = np.array([temps[0:5],
                                                [465, 450, 415, 360, 300]])
        self.SpeedTable["TI6-2-4-2"] = np.array([temps[0:10],
                                                [365, 355, 345, 335, 330, 320, 315, 310, 305, 300]])
        self.SpeedTable["A2886"] = np.array([temps[0:10],
                                            [445, 435, 430, 425, 423, 420, 420, 417, 410, 400]])
        self.SpeedTable["IN-100"] = np.array([temps[0:19],
                                                [435,438,440,441,443,445,448,449,450,452,453,450,425,332,285,281,252,231,220]])
        self.SpeedTable["IN-713LC"] = np.array([temps[0:18],
                                                [325,326,327,330,332,335,335,340,345,348,345,340,335,310,282,251,210,180]])
        self.SpeedTable["IN-718"] = np.array([temps[0:14],
                                                [482,481,480,479,475,471,470,470,469,468,440,375,311,275]])
        self.SpeedTable["INY-903"] = np.array([temps[0:14],
                                                [500,490,482,480,482,483,485,481,471,449,400,345,268,218]])
        self.SpeedTable["Waspaloy"] = np.array([temps[0:15],
                                                [439,445,440,438,440,440,441,442,445,470,442,395,330,289,255]])
        self.SpeedTable["HP9-4-30"] = np.array([temps[0:10],
                                                [545,541,538,536,528,515,500,471,450,431]])
        self.SpeedTable["Astroloy"] = np.array([temps[0:15],
                                                [510,510,508,502,501,500,499,498,492,471,458,408,301,298,260]])

    def returnTable(self, table, key = None):
        if table == 0:
            return self.EffTable
        else:
            return self.SpeedTable[key]
    def showKeyList(self):
        for elem in self.SpeedTable:
            print(elem)

    def EffLookup(self, StSpecSpeed):
        if StSpecSpeed > 10 or StSpecSpeed < 0.8:
            raise Exception("Stage Specific Speed Out Of Range")
        if StSpecSpeed >= 9:
            i = 17
            return self.EffTable[1][i] + (self.EffTable[1][i+1] - self.EffTable[1][i])/(self.EffTable[0][i+1] - self.EffTable[0][i]) * (StSpecSpeed - self.EffTable[0][i])
        else:
            for i in range(0,19):
                if self.EffTable[0][i] >= StSpecSpeed:
                    return self.EffTable[1][i] + (self.EffTable[1][i+1] - self.EffTable[1][i])/(self.EffTable[0][i+1] - self.EffTable[0][i]) * (StSpecSpeed - self.EffTable[0][i])


    def SpeedLookup(self, Material, Temp):
        if self.SpeedTable[Material][0][-1] < Temp:
            raise Exception("Temperature out of Range for Selected Material")
        for i in range(0,len(self.SpeedTable[Material][0])-1):
            if self.SpeedTable[Material][0][i] >= Temp:
                return self.SpeedTable[Material][1][i] + (self.SpeedTable[Material][1][i+1] - self.SpeedTable[Material][1][i])/(self.SpeedTable[Material][0][i+1] - self.SpeedTable[Material][0][i]) * (Temp - self.SpeedTable[Material][0][i])
        i = len(self.SpeedTable[Material][0])-2
        return self.SpeedTable[Material][1][i] + (self.SpeedTable[Material][1][i+1] - self.SpeedTable[Material][1][i])/(self.SpeedTable[Material][0][i+1] - self.SpeedTable[Material][0][i]) * (Temp - self.SpeedTable[Material][0][i])
