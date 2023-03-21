#!/usr/bin/env python
'''
# -*- coding: utf-8 -*-
Created on  2023/03/21 09:23:69
@author  Devial
'''

import math
import numpy as np


class LookUpTable():

    def __init__(self):
        self.genTables()

    def genTables(self):
        self.EffTable = np.array([[0.8, 0.85, 0.95, 1.05, 1.15, 1.3, 1.55, 1.7, 2, 2.1, 2.5, 3, 4, 5, 6, 7, 8, 9, 10],[0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,  0.77, 0.79, 0.8, 0.825, 0.833, 0.85, 0.85, 0.845, 0.837, 0.83, 0.825, 0.816]])
        self.SpeedTable = np.array([[0,0],[0,0]])
    def returnTable(self, table):
        if table == 0:
            return self.EffTable
        else:
            return self.SpeedTable

    def EffLookup(self, StSpecSpeed):
        if StSpecSpeed > 10 or StSpecSpeed < 0.8:
            raise Exception("Stage Specific Speed Out Of Range")
        if StSpecSpeed > 9:
            i = 17
            return self.EffTable[1][i] + (self.EffTable[1][i+1] - self.EffTable[1][i])/(self.EffTable[0][i+1] - self.EffTable[0][i]) * (StSpecSpeed - self.EffTable[0][i])
        else:
            for i in range(0,19):
                if self.EffTable[0][i] >= StSpecSpeed:
                    return self.EffTable[1][i] + (self.EffTable[1][i+1] - self.EffTable[1][i])/(self.EffTable[0][i+1] - self.EffTable[0][i]) * (StSpecSpeed - self.EffTable[0][i])