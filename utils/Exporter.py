# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 15:17:69 2022

@author: Devial
"""

import os.path
import numpy as np


class Exporter(object):
    def __init__(self, DataDir):
        self.DataDir = DataDir

    def ExportData(self, data, fName, fType, fDir = None):
  
        if fType[0] != ".":
            fType = "." + fType
        if fType != ".csv" and fType != ".txt" and fType != ".xlsx":
            print("Warning: Unknwon FileType")
            print("Do you wish to proceed ? [y/n]")
            InP = input()
            if InP == "n" or InP == "N":
                print("Enter New Filetype")
                self.ExportData(data, fName, input(), fDir)
                return 0
            elif InP != "y" and InP != "Y":
                print("unrecognised Input")
                self.ExportData(data, fName, fType, fDir)
                return 0
        fName += fType
        if fDir is None:
            fName = os.path.join(self.DataDir, fName)
        else:
            fName = os.path.join(fDir, fName)
        np.savetxt(fName, data, delimiter=",")

    def AppendData(self, data, fName, fDir = None):
        if fDir is None:
            fName = os.path.join(self.DataDir, fName)
        else:
            fName = os.path.join(fDir, fName)
        with open(fName,'a') as fd:
            fd.write(data)
    
