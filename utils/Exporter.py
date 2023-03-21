# -*- coding: utf-8 -*-

"""
# -*- coding: utf-8 -*-
Created on  2022/10/26 09:23:69
@author  Devial
"""

import os.path
import numpy as np
from datetime import datetime


class Exporter(object):
    def __init__(self, DataDir):
        self.DataDir = DataDir

    def ExportData(self, data, fType, fName = None, fDir = None):
        if not os.path.isdir(self.DataDir):
            return -1
        if fName == None:
            fName = "data" + datetime.now().strftime("%Y%m%d-%H_%M_%S")
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
        if fDir is None:
            fName = os.path.join(self.DataDir, fName)
        else:
            fName = os.path.join(fDir, fName)
        fName_Full = fName + fType
        i = 1
        while os.path.isfile(fName_Full):
            fName += "(" + str(i) + ")"
            fName_Full = fName + fType
            i += 1
        np.savetxt(fName_Full, data, delimiter=",", fmt = "%.4f")
        return 0

    def AppendData(self, data, fName, fDir = None):

        if fDir is None:
            fName = os.path.join(self.DataDir, fName)
        else:
            fName = os.path.join(fDir, fName)
        with open(fName,'a') as fd:
            fd.write(data)

    def setDataDir(self, DataDir):
        self.DataDir = DataDir
        if not os.path.isdir(self.DataDir):
            return -1
        
    def verifyDataDir(self):
        if not os.path.isdir(self.DataDir):
            return -1
        else:
            return 0
