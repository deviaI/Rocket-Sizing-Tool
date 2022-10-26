# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 15:17:69 2022

@author: Devial
"""

import os.path
import numpy as np


class Eporter(object):
    """
    NOT VERIFIED
    """
    def __init__(self, DataDir):
        self.DataDir = DataDir

    def ExportData(self, data, fName, fType, fDir = None):
        fName += fType
        if fType[0] != ".":
            fType = "." + fType
        if fType != ".csv" and fType != ".txt" and fType != ".xlsx":
            print("Warning: Unknwon FileType")
            print("Do you wish to proceed ? [y/n]")
            if input() == "n":
                print("Enter New Filetype")
                self.ExportData(data, fName, input(), fDir)
                return 0
        if fDir is None:
            fName = os.path.join(self.DataDir, fName)
        else:
            fName = os.path.join(fDir, fName)
        np.savetxt(fName, data, delimiter=",")

