#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on  2022/11/04 12:46:18
@author  Devial
'''

import tkinter as tk
from tkinter import messagebox as mb
import os.path

class GUI():
    def __init__(self, calculator, plotter, exporter):
        self.tools = {}
        self.tools["calculator"] = calculator
        self.tools["plotter"] = plotter
        self.tools["exporter"] = exporter
        self.root = 0
        self.entries = {}
        self.result = 0
        self.createMain()

    def createMain(self):
        Window =  tk.Tk()
        self.root = Window
        self.populateMain()

    def populateMain(self):
        self.root.title("Rocket Sizing Tool")
        self.root.geometry("390x400")
        label = tk.Label(self.root, text="Choose Calculation", font = ("comic sans", 30))
        label.grid(row=0,column=0,sticky="e")
        button = tk.Button(self.root, text = "Tsiolkowsky Calculator", command = self.Tsiolkowsky)
        button.config(height = 2, width=35)
        button.grid(row = 1, column = 0)
        button = tk.Button(self.root, text = "Delta V Calculator", command = self.delV)
        button.config(height = 2, width=35)
        button.grid(row = 2, column = 0)
        button = tk.Button(self.root, text = "Mass Calculator (Point)", command=self.Point)
        button.config(height = 2, width=35)
        button.grid(row = 3, column = 0)
        button = tk.Button(self.root, text = "Mass Calculator (Range)", command = self.Range)
        button.config(height = 2, width=35)
        button.grid(row = 4, column = 0)
        button = tk.Button(self.root, text = "Optimise Mass Ratio", command = self.Ratio)
        button.config(height = 2, width=35)
        button.grid(row = 5, column = 0)
        button = tk.Button(self.root, text = "Required Mission Propellant", command=self.Fuel)
        button.config(height = 2, width=35)
        button.grid(row = 6, column = 0)
    
    def Tsiolkowsky(self):
        self.root.destroy()
        Window = tk.Tk()
        self.root = Window
        self.result = tk.StringVar()
        self.result.set("")
        input_frm = tk.Frame(relief = tk.SUNKEN, borderwidth = 3)
        input_frm.pack()
        button_frm = tk.Frame()
        button_frm.pack()
        self.root.title("Tsiolkowsky Calculator")
        self.root.geometry("410x250")
        path = os.path.dirname(__file__)
        path = os.path.join(path, "back.png")
        back = tk.PhotoImage(file = path)
        button = tk.Button(input_frm, image = back, command = self.back, borderwidth=0)
        button.config(height = 40, width = 50)
        button.grid(row = 0, column = 0)
        self.addLabelColumn(0, 1, ( "Isp                  :",
                                    "Starting Mass:", 
                                    "Final Mass     :",
                                    "Result           :"),
                                    ("Arial Bold", 16), 
                                    frame = input_frm)
        entries = self.addEntryColumn(1, 1, 3, frame = input_frm)
        label = tk.Label(input_frm, textvariable = self.result, font= ("Arial Bold", 16))
        label.grid(row = 4, column = 1)
        self.entries["Isp"] = entries[0]
        self.entries["m0"] = entries[1]
        self.entries["mf"] = entries[2]
        button = tk.Button(button_frm, text = "Calculate", command = self.calcTsiolkowsky)
        button.pack(side = tk.RIGHT)
        button = tk.Button(button_frm, text = "Clear Inputs", command = self.clear)
        button.pack(side = tk.LEFT)
        self.run()

    def calcTsiolkowsky(self):
        isp = self.entries["Isp"].get() 
        m0  = self.entries["m0"].get()
        mf  = self.entries["mf"].get()
        isp = float(isp)
        m0 = float(m0)
        mf = float(mf)
        result = self.tools["calculator"].Tsiolkowsky(isp, m0, mf)
        
        self.result.set(str(round(result,7)) + " m/s")

    def clear(self):
        for entry in self.entries:
            self.entries[entry].delete(0, tk.END)
        self.result.set("")

    def back(self):
        self.root.destroy()
        self.createMain()
        self.entries = {}
        self.run()

    def addLabelColumn(self, _column, startrow, labels, _font, frame = None):
        _row = startrow
        if frame == None:
            frame = self.root
        for lbl in labels:
            label = tk.Label(frame, text = lbl, font = _font)
            label.grid(row = _row, column = _column)
            _row += 1

    def addEntryColumn(self, _column, startrow, num, frame = None):
        _row = startrow
        entries = []
        if frame == None:
            frame = self.root
        for i in range(0, num):
            entry = tk.Entry(frame, width = 50)
            entry.grid(row = _row, column = _column)
            entries.append(entry)
            _row += 1
        return entries

    def delV(self):
        self.placeholder()
    def Point(self):
        self.placeholder()
    def Range(self):
        self.placeholder()
    def Fuel(self):
        self.placeholder()
    def Ratio(self):
        self.placeholder()
    def placeholder(self):
        msg_box = mb.showerror(title="Error", message="Function not yet Implemented")
    def run(self):
        self.root.mainloop()
    
