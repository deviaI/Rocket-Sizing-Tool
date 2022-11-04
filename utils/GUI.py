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
        self.data = {}
        self.result = 0
        self.activeWindow = 0
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
        self.root.title("Tsiolkowsky Calculator")
        self.root.geometry("410x250")
        input_frm, button_frm = self.basicWindowSetup()
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
        button = tk.Button(button_frm, text = "Calculate", command = self.Tsiolkowsky_Calc)
        button.pack(side = tk.RIGHT)
        button = tk.Button(button_frm, text = "Clear Inputs", command = self.clear)
        button.pack(side = tk.LEFT)
        self.run()

    def Tsiolkowsky_Calc(self):
        isp = self.entries["Isp"].get() 
        m0  = self.entries["m0"].get()
        mf  = self.entries["mf"].get()
        if isp == "" or m0 == "" or mf == "":
            self.ErrorMsg("Missing Input Arguments")
            return -1
        isp = float(isp)
        m0 = float(m0)
        mf = float(mf)
        result = self.tools["calculator"].Tsiolkowsky(isp, m0, mf)
        
        self.result.set(str(round(result,7)) + " m/s")

    def delV(self):
        """
        inputs n, m, m_pl, isp
        """
        self.root.destroy()
        Window = tk.Tk()
        self.root = Window
        self.result = tk.StringVar()
        self.result.set("")
        input_frm, button_frm = self.basicWindowSetup(parent = Window)
        path = os.path.dirname(__file__)
        path = os.path.join(path, "back.png")
        back = tk.PhotoImage(file = path)
        button = tk.Button(input_frm, image = back, command = self.back, borderwidth=0)
        button.config(height = 40, width = 50)
        button.grid(row = 0, column = 0)
        self.root.title("Delta V Calculator")
        self.root.geometry("410x300")
        self.addLabelColumn(0, 1, ( "Num. of Stages:",
                                    "Total Launch Mass:", "(incl. Payload)",
                                    "Payload Mass:", 
                                    "First Engine Isp:",
                                    "Result           :"),
                                    ("Arial Bold", 16), 
                                    frame = input_frm)
        entries = self.addEntryColumn(1, 1, 2, frame = input_frm)
        label = tk.Label(input_frm, textvariable = self.result, font= ("Arial Bold", 16))
        label.grid(row = 6, column = 1)
        self.entries["num Stages"] = entries[0]
        self.entries["m0"] = entries[1]
        entries = self.addEntryColumn(1, 4, 2, frame = input_frm)
        self.entries["m_pl"] = entries[0]
        self.entries["Isp"] = (entries[1])
        button = tk.Button(button_frm, text = "Add Fuel Masses", command = self.delV_AddFuel)
        button.pack()
        button = tk.Button(button_frm, text = "Add Isps for later Stages", command = self.delV_AddIsp)
        button.pack()
        button = tk.Button(button_frm, text = "Calculate", command = self.delV_Calc)
        button.pack(side = tk.RIGHT)
        button = tk.Button(button_frm, text = "Clear Inputs", command = self.clear)
        button.pack(side = tk.LEFT)
        self.run()

    def delV_Calc(self):
        """
        inputs n, m, m_pl, isp
        """
        n = self.entries["num Stages"].get()
        m0 = self.entries["m0"].get()
        mpl = self.entries["m_pl"].get()
        if n == "" or m0 == "" or mpl == "" or self.entries["Isp"].get() == "":
            self.ErrorMsg("Missing Input Arguments")
            return -1
        n = int(n)
        m0 = float(m0)
        mpl = float(mpl)
        isp = []
        isp.append(float(self.entries["Isp"].get()))
        try:
            for val in self.data["Isps"]:
                isp.append(val)
        except KeyError:
            isp = isp[0]
        result = self.tools["calculator"].calcDelV(n, m0,  mpl, isp)
        self.result.set(str(round(result, 4)) + " m/s" + "\t")
        

    def delV_AddIsp(self):
        try:
            num_stages = self.entries["num Stages"].get()
        except KeyError:
            self.entries = self.data["temp"]
            num_stages = self.entries["num Stages"].get()
        try:
            num_stages = int(num_stages)
        except:
            num_stages = -1
        if  num_stages <=1:
            self.ErrorMsg("Adding further ISPs is only possible if the Number of Stages is at least 2")
            return -1
        Window = tk.Toplevel(self.root)
        self.activeWindow = Window
        window_height = 40 + num_stages*30
        size = "200x" + str(int(window_height))
        Window.geometry(size)
        input_frm, button_frm = self.basicWindowSetup(parent = Window)
        labels = []
        for i in range(0,num_stages-1):
            labels.append("ISP " + str(int(i+2)) + ": ")
        self.addLabelColumn(_column = 0,startrow = 1, labels = labels, frame = input_frm, _font = ("Arial Bold", 16))
        entries = self.addEntryColumn(1, 1, num_stages-1, frame = input_frm)
        self.data["temp"] = self.entries
        self.entries = {}
        for i in range(0, num_stages-1):
            key = "isp"+str(int(i+2))
            self.entries[key] = entries[i]
        button = tk.Button(button_frm, text = "Submit", command = self.delV_AddIsp_Submit)
        button.pack(side = tk.RIGHT)
        button = tk.Button(button_frm, text = "Clear Inputs", command = self.clear)
        button.pack(side = tk.LEFT)


    def delV_AddIsp_Submit(self):
        self.data["Isps"] = []
        for entry in self.entries:
            data = self.entries[entry].get()
            if data != "":
                try:
                    self.data["Isps"].append(float(self.entries[entry].get()))
                except:
                    self.ErrorMsg("Invalid ISP in at least one Field")
                    self.activeWindow.destroy()
                    self.entries = self.data["temp"]
                    self.delV_AddIsp()
                    return -1
        print("Isp List:")
        print(self.data["Isps"])
        self.entries = self.data["temp"]
        self.activeWindow.destroy()
        


    def delV_AddFuel(self):
        self.placeholder()

    def ErrorMsg(self, text):
        msg_box = mb.showerror(title="Error", message=text)

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


    def basicWindowSetup(self, parent = None):
        if parent == None:
            parent = self.root
        input_frm = tk.Frame(parent, relief = tk.SUNKEN, borderwidth = 3)
        input_frm.pack()
        button_frm = tk.Frame(parent)
        button_frm.pack()
        return input_frm, button_frm


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
    
