#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on  2022/11/04 12:46:18
@author  Devial
'''

import tkinter as tk
from tkinter import messagebox as mb
import os.path
from tkinter import filedialog as fd

class GUI():
    def __init__(self, calculator, plotter, exporter, exeDir):
        self.tools = {}
        self.tools["calculator"] = calculator
        self.tools["plotter"] = plotter
        self.tools["exporter"] = exporter
        if exeDir == None:
            self.FileDir = os.path.join(os.path.dirname(__file__), "files")
        else:
            self.FileDir = os.path.join(exeDir, "files")
        if not os.path.isdir(self.FileDir):
            self.ErrorMsg("Fatal - Missing Files")
            return -1
        with open(os.path.join(self.FileDir, "HelpMessages.txt")) as f:
            self.help_msgs = f.readlines()
            f.close()
        try:
            with open(os.path.join(self.FileDir, ".datadir")) as f:
                self.DataDir = f.readlines()[0]
        except:
            self.DataDir = ""
        self.root = 0
        self.entries = {}
        self.data = {}
        self.result = 0
        self.img1 = None
        self.activeWindow = 0
        self.createMain()

    def createMain(self):
        Window =  tk.Tk()
        self.root = Window
        self.root.title("Rocket Sizing Tool")
        self.help = tk.PhotoImage(file = os.path.join(self.FileDir, "help.png"))
        self.root.geometry("390x390")
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
        button = tk.Button(self.root, text = "Ascent Calculation", command=self.Ascent)
        button.config(height = 2, width=35)
        button.grid(row = 7, column = 0)
        button = tk.Button(self.root, text = "Booster Calculations", command=self.Booster_Parent)
        button.config(height = 2, width=35)
        button.grid(row = 8, column = 0)
        button = tk.Button(self.root, image = self.help, command = self.help_msg)
        button.image = self.help
        button.config(width = 20, height = 20)
        button.grid(row = 0, column = 1)
        self.data["title"] = "Rocket Sizing Tool"
        self.data["message"] = "Tool for basic sizing of Rockets, using the Ziolkowsky Equation. \n To begin a new Calculation select one of the Options"
    
    def Booster_Parent(self):
        self.root.destroy()
        Window =  tk.Tk()
        self.root = Window
        self.root.title("Booster Calculations")
        self.help = tk.PhotoImage(file = os.path.join(self.FileDir, "help.png"))
        self.root.geometry("390x300")
        path = os.path.join(self.FileDir, "back.png")
        back = tk.PhotoImage(file = path)
        button = tk.Button(Window, image = back, command = self.back, borderwidth=0)
        button.image = back
        button.config(height = 40, width = 50)
        button.grid(row = 4, column = 0)
        path = os.path.join(self.FileDir, "help.png")
        help = tk.PhotoImage(file = path)
        button = tk.Button(Window, image = help, command = self.help_msg)
        button.image = help
        button.config(width = 20, height = 20)
        button.grid(row = 0, column = 1)
        label = tk.Label(self.root, text="Choose Calculation", font = ("comic sans", 30))
        label.grid(row=0,column=0,sticky="e")
        button = tk.Button(self.root, text = "Flexible Boosters & Core", command = self.FlexFlex)
        button.config(height = 2, width=35)
        button.grid(row = 1, column = 0)
        button = tk.Button(self.root, text = "Flexible Boosters & Fixed Core", command = self.FlexFix)
        button.config(height = 2, width=35)
        button.grid(row = 2, column = 0)
        button = tk.Button(self.root, text = "Fixed Boosters & Core", command=self.FixFix)
        button.config(height = 2, width=35)
        button.grid(row = 3, column = 0)

    def FlexFlex(self):
        self.placeholder()
    def FlexFix(self):
        self.placeholder()
    def FixFix(self):
        self.placeholder()

    def Ascent(self):
        self.root.destroy()
        Window = tk.Tk()
        self.root = Window
        self.result = tk.StringVar()
        self.result.set("")
        self.root.title("Tsiolkowsky Calculator")
        self.root.geometry("370x650")
        msg = ""
        i1 = self.help_msgs.index("6\n")
        i2 = self.help_msgs.index("7\n")
        for i in range(i1+1, i2):
            msg += self.help_msgs[i]
        input_frm, button_frm = self.basicWindowSetup(help_msg = msg)
        self.addLabelColumn(0, 1, ( "C1 [-] ",
                                    "Thrust [N]", 
                                    "Latitude: [°]",
                                    "C2 [-]",
                                    "Mass [kg]",
                                    "CD [-]",
                                    "A [m^2]",
                                    "Mass Flow [kg/s]",
                                    "",
                                    "Transition Mode",
                                    "dt [s]",
                                    "Cutoff Alt. [m]",
                                    "Cutoff Mass [kg]",
                                    "TVC Rate [°/s]",
                                    "Throttle Rate [1/s]",
                                    "G Limit [m/s^2]",
                                    ),
                                    ("Arial Bold", 16), 
                                    frame = input_frm)
        entries = self.addEntryColumn(1, 1, 9, frame = input_frm)
        self.entries["C1"] = entries[0]
        self.entries["Thrust"] = entries[1]
        self.entries["Beta"] = entries[2]
        self.entries["C2"] = entries[3]
        self.entries["Mass"] = entries[4]
        self.entries["CD"] = entries[5]
        self.entries["A"] = entries[6]
        self.entries["propburn"] = entries[7]
        self.entries["Grav_Trans"] = entries[8]
        entries = self.addEntryColumn(1, 11, 6, frame = input_frm)
        self.entries["dt"] = entries[0]
        self.entries["cutoff_h"] = entries[1]
        self.entries["cutoff_m"] = entries[2]
        self.entries["steer_rate"] = entries[3]
        self.entries["throt_rate"] = entries[4]
        self.entries["a_max"] = entries[5]
        self.entries["Beta"].insert(0, "0")
        self.entries["dt"].insert(0, "0.001")
        self.entries["steer_rate"].insert(0, "1")
        self.entries["throt_rate"].insert(0, "0.5")
        self.entries["a_max"].insert(0, "100")
        self.data["ButtonLbl"] = tk.StringVar()
        self.data["ButtonLbl"].set("Altitude")
        self.data["TransitionMode"] = -1
        button = tk.Button(input_frm, textvariable=  self.data["ButtonLbl"], command = self.Ascent_SwitchGT)
        button.grid(row = 10, column = 1)
        label = tk.Label(input_frm, textvariable = self.result, font= ("Arial Bold", 16))
        label.grid(row = 15, column = 1)
        TransLabel = tk.StringVar()
        TransLabel.set("Grav. Transition [m]")
        self.data["TransLabel"] = TransLabel
        label = tk.Label(input_frm, textvariable = self.data["TransLabel"], font= ("Arial Bold", 16))
        label.grid(row = 9, column = 0)
        button = tk.Button(button_frm, text = "Add List of Thrusts", command = self.Ascent_AddThrust)
        button.pack()
        button = tk.Button(button_frm, text = "Add List of CDs", command = self.Ascent_AddCDs)
        button.pack()
        button = tk.Button(button_frm, text = "Add List of Areas", command = self.Ascent_AddAreas)
        button.pack()
        button = tk.Button(button_frm, text = "Add List of Mass Flows", command = self.Ascent_AddMassFlows)
        button.pack()
        button = tk.Button(button_frm, text = "Preview Ascent Profile", command = self.Ascent_Preview_profile)
        button.pack()
        button = tk.Button(button_frm, text = "Calculate", command = self.Ascent_Calc)
        button.pack(side = tk.RIGHT)
        button = tk.Button(button_frm, text = "Clear Inputs", command = self.clear)
        button.pack(side = tk.LEFT)
        self.run()
    
    def Ascent_Preview_profile(self):
        C1 = self.entries["C1"].get()
        C2 = self.entries["C2"].get()
        Grav_Trans = self.entries["Grav_Trans"].get()
        h_cutoff = self.entries["cutoff_h"].get()
        if self.data["TransitionMode"] == 1:
            self.ErrorMsg("Preview only available in Altitude Transition Mode")
            return -1
        if C1 == "" or C2 == "" or Grav_Trans == "":
            self.ErrorMsg("Profile Preview requires C1, C2 and Transition Altitude")
            return -1
        C1 = float(C1)
        C2 = int(C2)
        Grav_Trans = int(Grav_Trans)
        if h_cutoff != "":
            h_cutoff = int(h_cutoff)
        else:
            h_cutoff = 200000
        h, x = self.tools["calculator"].gen_ascent_path_preview(C1, C2, Grav_Trans, h_cutoff)
        x = [val*1e-3 for val in x]
        karman_line = [100000 for k in x]
        self.tools["plotter"].plot2D(x, h, dataY_List = [karman_line], xlim = [-max(x)*0.1, max(x)], xLab = "Downrange [km]", yLab = "Altitude [m]", data_Labels = ["Ascent Profile" ,"Karman Line"]) 

    def Ascent_SwitchGT(self):
        self.data["TransitionMode"] *= -1
        if self.data["TransitionMode"] == -1:
            self.data["ButtonLbl"].set("Altitude")
            self.data["TransLabel"].set("Grav. Transition [m]")
        elif self.data["TransitionMode"] == 1:
            self.data["ButtonLbl"].set("Velocity")
            self.data["TransLabel"].set("Grav. Transition [m/s]")

    def Ascent_AddThrust(self):
        self.placeholder()
    def Ascent_AddCDs(self):
        self.placeholder()
    def Ascent_AddAreas(self):
        self.placeholder()
    def Ascent_AddMassFlows(self):
        self.placeholder()
    def Ascent_Calc(self):
        self.placeholder()
    
    def Tsiolkowsky(self):
        self.root.destroy()
        Window = tk.Tk()
        self.root = Window
        self.result = tk.StringVar()
        self.result.set("")
        self.root.title("Tsiolkowsky Calculator")
        self.root.geometry("370x230")
        msg = ""
        i1 = self.help_msgs.index("0\n")
        i2 = self.help_msgs.index("1\n")
        for i in range(i1+1, i2):
            msg += self.help_msgs[i]
        input_frm, button_frm = self.basicWindowSetup(help_msg = msg)
        self.addLabelColumn(0, 1, ( "Isp                  :",
                                    "Starting Mass:", 
                                    "Final Mass     :",
                                    "Result            :"),
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
        try:
            result = self.tools["calculator"].Tsiolkowsky(isp, m0, mf)
        except Exception as e:
            self.ErrorMsg(str(e))
        self.result.set(str(round(result,7)) + " m/s")

    def delV(self):
        """
        inputs n, m, m_pl, isp
        """
        self.root.destroy()
        Window = tk.Tk()
        self.root = Window
        self.result = tk.StringVar()
        msg = ""
        i1 = self.help_msgs.index("1\n")
        i2 = self.help_msgs.index("2\n")
        for i in range(i1+1, i2):
            msg += self.help_msgs[i]
        input_frm, button_frm = self.basicWindowSetup(parent = Window, help_msg = msg)
        self.root.title("Delta V Calculator")
        self.root.geometry("390x330")
        self.addLabelColumn(0, 1, ( "Num. of Stages     :",
                                    "Launch Mass :", "(excl. Payload)",
                                    "Payload Mass  :", 
                                    "First Engine Isp:",
                                    "Result                  :"),
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
        button = tk.Button(button_frm, text = "Add later Stage Isps", command = self.AddIsp)
        button.pack()
        button = tk.Button(button_frm, text = "Specify List of Stage Masses", command = self.delV_AddStageMasses)
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
        mpl = float(mpl)
        isp = []
        m_f = []
        m_s = []
        isp.append(float(self.entries["Isp"].get()))
        if "Isps" in self.data:
            for val in self.data["Isps"]:
                isp.append(val)
        else:
            isp = isp[0]
        if "fuels" in self.data:
            for val in self.data["fuels"]:
                m_f.append(val)
        if "masses" in self.data:
            for val in self.data["masses"]:
                m_s.append(val)
        else:
            m_s = float(m0)
        try:
            if len(m_f) == 0:
                result = self.tools["calculator"].calcDelV(n, m_s,  mpl, isp)
            else:
                result = self.tools["calculator"].calcDelV(n, m_s, mpl, isp, m_f = m_f)
        except Exception as e:
            self.ErrorMsg(str(e))
        self.result.set(str(round(result, 4)) + "m/s" + "\t")
        
    def delV_AddFuel(self):
        try:
            num_stages = self.entries["num Stages"].get()
        except KeyError:
            self.entries = self.data["temp"]
            num_stages = self.entries["num Stages"].get()
        try:
            num_stages = int(num_stages)
        except:
            num_stages = -1
        if  num_stages <=0:
            self.ErrorMsg("Number of Stages must be entered first")
            return -1
        self.ListInputWindow(Field_Name = "Fuel, Stage ", num_inputs = num_stages, key_Name = "fuels")
    
    def delV_AddStageMasses(self):
        try:
            num_stages = self.entries["num Stages"].get()
        except KeyError:
            self.entries = self.data["temp"]
            num_stages = self.entries["num Stages"].get()
        try:
            num_stages = int(num_stages)
        except:
            num_stages = -1
        if  num_stages <=0:
            self.ErrorMsg("Number of Stages must be entered first")
            return -1
        self.ListInputWindow(Field_Name= "Mass of Stage ", key_Name="masses", num_inputs = num_stages)

    def Point(self):
        """
        inputs n, isp, m_pl, mu = 0.12, delv = 9000, limit = 1e6
        opt: Size Fact.
        """
        self.root.destroy()
        Window = tk.Tk()
        self.root = Window
        self.result = [tk.StringVar(), tk.StringVar()]
        self.result[0].set("")
        self.result[1].set("")
        msg = ""
        i1 = self.help_msgs.index("2\n")
        i2 = self.help_msgs.index("3\n")
        for i in range(i1+1, i2):
            msg += self.help_msgs[i]
        input_frm, button_frm = self.basicWindowSetup(parent = Window, help_msg = msg)
        self.root.title("Mass Calculator (Point)")
        self.root.geometry("390x380")
        self.addLabelColumn(0, 1, ["Number of Stages", 
                                    "First Engine Isp",
                                    "Payload Mass",
                                    "Structural Factor",
                                    "target delta v",
                                    "divergence limit",
                                    "Rocket Mass",
                                    "Stage Rel. Mass",
                                    "Req. tot. Prop."], _font = ("Arial Bold", 16), frame = input_frm)
        label = tk.Label(input_frm, textvariable = self.result[0], font= ("Arial Bold", 16))
        label.grid(row = 7, column = 1)
        label = tk.Label(input_frm, textvariable = self.result[1], font= ("Arial Bold", 16))
        label.grid(row = 9, column = 1)
        entries = self.addEntryColumn(1,1, 6, frame = input_frm)
        entry = tk.Entry(input_frm, width = 30)
        entry.grid(row = 8, column = 1)
        self.entries["num Stages"] = entries[0]
        self.entries["Isp"] = entries[1]
        self.entries["m_pl"] = entries[2]
        self.entries["mu"] = entries[3]
        self.entries["Delta V"] = entries[4]
        self.entries["Limit"] = entries[5]
        self.entries["Size Fac"] = entry
        self.entries["Delta V"].insert(0, "9000")
        self.entries["mu"].insert(0, "0.12")
        self.entries["Limit"].insert(0, "1000000")
        self.entries["Size Fac"].insert(0, "Optional")
        button = tk.Button(button_frm, text = "Add Isp List", command = self.AddIsp)
        button.pack()
        button = tk.Button(button_frm, text = "Calculate", command = self.Point_Calc)
        button.pack(side = tk.RIGHT)
        button = tk.Button(button_frm, text = "Reset Inputs", command = self.Point)
        button.pack(side = tk.LEFT)
        self.run()

    def AddIsp(self):
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
        self.ListInputWindow(Field_Name= "Isp, Stage ", key_Name="Isps", num_inputs = num_stages-1, startoffset=2)

    def Point_Calc(self):
        n = self.entries["num Stages"].get()
        isp = self.entries["Isp"].get()
        mpl = self.entries["m_pl"].get()
        mu = self.entries["mu"].get()
        delv = self.entries["Delta V"].get()
        limit = self.entries["Limit"].get()
        mu = self.entries["mu"].get()
        size_fac = self.entries["Size Fac"].get()
        if n == "" or isp == "" or mpl == "" or mu == "" or delv == "" or limit == "":
            self.ErrorMsg("Missing Input Arguments")
            return -1
        if size_fac == "" or size_fac == "Optional":
            size_fac = -1
        try:
            n = int(n)
            isp = float(isp)
            mpl = float(mpl)
            mu = float(mu)
            delv = float(delv)
            limit = float(limit)
            mu = float(mu)
            size_fac = float(size_fac)
        except:
            self.ErrorMsg("Invalid Input")
            return -1
        if "Isps" in self.data:
            isp = [isp]
            for val in self.data["Isps"]:
                isp.append(val)
        try:
            if size_fac != -1:
                result = self.tools["calculator"].calcPoint(n, isp, mpl, mu, delv, limit, size_fac = size_fac)
            else:
                result = self.tools["calculator"].calcPoint(n, isp, mpl, mu, delv, limit)
        except Exception as e:
            self.ErrorMsg(str(e))
        if result[0] > 1e50:
            for result in self.result:
                result.set("Calculation Diverged")
        self.result[0].set(str(round(result[0], 3)) + "kg")
        self.entries["Size Fac"].delete(0, tk.END)
        self.entries["Size Fac"].insert(0, str(result[1]*100)+ "%")
        self.result[1].set(str(round(sum(result[2]), 3)) + "kg")

    def Range(self):
        """
        inputs n, isp, m_pl, mu = 0.12, delv = 9000, limit = 1e6
        opt: Size Fact.
        """
        self.root.destroy()
        Window = tk.Tk()
        self.root = Window
        msg = ""
        i1 = self.help_msgs.index("3\n")
        i2 = self.help_msgs.index("4\n")
        for i in range(i1+1, i2):
            msg += self.help_msgs[i]
        input_frm, button_frm = self.basicWindowSetup(parent = Window, help_msg = msg)
        self.data["input Frame"] = input_frm
        self.root.title("Mass Calculator (Range)")
        self.root.geometry("390x300")
        self.addLabelColumn(0, 1, ["Number of Stages", 
                                    "Payload Mass",
                                    "target delta v",
                                    "divergence limit",
                                    "Stage Rel. Mass",
                                    "Number of Steps"], _font = ("Arial Bold", 16), frame = input_frm)
        entries = self.addEntryColumn(1,1, 6, frame = input_frm)
        self.entries["num Stages"] = entries[0]
        self.entries["m_pl"] = entries[1]
        self.entries["Delta V"] = entries[2]
        self.entries["Limit"] = entries[3]
        self.entries["Size Fac"] = entries[4]
        self.entries["Num Steps"] = entries[5]

        self.entries["Delta V"].insert(0, "9000")
        self.entries["Limit"].insert(0, "1000000")
        self.entries["Size Fac"].insert(0, "Optional")
        self.entries["Num Steps"].insert(0, "100")
        top_frame = tk.Frame(button_frm, borderwidth = 3)
        top_frame.pack()
        bottom_frame = tk.Frame(button_frm, borderwidth = 3)
        bottom_frame.pack()
        button = tk.Button(top_frame, text = "Fix Isp", command = self.Range_FixIsp, bg = "#81fcf4")
        button.pack(side = tk.LEFT)
        button = tk.Button(top_frame, text = "Fix Mu", command = self.Range_FixMu, bg = "#81fcf4")
        button.pack(side = tk.RIGHT)
        button = tk.Button(bottom_frame, text = "Calculate", command = self.Range_Calc)
        button.pack(side = tk.RIGHT)
        button = tk.Button(bottom_frame, text = "Reset Inputs", command = self.Range)
        button.pack(side = tk.LEFT)
        self.run()
    
    def Range_FixMu(self):
        self.root.geometry("390x360")
        label = tk.Label(self.data["input Frame"], text = "Mu", font = ("Arial Bold", 16))
        label.grid(row = 7, column = 0)
        label = tk.Label(self.data["input Frame"], text = "Isp Start, Isp End", font = ("Arial Bold", 16))
        label.grid(row = 8, column = 0)
        range_frame = tk.Frame(self.data["input Frame"], borderwidth = 3)
        range_frame.grid(row = 8, column = 1)
        entry = tk.Entry(self.data["input Frame"], width = 30)
        entry.grid(row = 7, column = 1)
        self.entries["mu"] = entry
        if "Isp" in self.entries:
            del self.entries["isp"]
        entry = tk.Entry(range_frame, width = 15)
        entry.pack(side = tk.LEFT)
        self.entries["Range Lower"] = entry
        entry = tk.Entry(range_frame, width = 15)
        entry.pack(side = tk.RIGHT)
        self.entries["Range Upper"] = entry
        self.data["Range Var"] = "Isp"

    def Range_FixIsp(self):
        self.root.geometry("390x360")
        label = tk.Label(self.data["input Frame"], text = "Isp", font = ("Arial Bold", 16))
        label.grid(row = 7, column = 0)
        label = tk.Label(self.data["input Frame"], text = "Mu Start, Mu End", font = ("Arial Bold", 16))
        label.grid(row = 8, column = 0)
        range_frame = tk.Frame(self.data["input Frame"], borderwidth = 3)
        range_frame.grid(row = 8, column = 1)
        entry = tk.Entry(self.data["input Frame"], width = 30)
        entry.grid(row = 7, column = 1)
        self.entries["Isp"] = entry
        if "mu" in self.entries:
            del self.entries["mu"]
        entry = tk.Entry(range_frame, width = 15)
        entry.pack(side = tk.LEFT)
        self.entries["Range Lower"] = entry
        entry = tk.Entry(range_frame, width = 15)
        entry.pack(side = tk.RIGHT)
        self.entries["Range Upper"] = entry
        self.data["Range Var"] = "Mu"
    
    def Range_Calc(self):
        try:
            mu = self.entries["mu"].get()
            isp = -1
        except:
            try:
                isp = self.entries["Isp"].get()
                mu = -1
            except:
                self.ErrorMsg("Select Fix Variable First")
                return -1
        n = self.entries["num Stages"].get()
        mpl = self.entries["m_pl"].get()
        delv = self.entries["Delta V"].get()
        limit = self.entries["Limit"].get()
        size_fac = self.entries["Size Fac"].get()
        num_steps = self.entries["Num Steps"].get()
        try: 
            lower = self.entries["Range Lower"].get()
            upper = self.entries["Range Upper"].get()
        except KeyError:
            self.ErrorMsg("Select Fix Variable First")
            return -1
        if n == "" or isp == "" or mpl == "" or mu == "" or delv == "" or limit == "" or lower == "" or upper == "" or num_steps == "":
            self.ErrorMsg("Missing Input Arguments")
            print(n)
            print(isp)
            print(mpl)
            print(mu)
            print(delv)
            print(limit)
            print(lower)
            print(upper)
            print(num_steps)
            return -1
        if size_fac == "" or size_fac == "Optional":
            size_fac = -1
        try:
            n = int(n)
            isp = float(isp)
            mpl = float(mpl)
            mu = float(mu)
            delv = float(delv)
            limit = float(limit)
            mu = float(mu)
            size_fac = float(size_fac)
            lower = float(lower)
            upper = float(upper)
            num_steps = int(num_steps)
        except:
            self.ErrorMsg("Invalid Input")
            return -1
        if self.data["Range Var"] == "Isp":
            isp = (lower, upper)
        elif self.data["Range Var"] == "Mu":
            mu = (lower, upper)
        try:
            if size_fac != -1:
                result = self.tools["calculator"].calcRange(n, self.data["Range Var"], isp, mpl, mu, delv, limit, num_steps, size_fac = size_fac)
            else:
                result = self.tools["calculator"].calcRange(n, self.data["Range Var"], isp, mpl, mu, delv, limit, num_steps)
        except Exception as e:
            self.ErrorMsg(str(e))
        self.tools["plotter"].plot2D(dataX = result[1], dataY = result[0], yLab = "Mass", xLab = self.data["Range Var"], show = 1, savefile = 0)
        self.tools["plotter"].plot2D(dataX = result[1], dataY = result[2], yLab = "Optimal Rel. Stage Sizing", xLab = self.data["Range Var"], show =1, savefile = 0)
        answer = mb.askquestion("Mass Calculator (Range)", "Save generated data and plots ?")
        if answer == "yes":
            if os.path.isdir(self.DataDir):
                datadir = fd.askdirectory(initialdir = self.DataDir)
            else:
                datadir = fd.askdirectory()
            try:
                with open(os.path.join(self.FileDir, ".datadir"), "w") as f:
                    f.write(str(datadir))
                    f.close()
            except:
                pass
            self.tools["exporter"].setDataDir(datadir)
            self.tools["plotter"].setDataDir(datadir)
            ret = self.tools["plotter"].plot2D(dataX = result[1], dataY = result[0], yLab = "Mass", xLab = self.data["Range Var"], show = 0, savefile = 1)
            ret += self.tools["plotter"].plot2D(dataX = result[1], dataY = result[2], yLab = "Optimal Rel. Stage Sizing", xLab = self.data["Range Var"], show =0, savefile = 1)
            ret += self.tools["exporter"].ExportData(data = result, fType = ".csv")
            if ret != 0:
                self.ErrorMsg("Failed to save Data in" + str(datadir))
            self.Range()
        else:
            self.Range()

    def Ratio(self):
        """
        Inputs n, m_0, isp, m_pl, mu
        """
        self.root.destroy()
        Window = tk.Tk()
        self.root = Window
        self.result = tk.StringVar()
        msg = ""
        i1 = self.help_msgs.index("4\n")
        i2 = self.help_msgs.index("5\n")
        for i in range(i1+1, i2):
            msg += self.help_msgs[i]
        input_frm, button_frm = self.basicWindowSetup(parent = Window, help_msg = msg)
        self.data["input Frame"] = input_frm
        self.root.title("Mass Calculator (Range)")
        self.root.geometry("390x300")
        self.addLabelColumn(0, 1, ["Number of Stages", 
                                    "Launch Mass",
                                    "Isp, Stage 1",
                                    "Payload Mass",
                                    "Structure Factor",
                                    "Result"
                                    ], _font = ("Arial Bold", 16), frame = input_frm)
        entries = self.addEntryColumn(1,1, 5, frame = input_frm)
        self.entries["num Stages"] = entries[0]
        self.entries["m"] = entries[1]
        self.entries["Isp"] = entries[2]
        self.entries["m_pl"] = entries[3]
        self.entries["Mu"] = entries[4]
        self.entries["Mu"].insert(0, "0.12")
        label = tk.Label(input_frm, textvariable = self.result,font = ("Arial Bold", 16))
        label.grid(row = 6, column = 1)
        button = tk.Button(button_frm, text = "Add later Stage ISPs", command = self.AddIsp)
        button.pack()
        button = tk.Button(button_frm, text = "Calculate", command = self.Ratio_Calc)
        button.pack(side = tk.RIGHT)
        button = tk.Button(button_frm, text = "Reset Inputs", command = self.Ratio)
        button.pack(side = tk.LEFT)
    
    def Ratio_Calc(self):
        n = self.entries["num Stages"].get()
        mpl = self.entries["m_pl"].get()
        Isp = self.entries["Isp"].get()
        m = self.entries["m"].get()
        mu = self.entries["Mu"].get()
        if n == "" or mpl == "" or Isp == "" or m == "" or mu == "":
            self.ErrorMsg("Missing Input Arguments")
            return -1
        try:
            n = int(n)
            mpl = float(mpl)
            Isp = [float(Isp)]
            m = float(m)
            mu = float(mu)
        except:
            self.ErrorMsg("Invalid Input Argument")
        if "Isps" in self.data:
            for val in self.data["Isps"]:
                Isp.append(val)
        try:
            result = self.tools["calculator"].optimiseMassRatio(n, m, Isp, mpl, mu)["Optimal Stage Sizing Factor"]
        except Exception as e:
            self.ErrorMsg(str(e))
        self.result.set(str(round(result*100,0)) + "%" )

    def Fuel(self):
        """
        n, isp, m_pl, mu = 0.12, delv = 9000, limit = 1e6
        size_fac: relative stage to stage mass
                        If not given, it is calculated to be optimal
            mix_rat: Mass mixture ratio O/F 
            dens: density of propelant [dens Oxidzer, dens Fuel]
            fueltype: fuel type. If handed, mix_rat and dens are taken from typical values
                    Recognised Values: HydroLox, KeroLox, MethaLox
                    // if both fueltype and mix_rat are handed, the specified mixture ratio will be used
        """
        self.root.destroy()
        Window = tk.Tk()
        self.root = Window
        self.result = tk.StringVar()
        msg = ""
        i1 = self.help_msgs.index("5\n")
        i2 = self.help_msgs.index("6\n")
        for i in range(i1+1, i2):
            msg += self.help_msgs[i]
        input_frm, button_frm = self.basicWindowSetup(parent = Window, help_msg=msg)
        self.data["input Frame"] = input_frm
        self.root.title("Mass Calculator (Range)")
        self.root.geometry("390x500")
        self.addLabelColumn(0, 1, ["Number of Stages", 
                                    "Isp, Stage 1",
                                    "Payload Mass",
                                    "Structure Factor",
                                    "Target Delta V",
                                    "Stage Rel. Mass",
                                    "Convergence Limit",
                                    "Mixture Ratio",
                                    "Densities (Ox, Fu)",
                                    "Fueltype",
                                    "Oxidiser",
                                    "Fuel",
                                    "Total"
                                    ], _font = ("Arial Bold", 16), frame = input_frm)
        dens_frame = tk.Frame(input_frm, borderwidth= 3)
        entry = tk.Entry(dens_frame, width = 10)
        entry.grid(row = 0, column = 0)
        self.entries["Density Ox"] = entry
        entry = tk.Entry(dens_frame, width = 10)
        entry.grid(row = 0, column = 1)
        self.entries["Density Fu"] = entry
        dens_frame.grid(row = 9, column = 1)

        entries = self.addEntryColumn(1,1, 8, frame = input_frm)
        self.entries["num Stages"] = entries[0]
        self.entries["Isp"] = entries[1]
        self.entries["m_pl"] = entries[2]
        self.entries["Mu"] = entries[3]
        self.entries["Delta V"] = entries[4]
        self.entries["Size Fac"] = entries[5]
        self.entries["Limit"] = entries[6]
        self.entries["Mix Ratio"] = entries[7]



        entries = self.addEntryColumn(1, 10, 1, frame = input_frm)
        self.entries["Fueltype"] = entries[0]
        self.entries["Mu"].insert(0, "0.12")
        self.entries["Mix Ratio"].insert(0,"Optional")
        self.entries["Density Ox"].insert(0, "Optional")
        self.entries["Density Fu"].insert(0, "Optional")
        self.entries["Fueltype"].insert(0, "Optional")
        self.entries["Size Fac"].insert(0, "Optional")
        self.entries["Limit"].insert(0, "1000000")

        self.result = [tk.StringVar(), tk.StringVar(), tk.StringVar()]
        self.result[0].set("")
        self.result[1].set("")
        self.result[2].set("")
        label = tk.Label(input_frm, textvariable=self.result[0], font = ("Arial", 16))
        label.grid(row = 11, column = 1)
        label = tk.Label(input_frm, textvariable=self.result[1], font = ("Arial", 16))
        label.grid(row = 12, column = 1)
        label = tk.Label(input_frm, textvariable=self.result[2], font = ("Arial", 16))
        label.grid(row = 13, column = 1)
        button = tk.Button(button_frm, text = "Add later Stage Isps", command = self.AddIsp)
        button.pack()
        button = tk.Button(button_frm, text = "Calculate", command = self.Fuel_Calc)
        button.pack(side = tk.RIGHT)
        button = tk.Button(button_frm, text = "Reset Inputs", command = self.Fuel)
        button.pack(side = tk.LEFT)

    def Fuel_Calc(self):

        n = self.entries["num Stages"].get()
        isp = self.entries["Isp"].get()
        mpl = self.entries["m_pl"].get()
        mu = self.entries["Mu"].get()
        delv = self.entries["Delta V"].get()
        limit = self.entries["Limit"].get()
        size_fac = self.entries["Size Fac"].get()
        Mix_Rat = self.entries["Mix Ratio"].get()
        Dens = (self.entries["Density Ox"].get(), self.entries["Density Fu"].get())
        fueltype = self.entries["Fueltype"].get()

        if n == "" or isp == "" or mpl == "" or mu == "" or delv == "" or limit == "":
            self.ErrorMsg("Missing Input Arguments")
            return -1
        if size_fac == "" or size_fac == "Optional":
            size_fac = -1
        if Mix_Rat == "" or Mix_Rat == "Optional":
            Mix_Rat = -1
        if Dens[0] == "" or Dens[1] == "" or Dens[0] == "Optional" or Dens[1] == "Optional":
            Dens = (-1, -1)
        if fueltype == "" or fueltype == "Optional":
            fueltype = -1

        n = int(n)
        isp = [float(isp)]
        mpl = float(mpl)
        mu = float(mu)
        delv = float(delv)
        limit = float(limit)
        size_fac = float(size_fac)
        Mix_Rat = float(Mix_Rat)
        Dens = (float(Dens[0]), float(Dens[1]))
        if "Isps" in self.data:
            for val in self.data["Isps"]:
                isp.append(val)
        if len(isp) == 1:
            isp = isp[0]
        try:
            result = self.tools["calculator"].calcReqFuel(n, isp, mpl, mu, delv, limit, size_fac = size_fac, mix_rat = Mix_Rat, dens = Dens, fueltype = fueltype)
        except Exception as e:
            self.ErrorMsg(str(e))
        if result[1] != -1:
            if result[3] != -1:
                Ox_Res = str(round(result[1],2)) + "kg, " + str(round(result[3], 2)) + "l"
                Fu_Res = str(round(result[2],2)) + "kg, " + str(round(result[4], 2)) + "l"
            else:
                Ox_Res = str(round(result[1],2)) + "kg"
                Fu_Res = str(round(result[2],2)) + "kg"
        else:
            Ox_Res = ""
            Fu_Res = ""
        Tot_Res = str(round(result[0], 2)) + "kg"
        self.result[0].set(Ox_Res)
        self.result[1].set(Fu_Res)
        self.result[2].set(Tot_Res)

    def ErrorMsg(self, text):
        msg_box = mb.showerror(title="Error", message=text)

    def clear(self):
        try:
            for entry in self.entries:
                self.entries[entry].delete(0, tk.END)
        except: 
            self.entries = self.data["temp"]
            for entry in self.entries:
                self.entries[entry].delete(0, tk.END)
        try:
            for result in self.result:
                result.set("")
        except:
            self.result.set("")
        temp = None
        W_Dat = None
        if "temp" in self.data:
            temp = self.data["temp"]
        if "Window Data" in self.data:
            W_Dat = self.data["Window Data"]
        self.data = {}
        if temp != None:
            self.data["temp"] = temp
        if W_Dat != None:
            self.data["Window Data"] = W_Dat

    def back(self):
        self.data = {}
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
            entry = tk.Entry(frame, width = 20)
            entry.grid(row = _row, column = _column)
            entries.append(entry)
            _row += 1
        return entries

    def basicWindowSetup(self, parent = None, back_button = True, help_msg = "", help_button = True):
        if parent == None:
            parent = self.root
        self.data["message"] = help_msg
        input_frm = tk.Frame(parent, relief = tk.SUNKEN, borderwidth = 3)
        input_frm.pack()
        button_frm = tk.Frame(parent)
        button_frm.pack()
        if back_button:
            path = os.path.join(self.FileDir, "back.png")
            back = tk.PhotoImage(file = path)
            button = tk.Button(input_frm, image = back, command = self.back, borderwidth=0)
            button.image = back
            button.config(height = 40, width = 50)
            button.grid(row = 0, column = 0)
        if help_button:
            path = os.path.join(self.FileDir, "help.png")
            help = tk.PhotoImage(file = path)
            button = tk.Button(input_frm, image = help, command = self.help_msg)
            button.image = help
            button.config(width = 20, height = 20)
            button.grid(row = 0, column = 1)
        return input_frm, button_frm

    def ListInputWindow(self, Field_Name, num_inputs, key_Name = None, startoffset = 1):
        Window = tk.Toplevel(self.root)
        self.activeWindow = Window
        window_height = 40 + num_inputs*30
        size = "250x" + str(int(window_height))
        Window.geometry(size)
        input_frm, button_frm = self.basicWindowSetup(parent = Window, back_button = False, help_button = False)
        labels = []
        for i in range(0,num_inputs):
            labels.append(Field_Name + str(i+startoffset))
        self.addLabelColumn(_column = 0,startrow = 1, labels = labels, frame = input_frm, _font = ("Arial Bold", 16))
        entries = self.addEntryColumn(1, 1, num_inputs, frame = input_frm)
        self.data["temp"] = self.entries
        self.entries = {}
        for i in range(0, num_inputs):
            key = Field_Name+str(int(i+startoffset))
            self.entries[key] = entries[i]
            if key_Name in self.data:
                self.entries[key].insert(0, self.data[key_Name][i])
        self.data["Window Data"]  = [Field_Name, key_Name, num_inputs]
        button = tk.Button(button_frm, text = "Submit", command = self.ListInput_submit)
        button.pack(side = tk.RIGHT)
        button = tk.Button(button_frm, text = "Clear Inputs", command = self.clear)
        button.pack(side=tk.LEFT)
        return None

    def ListInput_submit(self):
        Key_Name = self.data["Window Data"][1]
        self.data[Key_Name] = []
        n = 0
        for entry in self.entries:
            n+=1
            data = self.entries[entry].get()
            if data != "":
                try:
                    self.data[Key_Name].append(float(self.entries[entry].get()))
                except:
                    self.ErrorMsg("Invalid Value in Field " + str(n))
                    self.activeWindow.destroy()
                    self.entries = self.data["temp"]
                    self.ListInputWindow(self.data["Window Data"][0], Key_Name, self.data["Window Data"][2])
                    return -1
        print(self.data["Window Data"][0] + "List:")
        print(self.data[Key_Name])
        self.entries = self.data["temp"]
        if "masses" in self.data:
            self.entries["m0"].delete(0, tk.END)
            self.entries["m0"].insert(0, sum(self.data["masses"]))
        if len(self.data[Key_Name]) == 0:
            del self.data[Key_Name]
        self.activeWindow.destroy()
        
    def help_msg(self):
        mb.showinfo("Help", self.data["message"])

    def placeholder(self):
        msg_box = mb.showerror(title="Error", message="Function not yet Implemented")
    
    def run(self):
        self.root.mainloop()