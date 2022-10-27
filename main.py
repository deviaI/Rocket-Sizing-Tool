from utils.Calculators import Calculator
from utils.Exporter import Exporter
import os.path

def main():


    base = os.path.dirname(__file__)
    z = Calculator(250, 9000)

    print(z.calcSinglePoint(2, [311, 343], 250))


    Exp = Exporter(os.path.join(base, "data"))
    results = z.f_reverse(2, 12500, 250, [311, 343]) #Calculate achievable delta v with a 12.5 ton, 2 stage Rocket with Isps of 311s and 343s and default structure factor of 12% for a 250kg payload
    print(results)
    
    # Generates and save a bunch of default comparative .csv files For a 250kg Payload and Delta V Target of 9 km/s:
    results = z.TwoDAlt(1, "Isp", 280, 100, 0.01, 0.15, 250, 9000, 1e6)                             #SSTO with 280s Isp Engine   
    fName =  "SSTO_ISP= 280s"
    fInfo = "ISP = 280s; Configuration: SSTO"
    Exp.ExportData(results, fName, "csv")
    Exp.AppendData(fInfo, fName + ".csv")

    results = z.TwoDAlt(1, "Isp", 311, 100, 0.01, 0.15, 250, 9000, 1e6)                             #SSTO with 311s Isp Engine
    fName =  "SSTO_ISP=311" 
    fInfo = "ISP = 311; Configuration: SSTO"
    Exp.ExportData(results, fName, "csv")
    Exp.AppendData(fInfo, fName + ".csv")

    results = z.TwoDAlt(2, "Isp", 280, 100, 0.01, 0.15, 250, 9000, 1e6)                             #2 Stage Rocket with 280s Isp Engines
    fInfo = "Configuration: 2-Stage, mass Stage 1 = 2 x mass Stage 2; Isp Stage 1=Isp Stage 2: 280s" 
    fName = "2Stage_ISP=280s"
    Exp.ExportData(results, fName, "csv")
    Exp.AppendData(fInfo, fName + ".csv")

    results = z.TwoDAlt(2, "Isp", 280, 100, 0.01, 0.15, 250, 9000, 1e6, isp_2 = 320)                #2 Stage Rocket with 280s Isp First Stage Engine, 320s Second Stage Engine
    fInfo = "Configuration: 2-Stage, mass Stage 1 = 2 x mass Stage 2; Isp Stage 1= 280s; Isp Stage 2 = 320s" 
    fName = "2Stage_ISP1==280s;ISP2=320s"   
    Exp.ExportData(results, fName, "csv")
    Exp.AppendData(fInfo, fName + ".csv")

    results = z.TwoDAlt(2, "Isp", 311, 100, 0.01, 0.15, 250, 9000, 1e6, isp_2 = 343)                #2 Stage Rocket with 311s Isp First Stage Engine, 343s Second Stage Engine [Rocket Lab Electron Engines]
    fInfo = "Configuration: 2-Stage, mass Stage 1 = 2 x mass Stage 2; Isp Stage 1= 311s; Isp Stage 2 = 343s" 
    fName = "2Stage_ISP1==311s;ISP2=343s"   
    Exp.ExportData(results, fName, "csv")
    Exp.AppendData(fInfo, fName + ".csv")

    results = z.TwoDAlt(3, "Isp", 280, 100, 0.01, 0.15, 250, 9000, 1e6)                             #3 Stage Rocket with 280s Isp Engines
    fInfo = "Configuration: 3-Stage, mass Stage 1 = 2 x mass Stage 2; Isp Stage 1=Isp Stage 2=Isp Stage 3=320s" 
    fName = "3Stage_ISP==320s"  
    Exp.ExportData(results, fName, "csv")
    Exp.AppendData(fInfo, fName + ".csv")

    results = z.TwoDAlt(3, "Isp", 311, 100, 0.01, 0.15, 250, 9000, 1e6, isp_2 = 343, isp_3 = 410)   #3 Stage Rocket with 311, 343s, 410s Isp Engines on respective Stages [Electron Engines + RL10 Upper Stage]
    fInfo = "Configuration: 3-Stage, mass Stage 1 = 2 x mass Stage 2; Isp Stage 1= 311s; Isp Stage 2 = 343s, Isp Stage 3 = 410s" 
    fName = "2Stage_ISP1==280s;ISP2=320s;ISP3=410s"  
    Exp.ExportData(results, fName, "csv")
    Exp.AppendData(fInfo, fName + ".csv")

    results = z.Optimised2Stage(311, 343, 100, 0.01, 0.15, 250, 9000, 1e6)                          #2 Stage Rocket with optimised First Stage/Second Stage Mass Ratio and Electron Engine Specs
    fName = "2-Stage Opt.;ISP1=311s;ISP2=343s"
    fInfo = " Configuration: 2-Stage optimised , Relative Stage Sizing = third column, Isp Stage 1=311s; Isp Stage 2=343s"
    Exp.ExportData(results, fName, "csv")
    Exp.AppendData(fInfo, fName + ".csv")


if __name__ == '__main__':

    main()
