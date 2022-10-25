import ziolkowsky as z

def main():


    results = z.f_reverse(2, 12500, 250, [311, 343]) #Calculate achievable delta v with a 12.5 ton, 2 stage Rocket with Isps of 311s and 343s and default structure factor of 12% for a 250kg payload
    print(results)
    # Generates a bunch of default comparative .csv files For a 250kg Payload and Delta V Target of 9 km/s:
    z.TwoDAlt(1, "Isp", 280, 100, 0.01, 0.15, 250, 9000, 1e6)                             #SSTO with 280s Isp Engine
    z.TwoDAlt(1, "Isp", 311, 100, 0.01, 0.15, 250, 9000, 1e6)                             #SSTO with 311s Isp Engine
    z.TwoDAlt(2, "Isp", 280, 100, 0.01, 0.15, 250, 9000, 1e6)                             #2 Stage Rocket with 280s Isp Engines
    z.TwoDAlt(2, "Isp", 280, 100, 0.01, 0.15, 250, 9000, 1e6, isp_2 = 320)                #2 Stage Rocket with 280s Isp First Stage Engine, 320s Second Stage Engine
    z.TwoDAlt(2, "Isp", 311, 100, 0.01, 0.15, 250, 9000, 1e6, isp_2 = 343)                #2 Stage Rocket with 311s Isp First Stage Engine, 343s Second Stage Engine [Rocket Lab Electron Engines]
    z.TwoDAlt(3, "Isp", 280, 100, 0.01, 0.15, 250, 9000, 1e6)                             #3 Stage Rocket with 280s Isp Engines
    z.TwoDAlt(3, "Isp", 311, 100, 0.01, 0.15, 250, 9000, 1e6, isp_2 = 343, isp_3 = 410)   #3 Stage Rocket with 311, 343s, 410s Isp Engines on respective Stages [Electron Engines + RL10 Upper Stage]
    z.Optimised2Stage(311, 343, 100, 0.01, 0.15, 250, 9000, 1e6)                          #2 Stage Rocket with optimised First Stage/Second Stage Mass Ratio and Electron Engine Specs

if __name__ == '__main__':

    main()
