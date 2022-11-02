import matplotlib.pyplot as plt
import os.path
from datetime import datetime
import numpy as np

class Plotter(object):
    def __init__(self, DataDir):
        self.DataDir = DataDir

    def hist2D(self, dataX, dataY, show = 1, savefile = 0, filename = None):
        """
        Method for plotting 2D Histogram of Data

        Inputs:
            DataX: X-Axis data
            DataY: Y-Axis data
        
        Optional:
            show: Flag for displaying the plot
            savefile: Flag for saving the plot to file
            filename: Filename. If no filename is specifified, it is set to be "plot_yyyymmdd-hh_mm_ss"

        Returns
            None
        """
        if len(dataX) != len(dataY):
            #Check for appropriate dimensions (Matplotlib Error message for mismtached dimensions not useful)
            raise ValueError("X and Y must have same length")
        #Plot the histogram
        plt.hist2d(self, dataX, dataY,  bins=(34,20), cmap=plt.cm.jet)
        if show:
            #Show the Plot
            plt.show()
        if savefile:
            #Determine the Base directory
            if filename != None:
                #Create path by joining base, data and filename
                filename = os.path.join(self.DataDir, filename)
                plt.savefig(filename)
            else:
                #if no filename specified, filename becomes current date and time
                filename = "plot" + datetime.now().strftime("%Y%m%d-%H_%M_%S")
                filename = os.path.join(self.DataDir, filename)
                plt.savefig(filename)


    def colourMap(self, dataX, dataY, dataZ, show = 1, savefile = 0, filename = None):
        """
        Function for creating 2D colourmap plot
        
        Inputs:
            dataX: X-Axis
            dataY: Y-Axis
            dataZ: mapped data
            
            Such that dataZ is mapped to dataX and dataY: dataZ[y, x] = f(dataX[x], dataY[y])
        Optional:
            show: Flag for displaying the plot
            savefile: Flag for saving the plot to file
            filename: Filename. If no filename is specifified, it is set to be "plot_yyyymmdd-hh_mm_ss"

        returns
            None

        """
        if len(dataX) != len(dataZ[0]) or len(dataY) != len(dataZ[2]):
            #Check for Data dimensions
            raise TypeError("Dimensions of DataZ must be (len(x), len(y))")
        #create the plot

        #Take 1-D input data dataX and dataY and tranform it into 2-D grid such that the gridpoints
        # are defined in accordance with pcolormesh documentation:
        #(X[i+1, j], Y[i+1, j])       (X[i+1, j+1], Y[i+1, j+1])
        #                  +-----+
        #                  |     |
        #                  +-----+
        #(X[i, j], Y[i, j])       (X[i, j+1], Y[i, j+1])
        #(see pcolormesh documentation)
        X = np.zeros((len(dataY), len(dataX)))
        Y = np.zeros((len(dataY), len(dataX)))
        for i in range(0, len(dataY)):
            for n in range(0, len(dataX)):
                X[i, n] = dataX[n] #X maps the columns of dataZ
                Y[i, n] = dataY[i] #Y maps the lines of dataZ
        fig, ax = plt.subplots()
        ax.pcolormesh(X, Y, dataZ,edgecolors = "black", cmap = "YlGn")
        #Note: Saveing MUST happen before displaying, since displaying overwrites fig
        #out of sequence ops cause a blank image to be saved instead of the generated plot
        if savefile:
            #Determine the Base directory
            if filename != None:
                #Create path by joininh base, data and filename
                filename = os.path.join(self.DataDir, filename)
                plt.savefig(filename)
            else:
                #if no filename specified, filename becomes current date and time
                filename = "plot" + datetime.now().strftime("%Y%m%d-%H_%M_%S")
                filename = os.path.join(self.DataDir, filename)
                plt.savefig(filename)
        if show:
            #Show the Plot
            plt.show()
        plt.close("all")

    def contour(self, dataX, dataY, dataZ, show = 1, savefile = 0, filename = None):
        """
        Function for creating 2D contour plot
        
        Inputs:
            dataX: X-Axis
            dataY: Y-Axis
            dataZ: mapped data
            
            Such that dataZ is mapped to dataX and dataY: dataZ[y, x] = f(dataX[x], dataY[y])
        Optional:
            show: Flag for displaying the plot
            savefile: Flag for saving the plot to file
            filename: Filename. If no filename is specifified, it is set to be "plot_yyyymmdd-hh_mm_ss"

        returns
            None

        """
        if len(dataX) != len(dataZ[0]) or len(dataY) != len(dataZ[2]):
            #Check for Data dimensions
            raise TypeError("Dimensions of DataZ must be (len(x), len(y))")
        #create the plot

        #Take 1-D input data dataX and dataY and tranform it into 2-D grid such that the gridpoints
        # are defined as:
        #(X[i+1, j], Y[i+1, j])       (X[i+1, j+1], Y[i+1, j+1])
        #                  +-----+
        #                  |     |
        #                  +-----+
        #(X[i, j], Y[i, j])       (X[i, j+1], Y[i, j+1])
        #(see pyplot.pcolormesh documentation)
        X = np.zeros((len(dataY), len(dataX)))
        Y = np.zeros((len(dataY), len(dataX)))
        for i in range(0, len(dataY)):
            for n in range(0, len(dataX)):
                X[i, n] = dataX[n] #X maps the columns of dataZ
                Y[i, n] = dataY[i] #Y maps the lines of dataZ
        fig, ax = plt.subplots()
        ax.contourf(X, Y, dataZ, levels= 100)
        #Note: Saveing MUST happen before displaying, since displaying overwrites fig
        #out of sequence ops cause a blank image to be saved instead of the generated plot
        if savefile:
            #Determine the Base directory
            base = os.path.dirname(__file__) #File Directory (utils)
            base = os.path.split(base)[0]    #Splits directory at utils, returns project dir
            if filename != None:
                #Create path by joininh base, data and filename
                filename = os.path.join(self.DataDir, filename)
                plt.savefig(filename)
            else:
                #if no filename specified, filename becomes current date and time
                filename = "plot_" + datetime.now().strftime("%Y%m%d-%H_%M_%S")
                filename = os.path.join(self.DataDir, filename)
                plt.savefig(filename)
        if show:
            #Show the Plot
            plt.show()
        #close all plots
        plt.close("all")


    def plot2D(self, dataX, dataY,xLab = None, yLab = None, show = 1, savefile = 0, filename = None, dataY_List = None):
        """
        Function for creating 2D colourmap plot
        
        Inputs:
            dataX: X-Axis
            dataY: Y-Axis

        Optional:
            dataY_List: List of Y datasets for comparative plotting
            show: Flag for displaying the plot
            savefile: Flag for saving the plot to file
            filename: Filename. If no filename is specifified, it is set to be "plot_yyyymmdd-hh_mm_ss"

        returns
            None
        """
        fig, ax = plt.subplots()
        ax.plot(dataX, dataY, linewidth=2.0)
        ax.set_xlim(dataX.min(), dataX.max())
        #ax.set_ylim(0, 1000000)
        if xLab != None:
            ax.set_xlabel(xLab)
        if yLab != None:
            ax.set_ylabel(yLab)
        ax.set_ylabel(yLab)
        # ax.set_yscale("log")
        ax.grid(True)
        if dataY_List != None:
            for data in dataY_List:
                ax.plot(dataX, data, linewidth=2.0)
        if savefile:
            #Determine the Base directory
            base = os.path.dirname(__file__) #File Directory (utils)
            base = os.path.split(base)[0]    #Splits directory at utils, returns project dir
            if filename != None:
                #Create path by joining base, data and filename
                filename = os.path.join(base, "data", filename)
                plt.savefig(filename)
            else:
                #if no filename specified, filename becomes current date and time
                filename = "plot_" + datetime.now().strftime("%Y%m%d-%H_%M_%S")
                filename = os.path.join(base, "data", filename)
                plt.savefig(filename)
        if show:
            #Show the Plot
            plt.show()
        #close all plots
        plt.close("all")