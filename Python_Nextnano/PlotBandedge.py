import nextnanopy as nn
import matplotlib.pyplot as plt
import pandas as pd
import os
import warnings
import numpy as np
import itertools
import scipy
from scipy.interpolate import interp2d
import re
from datetime import datetime
import itertools
from itertools import chain, combinations
import sys, os
import time


def plot_output(FileList, Note='', Xlim=None, Ylim=None, combine1D=0, SHOW=1, SAVE=0,F1D = 1):
    for FolderName in FileList:
        FolderName = winapi_path(FolderName)
        TwoD_PotentialData = nn.DataFile(FolderName, product='nextnano++')

        try:
            x = TwoD_PotentialData.coords['x']
            y = TwoD_PotentialData.coords['y']
            z = TwoD_PotentialData.coords['z']
        except:
            try:
                FolderName_Plot = FolderName.replace('.dat', '.png')
                FolderName_Plot = FolderName_Plot.replace('.fld', '.png')
                x = TwoD_PotentialData.coords['x']
                y = TwoD_PotentialData.coords['y']
            except:
                try:
                    FolderName_Plot = FolderName.replace('.dat', Note + '.png')
                    FolderName_Plot = FolderName_Plot.replace('.fld', Note + '.png')
                    x = TwoD_PotentialData.coords['x']
                    for y in TwoD_PotentialData.variables:
                        if y.name == 'Gamma':
                            X = x.value
                            Y = y.value
                        if y.name == 'electron_Fermi_level':
                            X_mu = x.value
                            Y_mu = y.value
                        if y.name == 'Electron_density':
                            X = x.value
                            Y = y.value
                        plt.figure()
                        plt.plot(x.value, y.value)

                        plt.xlabel(x.label)
                        plt.ylabel(y.label)
                        if Xlim != None:
                            plt.xlim(Xlim)
                            A = x.value
                            min_index1 = np.argmin(abs(x.value - Xlim[0]))
                            min_index2 = np.argmin(abs(x.value - Xlim[1]))
                            plt.ylim((min(y.value[min_index1:min_index2]), max(y.value[min_index1:min_index2])))
                        if SAVE == 1:
                            plt.savefig(FolderName_Plot)
                        if SHOW == 1:
                            plt.show()
                except:

                    print('-------------------------------------------------------------------------------------------')


            else:
                for z in TwoD_PotentialData.variables:
                    plt.figure()
                    fig, ax = plt.subplots(1)
                    pcolor = ax.pcolormesh(x.value, y.value, z.value.T, shading='auto')
                    if z.name == 'Gamma':
                        X = x.value
                        Y = y.value
                        Z = z.value.T
                    cbar = fig.colorbar(pcolor)
                    cbar.set_label(z.label)
                    ax.set_xlabel(x.label)
                    ax.set_ylabel(y.label)
                    fig.tight_layout()
                    plt.xlim((-1600, 2400))
                    if SAVE == 1:
                        fig.savefig(FolderName_Plot)
                    if SHOW == 1:
                        plt.show()
        # else:
        #     warnings.warn('3D data is not plotted')
    plt.close('all')
    if combine1D == 1:
        plt.show()
    if F1D == 1:
        return X,Y
    else:
        return X,Y,Z
def winapi_path(dos_path, encoding=None):
    if (not isinstance(dos_path, str) and encoding is not None):
        dos_path = dos_path.decode(encoding)
    path = os.path.abspath(dos_path)
    if path.startswith(u"\\\\"):
        return u"\\\\?\\UNC\\" + path[2:]
    return u"\\\\?\\" + path


def LinePlot2Origin(filefolder = r'C:\Users\li244\OneDrive\Documents\nextnano\Output\F1D_S2DEG3',
                    filename='bandedges.dat',savename='BandEdge.txt',Title = ["x(nm)", "Bandedge(eV)"],F1D = 1):

    My_datafolder = nn.DataFolder(filefolder)
    Potential = My_datafolder.find(filename, deep=True)

    if F1D == 1:
        X,Y = plot_output(Potential,SHOW = 1)
        plt.plot(X, Y)

        plt.show()
        A = np.trapz(Y, X) #1e11 cm^-2
        print('Electron density:'+"{:.3e}".format(A*1e11)+'cm^-2')
        Data = np.vstack((X, Y))

        Data = np.vstack((Title, Data.T))
        DataF = pd.DataFrame(Data)
        DataF.to_csv(savename, index=False, header=False)
    else:
        X,Y,Z = plot_output(Potential,SHOW = 0,F1D = F1D)
        Data = np.vstack((X.T, Z))
        Y = list(Y)
        Y.insert(0, 0)
        Y = np.array(Y)
        Data = np.vstack((Y, Data.T))
        DataF = pd.DataFrame(Data)
        DataF.to_csv(savename, index=False, header=False)
# LinePlot2Origin()


Filefolder = r'C:\Users\li244\OneDrive\Documents\nextnano\Output\2023Y05M04D-10h50m13s'

LinePlot2Origin( filefolder = Filefolder,filename='density_electron.dat',savename='Density.txt',Title = ["x(nm)", "Density(1e18/cm^3)"])

