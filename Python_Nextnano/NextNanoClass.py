import nextnanopy as nn
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
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
# Nextnanopy version: 0.1.12
# Nextnano version: 2021_03_15
def TimeFormat(sec):
    if sec > 3600:
        Resulth = str(np.round((sec - sec % 3600) / 3600, 0)) +'h'
        Resultm = str(np.round((sec % 3600 - (sec % 3600) % 60) / 60, 0)) +'m'
        Results = ''
    elif sec > 60:
        Resulth = ''
        Resultm = str(np.round((sec - sec % 60) / 60, 0)) +'m'
        Results = str(np.round(sec % 60, 0)) +'s'
    else:
        Resulth = ''
        Resultm = ''
        Results = str(np.round(sec, 0)) +'s'
    return Resulth+Resultm+Results

def save_to_txt(filename, content):
    with open(filename, 'w') as file:
        file.write(content)
def winapi_path(dos_path, encoding=None):
    if (not isinstance(dos_path, str) and encoding is not None):
        dos_path = dos_path.decode(encoding)
    path = os.path.abspath(dos_path)
    if path.startswith(u"\\\\"):
        return u"\\\\?\\UNC\\" + path[2:]
    return u"\\\\?\\" + path
# Disable the print output of Python
def blockPrint():
    sys.stdout = open(os.devnull, 'w')
# Enable the print output
def enablePrint():
    sys.stdout = sys.__stdout__

def RecordPyFile(FilePath,FileName):
    with open(winapi_path(FilePath+ '\\'+FileName), 'w') as f2:
        with open(__file__, 'r') as f1:
            contents = f1.read()

            f2.write(contents)
            f2.write('\n')
    f1.close()
    f2.close()
#Get the time

# Datalist to store the carrier dentisy information and variable tuned

e = 1.60217663e-19
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

class NN_inte():

    def __init__(self,OutputFilePath,DiskLoc='D:\\',NNLoc = 'Nextnano\\',NNVer = '2021_03_15',
                 Overwrite = 1,AllVarDis = {},PLOT = False,
                 ImportInitialPotential = 1,Block = True,InputFileName = 'S_NG',InputFilePathFolder = '' ,
                 LicenseFilePath = 'C:\\Users\\li244\\OneDrive\\Documents\\nextnano\\License',FirstImportFile='',file_content = '',SweepVar = {}):

        self.PLOT =PLOT
        self.mainfilecontent = file_content
        self.FirstImportFile = FirstImportFile
        self.NNVer = NNVer
        self.PyFileRec = "PythonInterfaceFile"
        self.NextNanoFilePath = DiskLoc + NNLoc
        self.InputFileName = InputFileName
        self.OutputFilePath = OutputFilePath
        self.InputFileExtension = '.in'
        self.ExeFilePath = self.NextNanoFilePath + NNVer
        self.LicenseFilePath = LicenseFilePath
        self.ConfigFilePath = self.OutputFilePath + "\\Config\\"
        self.ConfigFile = self.ConfigFilePath +"MyConf.nextnanopy-config"
        self.InputFilePathFolder = InputFilePathFolder
        self.InputFilePath = self.InputFilePathFolder + "\\" + self.InputFileName + self.InputFileExtension



        if not os.path.exists(self.ConfigFilePath):
            os.makedirs(self.ConfigFilePath)


        self.NNConf()

        self.BLOCK = Block
        self.Overwrite = Overwrite
        self.SweepVar = SweepVar
        total_permutations = math.prod(len(v) for v in SweepVar.values())
        self.ImportInitialPotential = ImportInitialPotential
        keys, values = zip(*AllVarDis.items())
        self.permutations_dicts = [dict(zip(keys, v)) for v in itertools.product(*values)]
        self.TotalRunTime = len(self.permutations_dicts) * total_permutations
        self.var_sweep_self()
    def NNConf(self):
        nn.config.to_default()
        nn.config.set('nextnano++', 'exe', self.ExeFilePath + "\\nextnano++\\bin 64bit\\nextnano++_Intel_64bit.exe")
        nn.config.set('nextnano++', 'license', self.LicenseFilePath + '\\License_nnp.lic')
        nn.config.set('nextnano++', 'database', self.ExeFilePath + '\\nextnano++\\Syntax\\database_nnp.in')
        nn.config.set('nextnano++', 'outputdirectory', self.OutputFilePath)
        nn.config.set('nextnano++', 'threads', 1)

        nn.config.set('nextnano3', 'exe', self.ExeFilePath + '\\nextnano3\\Intel 64bit\\nextnano3_Intel_64bit.exe')
        nn.config.set('nextnano3', 'license', self.LicenseFilePath + '\\License_nnp.lic')
        nn.config.set('nextnano3', 'database', self.ExeFilePath + '\\nextnano3\\Syntax\\database_nn3.in')
        nn.config.set('nextnano3', 'outputdirectory', self.OutputFilePath)

        nn.config.set('nextnano.MSB', 'database',
                      self.NextNanoFilePath + "\\"+self.NNVer+"\\nextnano.MSB\\Syntax\\Materials.xml")

        nn.config.save(self.ConfigFile)

    def save_input_tofolder(self):


        self.my_sweep.DateTime = self.Date+'-'+self.Time
        vars = ''
        for i in self.my_sweep.var_sweep.keys():
            vars+=('__'+i)
        name_of_file = self.my_sweep.filename_only
        output_directory = self.my_sweep.config.get(section = self.my_sweep.product,option = 'outputdirectory')

        name = (name_of_file + '_sweep' + vars)
        if self.Overwrite==1:
            directory = nn.utils.misc.mkdir_if_not_exist(self.OutputFilePath + "\\INPUTFILES\\")
        else:
            directory = nn.utils.misc.mkdir_even_if_exists(self.OutputFilePath, "INPUTFILES")
        iteration_combinations = list(itertools.product(*self.my_sweep.var_sweep.values()))

        for combination in iteration_combinations:

            filename_end = '__'
            inputfile = nn.InputFile(fullpath = self.my_sweep.fullpath, configpath = self.my_sweep.configpath)
            if self.GeneralSetVarName["M1D"] == 0:
            # if self.Dim == 3:
                inputfile.set_variable("ImportFormat", "AVS")
            else:
                inputfile.set_variable("ImportFormat", "DAT")
            if self.GeneralSetVarName != None:
                CountSetVar = 0
                Flag = 0
                for SetVarName in self.GeneralSetVarName:
                    if CountSetVar <= 8:
                        Flag = Flag + self.GeneralSetVarName[SetVarName]*2**(7-CountSetVar)
                        CountSetVar += 1
                        VarDir = 'Swp_' + str(int(Flag))

                    Prevalue = inputfile.get_variable(SetVarName).value
                    inputfile.set_variable(SetVarName, self.GeneralSetVarName[SetVarName], comment='This is modified from '+str(Prevalue)+' to '+str(self.GeneralSetVarName[SetVarName]))
                    if CountSetVar > 8:
                        if self.GeneralSetVarName[SetVarName] >1e8 or (0<self.GeneralSetVarName[SetVarName] <1e-8):
                            VarDir = VarDir+'_'+"{:.2e}".format(self.GeneralSetVarName[SetVarName])
                        else:
                            VarDir = VarDir + '_' + str(self.GeneralSetVarName[SetVarName])
            for var_name, var_value in zip(self.my_sweep.var_sweep.keys(), combination):
                inputfile.set_variable(var_name, var_value, comment='THIS VARIABLE IS UNDER SWEEP')
                filename_end += '{}_{}'.format(var_name, var_value)
            inputfile.save(directory+"\\"+self.my_sweep.DateTime+"\\"+self.InputFileName + filename_end + '.in', overwrite = self.Overwrite)
            self.my_sweep.input_files.append(inputfile)
            self.my_sweep.VarDir = VarDir
            self.my_sweep.sweep_output_directory_user = self.OutputFilePath+"\\"+self.my_sweep.DateTime
            # nnswp.prepare_output(Overwrite)
            # output_directory = nnswp.sweep_output_directory

    def clear_line(self,n=1):
        LINE_UP = '\033[1A'
        LINE_CLEAR = '\x1b[2K'
        for i in range(n):
            print(LINE_UP, end=LINE_CLEAR)

    def excute_custumise_sweep(self, ImportFileName='potential', PreviousPotentialFile=''):
        self.my_sweep.prepare_output(self.Overwrite)
        # output_directory = nnswp.sweep_output_directory
        output_directory = self.my_sweep.sweep_output_directory_user
        if not self.my_sweep.input_files:
            warnings.warn('Nothing was executed in sweep! Input files to execute were not created.')
        elapsed_total = 0
        for Inputfile in self.my_sweep.input_files:
            if self.CountRunTime == 0:
                self.t_start = time.time()
            Inputfile.set_variable('Import', value=self.ImportInitialPotential)
            if self.ImportInitialPotential == 1:

                if PreviousPotentialFile != '':
                    Inputfile.set_variable('ImportPotentialFile', value="\"" + str(PreviousPotentialFile + "\""))
                    Inputfile.save(overwrite=self.Overwrite)
                    print('         Previous File imported')

            # if self.CountRunTime>1:
            #     self.clear_line(n = 6)
            print('Running:                 The <' + str(self.CountRunTime+1) + '> runs')


            if self.BLOCK:
                blockPrint()
            Inputfile.execute(outputdirectory=output_directory)
            self.CountRunTime = self.CountRunTime + 1
            self.GlobalCountRunTime = self.GlobalCountRunTime+1
            if self.BLOCK:
                enablePrint()
            elapsed = np.round(time.time() - self.t_start, 2)
            if self.CountRunTime > 0:
                self.t_start = time.time()
            elapsed_total = elapsed_total + elapsed
            Elapsed = TimeFormat(elapsed)
            LeftRuns = np.round(self.TotalRunTime - self.GlobalCountRunTime, 0) # need to update if the simulation time is more then 24 hours, it will reset time
            TimeSpend = np.round(time.time() - self.GlobalStartTime, 2)

            print(
                  'TIME-ONESWP:             ' + Elapsed + '(TOTAL: ' + TimeFormat(TimeSpend) + ')')
            print('RUNTIME LEFT:            ' + str(LeftRuns) )
            print('EST TIME LEFT:           ' + TimeFormat(LeftRuns * elapsed_total / self.GlobalCountRunTime) )

            output_directory_Folder = Inputfile.folder_output
            My_datafolder = nn.DataFolder(output_directory_Folder)

            Potential = My_datafolder.find('potential', deep=True)
            band = My_datafolder.find('bandedge', deep=True)
            density = My_datafolder.find('integrated_density_electron', deep=True)
            IntegratedDensity = My_datafolder.find('integrated_density_electron.dat', deep=True)
            for FN in IntegratedDensity:
                FN = winapi_path(FN)
                TwoD_PotentialData = nn.DataFile(FN, product='nextnano++')
                self.data_Frame_row = []
                if self.TitleRun:
                    self.data_Frame_row.append('Output file path')
                    for data in self.GeneralSetVarName:
                        self.data_Frame_row.append(data)
                    for data in TwoD_PotentialData.data:
                        self.data_Frame_row.append(data.label)
                    self.data_Frame.append(self.data_Frame_row)
                    self.data_Frame_row = []
                    self.TitleRun = False


                # FN_Excel = re.search(self.my_sweep.DateTime+'(.*)' + self.InputFileName, FN)
                self.data_Frame_row.append(output_directory_Folder)
                ZeroGate = False
                Great2 = False
                for data in self.GeneralSetVarName:
                    self.data_Frame_row.append(self.GeneralSetVarName[data])
                for data in TwoD_PotentialData.data:
                    self.data_Frame_row.append(float(data.value))
                    if data.label == 'Gate_bias (V)' and float(data.value) == 0:
                        ZeroGate = True
                    if ZeroGate:
                        if data.label == 'region_8 (carriers)' and float(data.value) > 11:
                            Great2 = True

                if Great2:
                    self.data_Frame_row.append('>11')
                else:
                    self.data_Frame_row.append('')
                self.data_Frame.append(self.data_Frame_row)

            df = pd.DataFrame(self.data_Frame)
            df.to_excel(output_directory + '\\' + self.InDensity_Excel, index=False)
            try:
                df.to_excel(output_directory + '\\Backup_' + self.InDensity_Excel, index=False)
            except:
                pass
            HalfMeasureSize = Inputfile.get_variable('HalfMesuSize').value
            Depth_2DEG = Inputfile.get_variable('THICK_Ge_QW').value
            self.MeasuVolume = (Depth_2DEG*4*HalfMeasureSize**2)*1e-21
            if self.PLOT:
                self.plot_output(Potential)
                self.plot_output(band, Note='_zoomed', Xlim=(-40, 10), Ylim=(-0.1, 0.5))
                self.plot_output(band)
                self.plot_output(density)
            if self.ImportInitialPotential == 1:
                if self.GeneralSetVarName["M1D"] == 0:
                # if self.Dim == 3:
                    PreviousPotentialFile = output_directory_Folder + "\\bias_00000\\" + ImportFileName + ".fld"
                else:

                    PreviousPotentialFile = output_directory_Folder + "\\bias_00000\\" + ImportFileName + ".dat"
            A = self.my_sweep.var_sweep
            Checker = 0
            for setV in A:
                if A[setV][0] != Inputfile.get_variable(setV).value:
                    Checker = 1
            if Checker == 0:
                self.my_sweep.FirstSwpOutPotential = PreviousPotentialFile

    def plot_output(self,FileList, Note='', Xlim=None, Ylim=None, combine1D=0, SHOW=0, SAVE=1):

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
                        for data in TwoD_PotentialData.data:
                            if data.label == 'Gate_bias (V)':
                                print(data.label + ':           ' + str(float(data.value)))
                            if data.label == 'region_8 (carriers)':

                                print('Density (cm^3)' + ':          ' + "{:.2e}".format(float(data.value)/float(self.MeasuVolume)))
                        print('-----------------------------------------------------------------------------------------------')


                else:
                    for z in TwoD_PotentialData.variables:
                        plt.figure()
                        fig, ax = plt.subplots(1)
                        pcolor = ax.pcolormesh(x.value, y.value, z.value.T, shading='auto')

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

    def var_sweep_self(self):

        VarCount = 0
        self.GlobalCountRunTime = 0
        for VarDis in self.permutations_dicts:
            self.TitleRun = True
            self.GlobalStartTime = time.time()
            now = datetime.now()
            self.Time = now.strftime("%Hh%Mm%Ss")
            self.Date = now.strftime("%YY%mM%dD")

            self.data_Frame = []
            self.data_Frame_row = []
            self.CountRunTime = 0
            self.InDensity_Excel = self.Date + '-' + self.Time + '.xlsx'

            self.GeneralSetVarName = VarDis
            if VarCount != 0:
                print('------------------------------ Different setting ---------------------------------')
            self.my_sweep = nn.Sweep(self.SweepVar, self.InputFilePath, self.ConfigFile)

            self.save_input_tofolder()
            if not os.path.exists(self.my_sweep.sweep_output_directory_user):
                os.makedirs(self.my_sweep.sweep_output_directory_user)
            RecordPyFile(self.my_sweep.sweep_output_directory_user , self.PyFileRec+".txt")
            save_to_txt(self.my_sweep.sweep_output_directory_user +"\\main.txt",self.mainfilecontent)

            self.excute_custumise_sweep( PreviousPotentialFile=self.FirstImportFile)

            for i in range(len(self.data_Frame_row)):
                self.data_Frame_row[i] = ''
            self.data_Frame.append(self.data_Frame_row)

            self.FirstImportFile = self.my_sweep.FirstSwpOutPotential
            # my_datafolder = nn.DataFolder(my_sweep2.sweep_output_directory_user)
            # FileList = my_datafolder.find('bandedge_x_1d_Gated.dat', deep=True)
            # plot_output(FileList, Note = '', Xlim = (-50,20),  combine1D = 1, SHOW = 0, SAVE = 0)
            # FileList = my_datafolder.find('bandedge_x_1d.dat', deep=True)
            # plot_output(FileList, Note='', Xlim=(-50,20),  combine1D=1, SHOW=0, SAVE=0)
            VarCount = VarCount +1
        print('------------------------------ All finished ---------------------------------')
