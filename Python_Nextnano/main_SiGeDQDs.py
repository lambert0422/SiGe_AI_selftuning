import NextNanoClass as NNC
import numpy as np
def read_self_content():
    with open(__file__, 'r') as file:
        content = file.read()  # Read the entire content of the file

    return content

file_content = read_self_content()


InFileName = 'SiGe-DQDs'

# Vg = np.round(np.arange(0.5,-1.5,-0.01),3)

# ---------------- On Dell workstation -------------------------
# Disk = 'C:\\'
# NextnanoLoc = 'Program Files\\nextnano\\'
# OutputFilePath = "D:\\Nextnano\\SimulationFile\\PyTuneVar_newGateSize"
# ---------------- On ROG  -------------------------------------
Disk = 'C:\\'
LicenseFilePath = "C:\\Users\li244\OneDrive\Documents\\nextnano\License"
InFilePath = 'C:\\Users\li244\OneDrive\Desktop2\PostDoc-CAM\Simulation\\Nextnano'
NextnanoLoc = 'Program Files\\nextnano\\'
OutputFilePath = "C:\\Users\li244\Downloads\\NextnanoOutputfile"
# -------------------- On Alienware ----------------------------
Disk = 'D:\\'
NextnanoLoc = 'Nextnano\\'
LicenseFilePath = "D:\\Nextnano\\License"
InFilePath = 'D:\\Onedrive\\Desktop2\PostDoc-CAM\Simulation\\Nextnano'
OutputFilePath = "D:\\Nextnano\\Qi_temp\\nextnano_output"

# --------------------------------------------------------------11
NextnanoVer = '2021_03_15'
FirstImportFile = ''

# FirstImportFile = r'D:\Nextnano\SimulationFile\PyTuneVar_newGateSize\2023y02m09d\Swp_7_0_1_3.2_1_1.5_0.001_0.001_7.50e+18_1.10e+19_1.50e+15_430_200_1800_2500_1200_5000_1400_2.30e+17\S_NG__Bias_Gate_-1.37_02h36m44s\bias_00000\bandedges.fld'


VarDis = {
          # The 8 variables will be gathered into binary flag for example 000000011 = 3 and stored into the file name

          "M1D":                        [0],        # Whether to do 3D or 1D simulation 1 for 1D 0 for 3D
          "NextNanoRun":                [0],        # Whether the nextnano input file is running in nextnanomat(default be 0 in python)
          "FullOutput":                 [0],        # Whether to write full output including heavy hole etc, waste memory
          "InSwp":                      [0],        # Whether to have inside bias sweep, default 0 in python
          "BKD_donor":                  [0],        # Whether the background impurity is donor or acceptor(default accepter)
          "BKD":                        [0],        # Whether to add background doping
          "SurfaceCharge":              [0],        # Whether to add surface charge
          "QuantumSelfConsistent":      [0],        # Whether to run Poisson-quantum self-consistent simulation
          "THICK_AlOx_INSU":            [10]


          }

SweepVar = {'Bias_Gate2': [0.1,0.2,0.5,1],
            'Bias_Gate4': [0.05,0.1,0.2,0.3,0.4],
            'Bias_Gate5': [0.02,0.05,0.1],}
A = NNC.NN_inte(DiskLoc=Disk,NNLoc=NextnanoLoc,InputFilePathFolder=InFilePath,NNVer=NextnanoVer,
            OutputFilePath =OutputFilePath,AllVarDis =VarDis ,Block=True,InputFileName=InFileName,
            LicenseFilePath = LicenseFilePath,FirstImportFile = FirstImportFile,SweepVar = SweepVar)


