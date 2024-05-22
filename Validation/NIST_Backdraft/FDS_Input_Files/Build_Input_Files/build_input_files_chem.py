# Based on the build_input_files script of the same purpose in Deep_Fuel_Beds.
# This script gene


import math
import pandas as pd
import os
import numpy as np

# read input data to dataframe
df = pd.read_csv('./Cases_Backdraft_Chem.csv', header=0)

# initalize matrix to write final dataframe
FINAL = []

# function gets rid of .0 s or replaces . with p
def rep_0(str_in):
    if str_in[2:]=='.0':
        str_out = str_in[:2]
    else:
        str_out = str_in.replace('.','p')
    return(str_out)

Grid_Size_Vect   = ['0p25cm_2432','0p5cm_608','1cm_304']
#Grid_Size_Vect   = ['1cm_304']

# read data from the Backdraft_Two_Zone_init.csv file
for Grid_Size in Grid_Size_Vect:
    Grid_Size_File   = "mesh_" + Grid_Size + ".fds"        
    for irow in df.index:
        Fuel        = str(df.loc[irow,'Fuel'])
        Reac_File   = Fuel.lower() + "_2eq_fast.fds"
        Fire_Size   = str(df.loc[irow,'Fire_Size_kW'])
        Fire_Size   = rep_0(Fire_Size)
        Flow_Time   = str(df.loc[irow,'Fuel_Flow_Time_s'])
        Init_File1   = str(df.loc[irow,'Init_File'])
        Init_File   = Fuel + Fire_Size + "kW_" + Flow_Time + "s_"+ Init_File1 + "_init.fds"

        # make the row for that file
        CHID = Fuel + Fire_Size + "kW_" + Flow_Time + "s_" + Init_File1 + '_' + Grid_Size
        filename = CHID + ".fds"
        TITLE = "NIST Backdraft experiment : " + Fuel + " - " + Fire_Size + "kW - " + Flow_Time + "s - Temp var " + Init_File1 + "."
        paramline = [filename] + [CHID] + [TITLE] + [Grid_Size_File] + [Init_File] + [Reac_File] + [Fuel.upper()]
        # add paramline to FINAL matrix for each irow
        FINAL = FINAL + [paramline]


# make the header list
topline = ["NIST_BackdraftTemplate_chem.fds"] + ["paramCHID"] + ["paramTITLE"]
for i in range(3,7):
    topline = topline + ["param"+str(i)]

# make fdout and wirte paramfile
dfout = pd.DataFrame(FINAL, columns=topline)
dfout.to_csv('paramfile.csv', index=False)

# build input files: run swaps.py
os.system('python ../../../../Utilities/Input_File_Tools/swaps.py')

# Move input files up one level to FDS_Input_Files
os.system('mv Propane*.fds ../')
