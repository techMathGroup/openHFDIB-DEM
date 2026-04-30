#!/usr/bin/python

#FILE DESCRIPTION=======================================================
# Python script to automatically generate cases for the settling test
# case
#
# the script is used to prepare a set of benchmarks that are to be ran
# on a cluster at it.cas.cz 
#

#LICENSE================================================================
#  case_constructor.py
#
#  Copyright 2025 Martin Isoz <isozm@it.cas.cz>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

#########DO NOT EDIT####################################################
import os
import math
import shutil as sh

#IMPORT BLOCK-CUSTOM====================================================

#CUSTOM FUNCTIONS=======================================================
def update_file(file_path,id_strings,replace_strings):
    with open(file_path, 'r') as file:
        data = file.readlines()
        
    for j in range(len(idStr)):
        for i in range(len(data)):
            fInd = data[i].find(idStr[j])
            if fInd>-1:
                data[i] = data[i][:fInd] + pVals[j] + '\n'              #always adds end of line
    
    with open(file_path, 'w') as file:
        file.writelines( data )

#=======================================================================
#							EDITABLE
#=======================================================================

#-----------------------------------------------------------------------
# I/O DATA
#-----------------------------------------------------------------------
baseCase = '../10_settling_base_case/'
outFolder= '../20_testCases/'
caseName = 'settling'

#-----------------------------------------------------------------------
# CASE/RUN SETTINGS
#-----------------------------------------------------------------------
coresComp   = [1,4,1]                                                   #we are working with simple decomposition
endTime     = 10                                                        #simulation endTime
wrInt       = 1e-1                                                      #simulation write interval
maxCo       = 0.8                                                       #maximum Courant number
initDeltaT  = 1e-6                                                      #initial simulation timestep

node_spec   = "kraken-x4"                                               #on which node I want to run the case

#-----------------------------------------------------------------------
# GEOMETRY DATA                                                         #case geometry is fixed
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# MESH DATA
#-----------------------------------------------------------------------
nCellsHor = 70                                                          #resolution in horizontal directions
nCellsVer = 2*nCellsHor                                                 #resolution in vertical direction
# Note (MI): vertical size of the domain is twice the horizontal size

#-----------------------------------------------------------------------
# SOLUTION SETTINGS
#-----------------------------------------------------------------------
# system/fvSolution -> PIMPLE
nOuterCorrectors = 3
nCorrectors      = 3
# constant/HFDIBDEMDict
nDirectForcingIterations = 10
immersedUTolerance       = 1e-4

#-----------------------------------------------------------------------
# CASE PARAMETERS
#-----------------------------------------------------------------------
nuF  = 1e-3
rhoF = 1e3
rhoS = 6*rhoF


#########PREFERABLY DO NOT EDIT#########################################
#-----------------------------------------------------------------------
# AUXILIARY COMPUTATIONS
#-----------------------------------------------------------------------
nCores      = math.prod(coresComp)                                      #number of cores to run the case on

job_name    = caseName + "_nuF_%g_rhoF_%g_rhoRatio_%f"%(nuF,rhoF,rhoS/rhoF)

caseName    = job_name + "/"

caseDir     = outFolder + caseName

#-----------------------------------------------------------------------
# GATHERING AUXILIARY SCRIPTS
#-----------------------------------------------------------------------
pyFolder = '../00_python_codes/'                                        #source directory for python codes
pyList   = [                                                            #python files to copy to the created folder
    'case_constructor',
]
# Note: these files are copied into folder caseName/ for
#       future reference
# Note: update these files if needed - in case of new code versions

shFolder = '../01_bash_codes/'                                          #source directory for bash codes
shList   = [
    'cleanAllCases',
    'queueAllCases',
]

#-----------------------------------------------------------------------
# DATA GATHERING
#-----------------------------------------------------------------------
if os.path.isdir(caseDir):                                              #ensure, that the caseDir is clear
    sh.rmtree(caseDir)

sh.copytree(baseCase,caseDir)

# -- copy the data from the folder with python codes
for pyCode in pyList:
    sh.copyfile(pyFolder + pyCode + '.py',caseDir + pyCode + '.py')
    
# -- copy the data from the folder with base codes
# Note (MI): we are using these codes to do batch job submitting and
#            cleaning
for shCode in shList:
    sh.copyfile(shFolder + shCode + '.sh',outFolder + shCode + '.sh')

#-----------------------------------------------------------------------
# ALLRUN-SLURM PREPARATION
#-----------------------------------------------------------------------
# Note (MI): if you want to run the case on a different cluster than
#            kraken@it.cas.cz, you need to modify "Allrun-slurm" in
#            the baseCase folder and also change this accordingly
pVals   = [
    "#SBATCH -J %s"%(job_name),
    "#SBATCH --ntasks=%d"%(nCores),
    "#SBATCH --nodelist=%s"%(node_spec),
]  
idStr   = [
    "#SBATCH -J AAAA",
    "#SBATCH --ntasks=BBBB",
    "#SBATCH --nodelist=CCCC",
]

update_file(caseDir + './Allrun-slurm',idStr,pVals)
    
#########DO NOT EDIT####################################################
    
#-----------------------------------------------------------------------
# BLOCKMESHDICT PREPARATION
#-----------------------------------------------------------------------
pVals   = [
    "nCellsHor %d;"%(nCellsHor),
    "nCellsVer %d;"%(nCellsVer),
]  
idStr   = [
    "nCellsHor AAAA;",
    "nCellsVer BBBB;",
]

update_file(caseDir + './system/blockMeshDict',idStr,pVals)

#-----------------------------------------------------------------------
# CONTROLDICT PREPARATION
#-----------------------------------------------------------------------
pVals   = [
    "endTime         %g;"%(endTime),
    "deltaT          %g;"%(initDeltaT),
    "writeInterval   %g;"%(wrInt),
    "maxCo           %g;"%(maxCo),
]  
idStr   = [
    "endTime         AAAA;",
    "deltaT          BBBB;",
    "writeInterval   CCCC;",
    "maxCo           DDDD;",
]

update_file(caseDir + './system/controlDict',idStr,pVals)

#-----------------------------------------------------------------------
# FVSOLUTION PREPARATION
#-----------------------------------------------------------------------
pVals   = [
    "nOuterCorrectors    %d;"%(nOuterCorrectors),
    "nCorrectors         %d;"%(nCorrectors),
]  
idStr   = [
    "nOuterCorrectors    AAAA;",
    "nCorrectors         BBBB;",
]

update_file(caseDir + './system/fvSolution',idStr,pVals)

#-----------------------------------------------------------------------
# DECOMPOSEPARDICT PREPARATION
#-----------------------------------------------------------------------
pVals   = [
    "numberOfSubdomains  %d;"%(nCores),
    "n           (%d %d %d);"%tuple(coresComp),
]  
idStr   = [
    "numberOfSubdomains  AAAA;",
    "n           (BBBB CCCC DDDD);",
]

update_file(caseDir + './system/decomposeParDict',idStr,pVals)

#-----------------------------------------------------------------------
# TRANSPORTPROPERTIES PREPARATION
#-----------------------------------------------------------------------
pVals   = [
    "nu             nu   [0 2 -1 0 0 0 0] %g;"%(nuF),
    "rho            rho  [1 -3 0 0 0 0 0] %g;"%(rhoF),
]  
idStr   = [
    "nu             nu   [0 2 -1 0 0 0 0] YYYY;",
    "rho            rho  [1 -3 0 0 0 0 0] XXXX;",
]

update_file(caseDir + './constant/transportProperties',idStr,pVals)

#-----------------------------------------------------------------------
# HFDIBDEMDICT PREPARATION
#-----------------------------------------------------------------------
pVals   = [
    "nDirectForcingIterations        %d;"%(nDirectForcingIterations),
    "immersedUTolerance              %g;"%(immersedUTolerance),
    "rho         rho [1 -3 0 0 0 0 0] %g;"%(rhoS),
]  
idStr   = [
    "nDirectForcingIterations        AAAA;",
    "immersedUTolerance              BBBB;",
    "rho         rho [1 -3 0 0 0 0 0] AAAA;",
]

update_file(caseDir + './constant/HFDIBDEMDict',idStr,pVals)


