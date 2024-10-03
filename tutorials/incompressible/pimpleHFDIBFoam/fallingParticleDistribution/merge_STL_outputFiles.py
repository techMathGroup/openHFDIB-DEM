#!/usr/bin/python3
#ReadME: This script merges the stl files from the bodiesInfo/ directory into a single file in the STLMerged/ directory.
#created by OStudenik 
import os
import numpy # as np
def canBeConvertedToFloat(input):
    try:
        float(input)
        return True
    except ValueError:
        return False

def getParticles_List(caseDir):
    "Returns a list with the particle numbers in the case directory."
    Directory = caseDir
    Strings   = list(set([numStr for numStr in os.listdir(Directory)]))
    Full_List = []
    for i in range(len(Strings)):
        string_to_save = caseDir+'/'+Strings[i]
        if(Strings[i][-4:] == '.stl'):
                Full_List.append(Strings[i][:-4])
    if(len(Full_List) > 0):            
        Full_List.sort(key = float)
    return Full_List 

Full_List = os.listdir('bodiesInfo/')
Full_List.sort(key = float)

if(not (os.path.isdir('STLMerged'))):
    os.system('mkdir STLMerged')
    time_iter = 0

    for item in Full_List:
        Full_List_II= getParticles_List('bodiesInfo/'+item+'/stlFiles')
        time_iter += 1
        for item_II in Full_List_II:
            with open('bodiesInfo/'+item+'/stlFiles/'+item_II+'.stl', 'r') as file:
                data = file.readlines()    
            with open('STLMerged/STL_Results'+str(time_iter).zfill(4)+'.stl', 'a+') as file:
                file.writelines('solid '+item_II+'.stl'+'\n')
                for i in range(1, len(data) -1):
                    file.writelines(data[i])
                file.writelines('endsolid '+item_II+'\n')
        print("-- Reading: ", item)

elif(os.path.isdir('STLMerged')):
    reduced_list = os.listdir('STLMerged/')
    reduced_list2 = os.listdir('bodiesInfo/')
    for item in reduced_list2:
        if canBeConvertedToFloat(item):
            if(float(item) == int(float(item))):
                reduced_list2[reduced_list2.index(item)] = int(float(item))
            else:
                reduced_list2[reduced_list2.index(item)] = float(item)
    reduced_list2.sort(key = float)
    reduced_list2 = [str(item) for item in reduced_list2]             

    print(" -- Warning: STLMerged/ already exists. The program will append the data to the existing files size: ",len(reduced_list))  
    print(" -- Simulation length: ",len(reduced_list2))  
    time_iter = len(reduced_list)

    print(reduced_list2[len(reduced_list):])

    for item in reduced_list2[len(reduced_list):]:
        Full_List_II= getParticles_List('bodiesInfo/'+item+'/stlFiles')
        time_iter += 1
        for item_II in Full_List_II:
            with open('bodiesInfo/'+item+'/stlFiles/'+item_II+'.stl', 'r') as file:
                data = file.readlines()    
            with open('STLMerged/STL_Results'+str(time_iter).zfill(4)+'.stl', 'a+') as file:
                file.writelines('solid '+item_II+'.stl'+'\n')
                for i in range(1, len(data) -1):
                    file.writelines(data[i])
                file.writelines('endsolid '+item_II+'\n')
        print("-- Reading: ", item)
