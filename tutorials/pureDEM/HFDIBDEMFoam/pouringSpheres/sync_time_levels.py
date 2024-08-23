#!usr/bin/python3

import os

def canBeConvertedToFloat(input):
    try:
        float(input)
        return True
    except ValueError:
        return False
    
def shouldBeConvertedInt(input):
    try:
        if(float(input) == int(input)):
            return True
        else:
            return False
    except ValueError:
        return False

def getAllTimeLevelDirs():
    allDirs = os.listdir()
    timeLevelDirs = []
    for dir in allDirs:
        if(canBeConvertedToFloat(dir)):
            if(shouldBeConvertedInt(dir)):
                timeLevelDirs.append(int(dir))
            else:
                timeLevelDirs.append(float(dir))
    timeLevelDirs.sort(key = float)
    return timeLevelDirs

# Rename all time level directories to integers
if(os.path.isdir('renamedDirs')):
    os.system("rm -r renamedDirs")
os.system("mkdir renamedDirs")
timeLevelDirs = getAllTimeLevelDirs()
print("This is sorted time level directories: ", timeLevelDirs)
print("Do you want to rename these directories to integers? (Y/n)")
answer = input()
if(answer == 'n' or answer == 'N'):
    print("Exiting...")
    exit()

i = 0 # adjust this to start from a different number
for dir in timeLevelDirs:
    os.system("mv " + str(dir) + " renamedDirs/" + str(i))
    i += 1
os.system("mv renamedDirs/* .")
os.system("rmdir renamedDirs")
