#!/usr/bin/env python

'''
Created on Sep 18, 2013

@author: kirill lykov
'''
import os
import argparse
from Atom2objTranslator import Atom2objTranslator
from splitObj import SplitObj

# There is two options to be used for the restart to data file convertation:
# 1. restart2data utility
# 2. ./lammps -r data.restart remap data.atom
programType = "restart2data"

# restart file name is like <core name>.restart.<iteration number>
# the method returns  <core name>.<iteration number>
def getCoreName(fileName):
    ind = fileName.find('.')
    itNumInd = fileName.rfind('.')
    assert(ind != -1 and itNumInd != -1)
    return fileName[:ind] + fileName[itNumInd:]
    
def getAtomName(fileName):
    res = getCoreName(fileName) + ".atom"
    return res

def getObjName(fileName, index):
    ind = fileName.find('.')
    assert(ind != -1)
    res = fileName[:ind + 1]
    res.replace('_', '-')
    return res + str(index) + ".obj"

def getIteration(fileName):
    itNumInd = fileName.rfind('.')
    assert(itNumInd != -1)
    return fileName[itNumInd+1:]

######## Main ########
print("restart2obj started")
parser = argparse.ArgumentParser(description='Transforms all restart files which are in the current directory to obj files.', 
                                 usage= './restart2obj.py')
# to be implemented
#parser.add_argument('-i','--indexStart', help='Outfiles first index', required=False)
parser.add_argument('-c','--cut', help='cut obj or not', required=False, default=False)
args = vars(parser.parse_args())

files = [f for f in os.listdir('.') if (os.path.isfile(f) and f.find("restart") != -1)]

filesSorted = list(tuple())
for f in files:
    it = int(getIteration(f))
    filesSorted.append((it, f))

filesSorted.sort(key=lambda x: x[0])
print("there are " + str(len(files)) + " restart files found")

# colors specified in mlt file, assumed that mlt file is called the <coreName>.mlt
# and one for all obj files
colors = ("red", "blue", "blue")

translator = Atom2objTranslator(0.1, 2, 2, colors)
for f,ind in zip(filesSorted, range(0,len(filesSorted))):
    atomFileName = getAtomName(f[1])
    if programType == "restart2data ":
        cmd = "restart2data " + f[1] + " " + atomFileName
    else:
        cmd = "./lammps -r " + f[1] + " remap " + atomFileName
    
    print(atomFileName + ": restart2atom")
    os.system(cmd)
    
    print(atomFileName + ": atom2obj")
    translator.run(atomFileName, getObjName(f[1], ind))
    
    os.remove(atomFileName)
