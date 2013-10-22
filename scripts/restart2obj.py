#!/usr/bin/env python
'''
(C) Copyright Kirill Lykov 2013.
Distributed under the GNU Software License (See accompanying file LICENSE)

Convert collection of restart files into object files, assumed that the restart2data 
is available
@author: kirill lykov
'''


import os
from src.Atom2objTranslator import Atom2objTranslator

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

def getObjName(fileName):
    res = getCoreName(fileName) + ".obj"
    return res 

files = [f for f in os.listdir('.') if (os.path.isfile(f) and f.find("restart") != -1)]
translator = Atom2objTranslator()
for f in files:
    atomFileName = getAtomName(f)
    cmd = "restart2data " + f + " " + atomFileName
    os.system(cmd)
    translator.run(atomFileName, getObjName(f))
    os.remove(atomFileName)