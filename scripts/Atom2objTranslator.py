'''
Created on Sep 18, 2013

@author: kirill lykov
'''
from itertools import dropwhile
from collections import OrderedDict

class Atom2objTranslator:
    
    # max area and edge length used to skip triangles which are separated by pbc
    scaleFactor = 1.0
    maxEdgeLegth2 = 2.0
    maxTriangleArea2 = 2.0
    colors = None
    
    def __init__(self, scaleFactor = None, maxEdgeLegth2 = None, maxTriangleArea2 = None, colors=None):
        if scaleFactor is not None:
            self.scaleFactor = scaleFactor
        if maxEdgeLegth2 is not None:
            self.maxEdgeLegth2 = maxEdgeLegth2
        if maxTriangleArea2 is not None:
            self.maxTriangleArea2 = maxTriangleArea2
        self.colors = colors
            
    def is_not_atom_section(self, s):
        ss = s.rstrip().split() + [""]
        return (ss[0] != "Atoms")
    
    def is_not_angle_section(self, s):
        ss = s.rstrip().split() + [""] # just to avoid check on size
        return (ss[0] != "Angles")
    
    def getAtomIndex(self, s):
        words = s.split()
        return int(words[0])
    
    def getAtomType(self, s):
        words = s.split()
        return int(words[2])
    
    def atom2vertex(self, s):
        words = s.split()
        if self.scaleFactor == 1:
            res = ' '.join(['v', words[3], words[4], words[5]])
        else:
            scaledVals = (self.scaleFactor * float(words[3]), \
                          self.scaleFactor * float(words[4]), \
                          self.scaleFactor * float(words[5]))
            res = ' '.join(['v', str(scaledVals[0]), str(scaledVals[1]), str(scaledVals[2])])
        return res
    
    def angle2face(self, s, old2new):
        words = s.split()
        v = []
        for i in range(2, 5):
            v.insert( i - 1, str(old2new[int(words[i])]) )
        res = ' '.join(['f', v[0], v[1], v[2]])
        return res
    
    def checkAngle(self, vertices, line):
        words = line.split()
        for i in range(2, 5):
            v = int(words[i])
            if (v not in vertices):
                return False
        return True
    
    def angleToVertices(self, line):
        words = line.split()
        res = list()
        for i in range(2, 5):
            v = int(words[i])
            res.append(v)
        return res
    
    def getPointFromLine(self, line):
        words = line.split()
        res = list()
        for i in range(3, 6):
            res.append(float(words[i]))
        return res
    
    #remove angles with too long edges or big area
    def angleIsBig(self, atomsDict, line):
        words = line.split()

        pointA = self.getPointFromLine(atomsDict[ int(words[2]) ])
        pointB = self.getPointFromLine(atomsDict[ int(words[3]) ])
        pointC = self.getPointFromLine(atomsDict[ int(words[4]) ])
        
        bma = list()
        lengthBma2 = 0.0
        cma = list()
        lengthCma2 = 0.0
        for i in range(0, 3):
            bma.append(pointB[i] - pointA[i])
            lengthBma2 = lengthBma2 + (pointB[i] - pointA[i])*(pointB[i] - pointA[i])
            cma.append(pointC[i] - pointA[i])
            lengthCma2 = lengthCma2 + (pointC[i] - pointA[i])*(pointC[i] - pointA[i])
        
        area2 = (bma[1]*cma[2] - bma[2]*cma[1])*(bma[1]*cma[2] - bma[2]*cma[1]) + \
        (bma[2]*cma[0] - bma[0]*cma[2])*(bma[2]*cma[0] - bma[0]*cma[2]) + \
        (cma[0]*bma[1] - bma[0]*cma[1])*(cma[0]*bma[1] - bma[0]*cma[1])
        
        if lengthBma2 > self.maxEdgeLegth2 or lengthCma2 > self.maxEdgeLegth2 or area2 > self.maxTriangleArea2:
            return True
        return False
    
    # get sorted lines from file according to the fist number (required by obj format)
    def getSortedAtomsLines(self, lines):
        atomsDict = {}
        alreadyParsed = False
       
        for line in lines:
            line = line.rstrip()
            line = line.rstrip()
            # skip Atoms and the first line after Atoms section
            if ((line.split() + [""])[0] == "Atoms"): continue
            if ((not line) and (not alreadyParsed)):
                alreadyParsed = True
                continue
            # break when came to the end of the section
            if (not line) and (alreadyParsed):
                break
            words = line.split()
            assert(len(words) > 0)
            n = int(words[0])
            atomsDict[n] = line;
        
        return OrderedDict(sorted(atomsDict.items(), key=lambda t: t[0]))
    
    #assign a color according to the type of atom. Assumed that one triangle has atoms of one type
    def getColorLine(self, line, atomInd2type):
        words = line.split()
        typeOfAtoms = (atomInd2type[int(words[1])] - 1) % len(self.colors)
        return "usemtl " + self.colors[typeOfAtoms] + "\n"
        
    def run(self, inputFileName, outputFileName):
        try:
            verticesFromAngles = set() # to avoid adding point which has no angles
            with open(inputFileName, 'r') as f:
                # parse Angles
                alreadyParsed = False    
                lines = dropwhile(self.is_not_angle_section, f)
                for line in lines:
                    line = line.rstrip()
                    # skip Atoms and the first line after Angles section
                    if ((line.split() + [""])[0] == "Angles"): continue
                    if ((not line) and (not alreadyParsed)):
                        alreadyParsed = True
                        continue
                    # break when came to the end of the section
                    if (not line) and (alreadyParsed):
                        break
                    
                    verticesFromAngles.update(self.angleToVertices(line))
            
            output = open(outputFileName, "w")
            output.write("# Generated by atom2obj script\n\n")
            if self.colors is not None:
                output.write("  mtllib " + inputFileName[:inputFileName.find('.')]+ ".mtl\n\n")
            
            with open(inputFileName, 'r') as f:
                
                # set of vertices to check that all vertices mentioned in Angles are in Atoms
                vertices = set()
                
                lines = dropwhile(self.is_not_atom_section, f)
                
                atomsDict = self.getSortedAtomsLines(lines)
                old2new = dict()
                atomInd2type = dict()
                # parse Atoms and write
                output.write("\n# Vertex list\n\n")
                newLineInd = 1
                for ind in atomsDict:
                    line = atomsDict[ind]
                    
                    vertices.add(self.getAtomIndex(line))
                    old2new[self.getAtomIndex(line)] = newLineInd
                    atomInd2type[newLineInd] = self.getAtomType(line)
                    
                    outLine = self.atom2vertex(line)
                    if ind in verticesFromAngles:
                        output.write(outLine + "\n")
                        newLineInd = newLineInd + 1
                
                # parse Angles
                output.write("\n# Triangle list\n\n")
                alreadyParsed = False    
                lines = dropwhile(self.is_not_angle_section, lines)
                for line in lines:
                    line = line.rstrip()
                    # skip Atoms and the first line after Angles section
                    if ((line.split() + [""])[0] == "Angles"): continue
                    if ((not line) and (not alreadyParsed)):
                        alreadyParsed = True
                        continue
                    # break when came to the end of the section
                    if (not line) and (alreadyParsed):
                        break
                    
                    if (not self.checkAngle(vertices, line)):
                        raise RuntimeError("There are vertices from Angles section which are not in Atoms section. Line content is \"" + line + "\"")
                    if (not self.angleIsBig(atomsDict, line)):
                        outLine = self.angle2face(line, old2new)
                        output.write(outLine + "\n")
                        if self.colors is not None:
                            outLineColor = self.getColorLine(outLine, atomInd2type)
                            output.write(outLineColor + "\n")
                    
        except IOError:
            print("Error: can\'t find file or read data")
        except RuntimeError as e:
            print(e.args)
