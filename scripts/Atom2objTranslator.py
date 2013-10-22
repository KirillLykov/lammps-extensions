'''
(C) Copyright Kirill Lykov 2013.
Distributed under the GNU Software License (See accompanying file LICENSE)

@author: kirill lykov
'''
from itertools import dropwhile
from collections import OrderedDict

class Atom2objTranslator:
    
    def is_not_atom_section(self, s):
        return (s.rstrip() != "Atoms")
    
    def is_not_angle_section(self, s):
        return (s.rstrip() != "Angles")
    
    def getAtomIndex(self, s):
        words = s.split()
        return int(words[0])
    
    def atom2vertex(self, s):
        words = s.split()
        res = ' '.join(['v', words[3], words[4], words[5]])
        return res
    
    def angle2face(self, s):
        words = s.split()
        res = ' '.join(['f', words[2], words[3], words[4]])
        return res
    
    def checkAngle(self, vertices, line):
        words = line.split()
        for i in range(1, 3):
            v = int(words[i])
            if (v not in vertices):
                return False
        return True
    
    # get sorted lines from file according to the fist number (required by obj format)
    def getSortedAtomsLines(self, lines):
        atomsDict = {}
        alreadyParsed = False
       
        for line in lines:
            line = line.rstrip()
            line = line.rstrip()
            # skip Atoms and the first line after Atoms section
            if (line == "Atoms"): continue
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
    
    def run(self, inputFileName, outputFileName):
        try:
            output = open(outputFileName, "w")
            output.write("# Generated by atom2obj script\n\n")
            
            with open(inputFileName, 'r') as f:
                
                # set of vertices to check that all vertices mentioned in Angles are in Atoms
                vertices = set()
                
                lines = dropwhile(self.is_not_atom_section, f)
                
                atomsDict = self.getSortedAtomsLines(lines)
                
                # parse Atoms and write
                output.write("\n# Vertex list\n\n")
                for ind in atomsDict:
                    line = atomsDict[ind]
                    
                    vertices.add(self.getAtomIndex(line) )
                    outLine = self.atom2vertex(line)
                    output.write(outLine + "\n")
                
                # parse Angles
                output.write("\n# Triangle list\n\n")
                alreadyParsed = False    
                lines = dropwhile(self.is_not_angle_section, lines)
                for line in lines:
                    line = line.rstrip()
                    # skip Atoms and the first line after Angles section
                    if (line == "Angles"): continue
                    if ((not line) and (not alreadyParsed)):
                        alreadyParsed = True
                        continue
                    # break when came to the end of the section
                    if (not line) and (alreadyParsed):
                        break
                    
                    if (not self.checkAngle(vertices, line)):
                        raise RuntimeError("There are vertices from Angles section which are not in Atoms section. Line content is \"" + line + "\"")
                    outLine = self.angle2face(line)
                    output.write(outLine + "\n")
                    
        except IOError:
            print("Error: can\'t find file or read data")
        except RuntimeError as e:
            print(e.args)
        finally:            
            f.close()
            output.close()
