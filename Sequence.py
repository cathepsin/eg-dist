import json

class ProteinSequence:
    #Store data about a protein from a .pdb and a .socket
    class Atom:
        def __init__(self, num, coord, tag, ele):
            self.num = num
            self.location = coord
            self.id = tag
            self.element = ele

    class AminoAcid:
        def __init__(self, numb, atomList, cName, ch):
            self.num = numb
            self.atoms = atomList
            self.code = cName
            self.chain = ch

    def __init__(self):
        self.sequence = []
        print("Making this guy!")

    def cleanLargeTag(self, str):
        #If an atomID is 4+ characters (as in ATOM#93 for PDB: 5u59), str.split() will not work as intended
        lst = str.split()



    def parsePDB(self, file):
        currNum = -1
        currName = ""
        atomGroup = []
        for line in file:
            if line.find("ATOM") == 0:
                spl = line.split()
                #Parse atomic data
                num = spl[1]
                if num == '93':
                    print("stop here")
                tag = spl[2]
                resNum = spl[5]
                if currNum == -1:
                    currNum = resNum
                if resNum != currNum:
                    #Found the next residue, so save previously recorded data as an AminoAcid instance
                    self.sequence.append(self.AminoAcid(currNum, atomGroup.copy(), residue, chain))
                    currNum = resNum
                    #Wipe stored data about the previously recorded AminoAcid and start over
                    atomGroup.clear()
                chain = spl[4]
                residue = spl[3]
                coordinates = [spl[6],spl[7],spl[8]]
                element = spl[11]
                #Save atom
                newAtom = self.Atom(num, coordinates, tag, element)
                atomGroup.append(newAtom)
        #Save the final AminoAcid instance
        self.sequence.append(self.AminoAcid(currNum, atomGroup.copy(), residue, chain))
        print("wow")