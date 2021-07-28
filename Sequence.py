class ProteinSequence:
    #Extract all ATOM lines from a .pdb (or .ent) file. Data is stored as Atoms, AminoAcids
    class Atom:#Class to store information about an individual atom
        def __init__(self, num, coord, tag, occ, bf, ele):
            self.num = num
            self.location = coord
            self.id = tag
            self.occupancy = occ
            self.bfactor = bf
            self.element = ele

    class AminoAcid:#Class to store information about a whole residue
        def __init__(self, numb, atomList, cName, ch):
            self.num = numb
            self.atoms = atomList
            self.residue = cName
            self.chain = ch

        def __lt__(self, other):
            return self.num < other.num

    def __init__(self):
        self.sequence = []

    #Split an ATOM line into a list
    def strToList(self, str):
        retLst = str.split()
        if len(retLst[2]) > 4:
            #If an ATOM tag is 4+ characters, str.split() will not work as intended
            retLst = self.cleanLargeTag(str)
        if len(retLst[9]) > 4:
            #If an ATOM bfactor is 4+ characters, str.spli() will not work as intended
            retLst = self.cleanLargeBFactor(str)
        return retLst

    #Takes the original ATOM line (string) and splits the tag and residue correctly
    def cleanLargeTag(self, str):
        lst = str.split()
        retlst = lst.copy()
        retlst.insert(2,lst[2][:4])
        retlst.insert(3,lst[2][4:])
        retlst.remove(retlst[4])
        return retlst

    #Takes the original ATOM line (string) and splits the occupancy and bfactor correctly
    def cleanLargeBFactor(self, str):
        lst = str.split()
        retlst = lst.copy()
        retlst.insert(9,lst[9][:4])
        retlst.insert(10,lst[9][4:])
        retlst.remove(retlst[11])
        return retlst

    # def GetAtomDist(self, atom1, atom2):
    #     x_term = (atom1.location[0] - atom2.location[0]) ** 2
    #     y_term = (atom1.location[1] - atom2.location[1]) ** 2
    #     z_term = (atom1.location[2] - atom2.location[2]) ** 2
    #     return math.sqrt(x_term + y_term + z_term)

    # def SplitChains(self):
    #     chainset = {}
    #     setLength = len(chainset)
    #     for val in self.sequence:
    #         prevLen = len(chainset)
    #         chainset.add(val.chain)

    #Parses a .pdb or .ent file. Returns an AA sequence.
    #Each residue is placed in the order they appear in the .pdb/.ent file
    def parsePDB(self, file):
        currResNum = -1
        currName = ""
        atomGroup = []
        for line in file:
            if line.find("ATOM") == 0:
                spl = self.strToList(line)
                resNum = spl[5]
                if currResNum == -1:
                    #Not all files start at residue 0 or 1. This ensures correct starting position
                    currResNum = resNum

                #Check if the current ATOM is on a new residue
                if resNum != currResNum:
                    #Found the next residue, so save previously recorded data as an AminoAcid
                    self.sequence.append(self.AminoAcid(currResNum, atomGroup.copy(), residue, chain))
                    currResNum = resNum
                    atomGroup.clear()
                #Parse ATOM data
                num = spl[1]
                tag = spl[2]
                residue = spl[3]
                chain = spl[4]
                coordinates = [float(spl[6]),float(spl[7]),float(spl[8])]
                occupancy = spl[9]
                bfactor = spl[10]
                element = spl[11]
                #Save atom
                newAtom = self.Atom(num, coordinates, tag, occupancy, bfactor, element)
                atomGroup.append(newAtom)
        #Save the final AminoAcid instance
        self.sequence.append(self.AminoAcid(currResNum, atomGroup.copy(), residue, chain))