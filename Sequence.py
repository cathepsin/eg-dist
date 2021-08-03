class ProteinSequence:
    #Extract all ATOM lines from a .pdb (or .ent) file. Data is stored as Atoms, AminoAcids
    class Atom:#Class to store information about an individual atom
        def __init__(self, num, coord, tag, res, ch, occ, bf, ele):
            self.num = num
            self.location = coord
            self.id = tag
            self.occupancy = occ
            self.bfactor = bf
            self.element = ele
            self.residue = res
            self.chain = ch

        def __repr__(self):
            return f"ATOM: NUMBER: {self.num}, TAG: {self.id}, RESIDUE: {self.residue}, CHAIN: {self.chain}, COORDINATES: {self.location}, OCCUPANCY: {self.occupancy}, B_FACTOR: {self.bfactor}, ELEMENT: {self.element}"
        def __str__(self):
            return f"{self.id} {self.residue} {self.location}"

    class AminoAcid:#Class to store information about a whole residue
        def __init__(self, numb, atomList, cName, ch):
            self.num = numb
            self.atoms = atomList
            self.residue = cName
            self.chain = ch
            self.rotation = dict()
            for atom in self.atoms:
                print("Residue: ", atom.residue, " and ", self.residue)
                if atom.residue != self.residue:
                    if not atom.residue in self.rotation:
                        self.rotation[atom.residue] = list()
                    self.rotation[atom.residue].append(atom)
            for key in self.rotation:
                for atom in self.rotation[key]:
                    self.atoms.remove(atom)
            self.specialAtoms = list()
            for atom in self.atoms:
                if atom.id.upper().find('X') != -1:
                    self.specialAtoms.append(atom)
            for atom in self.specialAtoms:
                self.atoms.remove(atom)

        def __lt__(self, other):
            return self.num < other.num

        def __repr__(self):
            return f"RESIDUE: {self.residue}, NUMBER: {self.num}, CHAIN: {self.chain}"

        def __str__(self):
            return f"{self.residue} {self.num} {self.chain}"

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

    #Parses a .pdb or .ent file. Returns an AA sequence.
    #Each residue is placed in the order they appear in the .pdb/.ent file
    def parsePDB(self, file):
        currResNum = -1
        atomGroup = []
        for line in file:
            if line.find("ATOM") == 0:
                spl = self.strToList(line)
                resNum = float(line[22:27])
                if currResNum == -1:
                    #Not all files start at residue 0 or 1. This ensures correct starting position
                    currResNum = resNum

                #Check if the current ATOM is on a new residue
                if resNum != currResNum:
                    #Found the next residue, so save previously recorded data as an AminoAcid
                    self.sequence.append(self.AminoAcid(currResNum, atomGroup.copy(), atomGroup[0].residue, chain))
                    currResNum = resNum
                    atomGroup.clear()
                #Parse ATOM data
                num = float(line[4:11])
                tag = line[11:16].strip()
                residue = line[16:20].strip()
                chain = line[20:22].strip()
                coordinates = [float(line[27:38]),float(line[38:46]),float(line[46:55])]
                occupancy = float(line[55:60])
                bfactor = float(line[60:73])
                element = line[73:].strip()
                #Save atom
                newAtom = self.Atom(num, coordinates, tag, residue, chain, occupancy, bfactor, element)
                atomGroup.append(newAtom)
                #print(newAtom)
        #Save the final AminoAcid instance
        self.sequence.append(self.AminoAcid(currResNum, atomGroup.copy(), atomGroup[0].residue, chain))


#ATOM   1895  CB  ILE B 954      32.149   7.481 198.133  1.00 89.00           C
#ATOM      1  N   VAL A   1       6.204  16.869   4.854  1.00 49.05           N
#0-4| 4-11 |11-16|16-20|20-22|22-27|27-38|38-46|46-55|  55-60|  60-73| 73-end