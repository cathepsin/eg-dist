import collections


class CentroidFinder:
    def __init__(self):
        # Using the .pdb/.ent tag encodings, all atoms besides alpha/beta carbons are considered when finding a centroid
        # *For Glycine, only the alpha carbon is considered
        # *For Alanine, only the beta carbon is considered
        self.AAs = {
            "SER": ['OG'],
            "CYS": ['SG'],
            "SEC": ['SE'],
            "GLY": ['CA'],
            "ALA": ['CB'],
            "THR": ['OG1', 'CG2'],
            "PRO": ['CG', 'CD'],
            "VAL": ['CG1', 'CG2'],
            "ASP": ['CG', 'OD1', 'OD2'],
            "ASN": ['CG', 'OD1', 'ND2'],
            "ILE": ['CG1', 'CG2', 'CD1'],
            "LEU": ['CG', 'CD1', 'CD2'],
            "MET": ['CG', 'SD', 'CE'],
            "LYS": ['CG', 'CD', 'CE', 'NZ'],
            "GLU": ['CG', 'CD', 'OE1', 'OE2'],
            "GLN": ['CG', 'CD', 'OE1', 'NE2'],
            "HIS": ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
            "ARG": ['CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
            "PHE": ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            "TYR": ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
            "TRP": ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
        }
        self.warning = ""

    def GetCentroid(self, residue):
        if residue.rotation:
            #TODO deal with rotamers
            self.warning += f"Residue {residue} has a decimal occupancy. Its centroid cannot be determined"
            return float('NaN')
        x = 0
        y = 0
        z = 0
        n = 0
        # Because, unfortunately, some pdb files have incorrectly labled residues, we must check (as in residue 1340 in 1PL5)
        # This causes the method to be relatively slow, but ensures correctness
        resName = self.CheckResidue(residue)
        if resName != residue.residue and resName != "Unknown residue":
            print("Uh oh! Somebody made a bad file! ", residue.residue, "-->", resName, "(", residue.num, ")")
        elif resName == "Unknown residue":
            #Either a mistake in the protein file or simplified.
            self.warning += f"Residue {residue} could not be identified from its ATOM listings. The centroid will be " \
                            f"determined using available ATOMS"
            print("Unknown residue")
            return self.UnknownCentroid(residue)
        for atom in residue.atoms:
            if atom.id in self.AAs[resName]:
                n = n + 1
                x = x + atom.location[0]
                y = y + atom.location[1]
                z = z + atom.location[2]
        return [x/n, y/n, z/n]

    #Unfortunately, some PDB files have incorrect data
    #Determine what residue res is from its ATOM data
    def CheckResidue(self, res):
        checkList, addList = self.GetLists(res)
        if collections.Counter(['N', 'CA', 'C', 'O']) == collections.Counter(checkList + addList):
            return 'GLY'
        if collections.Counter(['N', 'CA', 'C', 'O', 'CB']) == collections.Counter(checkList + addList):
            return 'ALA'
        for AA in self.AAs:
            if collections.Counter(self.AAs[AA] + ['N', 'CA', 'C', 'O', 'CB']) == collections.Counter(checkList + addList):
                return AA
        return "Unknown residue"

    #Get all atoms that are not hydrogen. If possible, don't include CA or CB
    def UnknownCentroid(self, res):
        checkList, addList = self.GetLists(res)
        print("here")
        x = 0
        y = 0
        z = 0
        n = 0
        for atom in res.atoms:
            if atom.id in checkList + addList and atom.id not in ['N', 'CA', 'C', 'O', 'CB']:
                #Use available ATOM information to make a best-informed centroid
                n += 1
                x += atom.location[0]
                y += atom.location[1]
                z += atom.location[2]
        try:
            return [x/n, y/n, z/n]
        except ZeroDivisionError:
            #In case of glycine or alanine
            for atom in res.atoms:
                if atom.id == 'CA':
                    ca = atom
                elif atom.id == 'CB':
                    return atom.location
            return ca.location

    def GetLists(self, res):
        checkList = []
        addList = []
        if len(res.rotation.keys()) > 0:
            keyList = list(res.rotation.keys())
            for atom in [atom for atom in res.rotation[keyList[0]] if "H" not in atom.element]:
                addList.append(atom.id)
        for atom in [atom for atom in res.atoms if "H" not in atom.element]:
            checkList.append(atom.id)
        return checkList, addList


