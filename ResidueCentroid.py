import collections
class CentroidFinder:
    def __init__(self):
        # Using the .pdb/.ent tag encodings, all atoms besides alpha/beta carbons are considered when finding a centroid
        # *For Glycine, only the alpha carbon is considered
        # *For Alanine, only the beta carbon is considered
        self.AAs = {
            "ARG": ['CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
            "HIS": ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
            "LYS": ['CG', 'CD', 'CE', 'NZ'],
            "ASP": ['CG', 'OD1', 'OD2'],
            "GLU": ['CG', 'CD', 'OE1', 'OE2'],
            "SER": ['OG'],
            "THR": ['OG1', 'CG2'],
            "ASN": ['CG', 'OD1', 'ND2'],
            "GLN": ['CG', 'CD', 'OE1', 'NE2'],
            "CYS": ['SG'],
            "SEC": ['SE'],
            "GLY": ['CA'],
            "PRO": ['CG', 'CD'],
            "ALA": ['CB'],
            "VAL": ['CG1', 'CG2'],
            "ILE": ['CG1', 'CG2', 'CD1'],
            "LEU": ['CG', 'CD1', 'CD2'],
            "MET": ['CG', 'SD', 'CE'],
            "PHE": ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            "TYR": ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
            "TRP": ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
        }
        self.stored = []

    def GetCentroid(self, residue):
        x = 0
        y = 0
        z = 0
        n = 0
        # Because, unfortunately, some pdb files have incorrectly labled residues, we must check (as in residue 1340 in 1PL5)
        #This causes the method to be relatively slow, but ensures correctness
        if residue.num == 428:
            print("stopping")
        resName = self.CheckResidue(residue)
        if resName != residue.residue and resName != "Unknown residue":
            print("Uh oh! Somebody made a bad file! ", residue.residue, "-->", resName , "(", residue.num,")")
        elif resName == "Unknown residue":
            print("Unknown residue")
            #TODO Add another check to get the correct residue. Until then, return -1
            return residue, [-1,-1,-1]
        for atom in residue.atoms:
            if atom.id in self.AAs[resName]:
                n = n + 1
                x = x + atom.location[0]
                y = y + atom.location[1]
                z = z + atom.location[2]
        x = x / n
        y = y / n
        z = z / n

        return residue, [x, y, z]

    def CheckResidue(self, res):
        checkList = []
        for atom in res.atoms:
            checkList.append(atom.id)
        if collections.Counter(['N','CA','C','O']) == collections.Counter(checkList):
            return 'GLY'
        if collections.Counter(['N','CA','C','O','CB']) == collections.Counter(checkList):
            return 'ALA'
        for AA in self.AAs:
            if collections.Counter(self.AAs[AA] + ['N','CA','C','O','CB']) == collections.Counter(checkList):
                return AA
        return "Unknown residue"


