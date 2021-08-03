import numpy as np

class Chain:
    def __init__(self, sequence):
        self.chains = self.MakeChains(sequence)

    def MakeChains(self, sequence):
        chains = dict()
        for residue in sequence.sequence:
            if not residue.chain in chains:
                chains[residue.chain] = list()
            chains[residue.chain].append(residue)
        for key in chains:
            chains[key].sort()
        return chains

    #Some .pdb files contain a translation matrix to generate a full coiled-coil
    def Symmetry(self, infile, sequence):
        foundRemark = False
        for line in infile:
            if line.find('REMARK') == 0:
                foundRemark = True
                mat = dict()
            elif foundRemark == True:
                break
            if line.find('REMARK 350   BIOMT1') != -1:
                while line.find("BIOMT") != -1:
                    spl = line.split()
                    row1 = [float(spl[4]),float(spl[5]),float(spl[6]),float(spl[7])]
                    spl = next(infile).split()
                    row2 = [float(spl[4]), float(spl[5]), float(spl[6]), float(spl[7])]
                    spl = next(infile).split()
                    row3 = [float(spl[4]), float(spl[5]), float(spl[6]), float(spl[7])]
                    mat[spl[3]] = np.array([row1, row2, row3])
                    line = next(infile)

                print("here")






        print("here we go...")
