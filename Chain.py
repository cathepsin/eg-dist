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

    # Some .pdb files contain a translation matrix to generate a full coiled-coil
    def Symmetry(self, infile, sequence):
        foundRemark = False
        for line in infile:
            if line.find('REMARK') == 0:
                foundRemark = True
                mat = dict()
            elif foundRemark == True:
                break
            copyChains = list()
            if line.find('REMARK 350 APPLY THE FOLLOWING TO CHAINS:') != -1:
                spl = line.split()
                for ch in reversed(spl):
                    if ch.find(':') != -1:
                        break
                    copyChains.append(ch)
            if line.find('REMARK 350   BIOMT1') != -1:
                while line.find("BIOMT") != -1:
                    spl = line.split()
                    row1 = [float(spl[4]), float(spl[5]), float(spl[6]), float(spl[7])]
                    spl = next(infile).split()
                    row2 = [float(spl[4]), float(spl[5]), float(spl[6]), float(spl[7])]
                    spl = next(infile).split()
                    row3 = [float(spl[4]), float(spl[5]), float(spl[6]), float(spl[7])]
                    row4 = [0, 0, 0, 1] #4th row of the identity matrix
                    mat[spl[3]] = np.array([row1, row2, row3, row4])
                    line = next(infile)
        self.matrix = mat
        self.copyChains = copyChains
