import numpy as np

class Chain:
    def __init__(self, sequence):
        self.chains = self.MakeChains(sequence)
        self.message = ""

    #Scan sequence and organize pointers to each residue into grouped chains
    def MakeChains(self, sequence):
        chains = dict()
        for residue in sequence.sequence:
            if not residue.chain in chains:
                chains[residue.chain] = dict()
            chains[residue.chain][residue.num] = residue
        return chains

    #Some .pdb files contain a translation matrix to generate symmetry pairs
    #Determine the translation matrix from REMARK 350 and store
    #TODO Check .mmol files to determine if this is even needed
    def Symmetry(self, infile, sequence):
        foundRemark = False
        mat = dict()
        copyChains = list()
        for line in infile:
            if line.find('REMARK') == 0:
                foundRemark = True
            elif foundRemark == True:
                break
            if line.find('REMARK 350 APPLY THE FOLLOWING TO CHAINS:') != -1:
                spl = line.split()
                for ch in reversed(spl):
                    if ch.find(':') != -1:
                        break
                    copyChains.append(ch)
            if line.find('REMARK 350   BIOMT1') != -1:
                #Parse matrix
                while line.find("BIOMT") != -1:
                    spl = line.split()
                    row1 = [float(spl[4]), float(spl[5]), float(spl[6]), float(spl[7])]
                    spl = next(infile).split()
                    row2 = [float(spl[4]), float(spl[5]), float(spl[6]), float(spl[7])]
                    spl = next(infile).split()
                    row3 = [float(spl[4]), float(spl[5]), float(spl[6]), float(spl[7])]
                    row4 = [0, 0, 0, 1] #4th row of the identity matrix
                    if not np.all(np.equal(np.array([row1, row2, row3, row4]), np.identity(4))):
                        mat[spl[3]] = np.array([row1, row2, row3, row4])
                    line = next(infile)
        self.matrix = mat
        self.copyChains = copyChains