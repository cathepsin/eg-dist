import Sequence

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
